#include <string>
#include <iostream>

#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include <AMReX_DataServices.H>
#include <WritePlotFile.H>
#include <AppendToPlotFile.H>

#include <AMReX_BLFort.H>

using namespace amrex;

extern "C" {
  void pushvtog(const int* lo,  const int* hi,
                const int* dlo, const int* dhi,
                Real* U, const int* Ulo, const int* Uhi,
                const int* nc);
  void gradient(const int* lo,  const int* hi,
                const Real* U, const int* Ulo, const int* Uhi,
                Real* G, const int* Glo, const int* Ghi,
                const Real* dx);
}

static
void 
print_usage (int,
             char* argv[])
{
  std::cerr << "usage:\n";
  std::cerr << argv[0] << " infile=<plotfilename> \n\tOptions:\n\tis_per=<L M N> gradVar=<name>\n";
  exit(1);
}

std::string
getFileRoot(const std::string& infile)
{
  vector<std::string> tokens = Tokenize(infile,std::string("/"));
  return tokens[tokens.size()-1];
}

int
main (int   argc,
      char* argv[])
{
  amrex::Initialize(argc,argv);
  {
  if (argc < 2)
    print_usage(argc,argv);

  ParmParse pp;

  if (pp.contains("help"))
    print_usage(argc,argv);

  std::string infile; pp.get("infile",infile);
  DataServices::SetBatchMode();
  Amrvis::FileType fileType(Amrvis::NEWPLT);
  DataServices dataServices(infile, fileType);
  if( ! dataServices.AmrDataOk()) {
    DataServices::Dispatch(DataServices::ExitRequest, NULL);
  }
  AmrData& amrData = dataServices.AmrDataRef();

  int finestLevel = amrData.FinestLevel();
  pp.query("finestLevel",finestLevel);
  int Nlev = finestLevel + 1;

  std::string gradVar = "temp"; pp.query("gradVar",gradVar);

  int idC = -1;
  const Vector<std::string>& plotVarNames = amrData.PlotVarNames();
  for (int i=0; i<plotVarNames.size(); ++i)
  {
    if (plotVarNames[i] == gradVar) idC = i;
  }
  if (idC<0) {
    Print() << "Cannot find required data in pltfile" << std::endl;
  }

  const int idCst = 0;
  const int nCompIn = idCst + 1;

  Vector<std::string> inVarNames(nCompIn);
  inVarNames[idCst] = plotVarNames[idC];
  Vector<int> destFillComps(nCompIn);
  for (int i=0; i<nCompIn; ++i) {
    destFillComps[i] = i;
  }

  const int idGr = nCompIn;
  const int nCompOut = idGr + BL_SPACEDIM +1 ; // 1 component stores the ||gradT||

  Vector<MultiFab*> state(Nlev);
  Vector<Geometry*> geoms(Nlev);
  const int nGrow = 1;

  RealBox rb(&(amrData.ProbLo()[0]),
             &(amrData.ProbHi()[0]));
  Vector<int> is_per(BL_SPACEDIM,1);
  pp.queryarr("is_per",is_per,0,BL_SPACEDIM);
  int coord = 0;

  for (int lev=0; lev<Nlev; ++lev) {

    const BoxArray ba = amrData.boxArray(lev);
    const DistributionMapping dmap(ba);
    state[lev] = new MultiFab(ba,dmap,nCompOut,nGrow);
    const Vector<Real>& delta = amrData.DxLevel()[lev];

    // Get input state data onto intersection ba
    const int myNComp = 1; // gonna need this for fortran calls
    
    //state[lev]->setBndry(0.0); // Initialize grow cells to 0....here don't care to get grads right on c-f

    Print() << "Doing FillVar at lev: " << lev << std::endl;
    amrData.FillVar(*state[lev],lev,inVarNames,destFillComps);
    geoms[lev] = new Geometry(amrData.ProbDomain()[lev],&rb,coord,&(is_per[0]));

    // Fix up grow cells.  Use extrap for guess
    const Box& dbox = amrData.ProbDomain()[lev];
    for (MFIter mfi(*state[lev]); mfi.isValid(); ++mfi)
    {
      FArrayBox& fab = (*state[lev])[mfi];
      const Box& box = mfi.validbox();
      pushvtog(BL_TO_FORTRAN_BOX(box),
               BL_TO_FORTRAN_BOX(dbox),
               BL_TO_FORTRAN_ANYD(fab),
               &myNComp);
    }
    
    // Fix up fine-fine and periodic
    state[lev]->FillBoundary(geoms[lev]->periodicity());

    for (MFIter mfi(*state[lev]); mfi.isValid(); ++mfi)
    {
      FArrayBox& fab = (*state[lev])[mfi];
      const Box& box = mfi.validbox();
      gradient(BL_TO_FORTRAN_BOX(box),
               BL_TO_FORTRAN_N_ANYD(fab,0),
               BL_TO_FORTRAN_N_ANYD(fab,idGr),
               &(delta[0]));
    }
  }

  Vector<std::string> nnames(nCompOut);
  for (int i=0; i<nCompIn; ++i) {
    nnames[i] = inVarNames[i];
  }
  nnames[idGr+0] = gradVar + "_gx";
  nnames[idGr+1] = gradVar + "_gy";
#if BL_SPACEDIM==3
  nnames[idGr+2] = gradVar + "_gz";
#endif
  nnames[idGr+BL_SPACEDIM] = "||grad"+ gradVar+ "||";
  bool verb=false;
  bool appendPlotFile=false; pp.query("appendPlotFile",appendPlotFile);

  if (appendPlotFile) {
    int nStateOut = nCompOut - nCompIn;
    Vector<MultiFab*> ostate(Nlev);
    for (int lev=0; lev<Nlev; ++lev) {
      const BoxArray ba = state[lev]->boxArray();
      const DistributionMapping dm(ba);
      ostate[lev] = new MultiFab(ba,dm,nStateOut,nGrow);
      MultiFab::Copy(*ostate[lev],*state[lev],nCompIn,0,nStateOut,0);
    }
    Vector<std::string> namesOut(nStateOut);
    for (int i=0; i<nStateOut; ++i) {
      namesOut[i] = nnames[nCompIn+i];
    }

    std::string newMFBaseName = "NEWDAT"; pp.query("newMFBaseName",newMFBaseName);
    std::string newHeaderName = "NewHeader"; pp.query("newHeaderName",newHeaderName);
    AppendToPlotFile(amrData,ostate,infile,namesOut,newMFBaseName,newHeaderName,verb);

    Print() << "...finished.  Note: to see new data, you must rename NewHeader in the" << std::endl;
    Print() << "              pltfile to Header (probably want to save the original first)" << std::endl;
  } else {
    std::string outfile(getFileRoot(infile) + "_gt"); pp.query("outfile",outfile);

    Print() << "Writing new data to " << outfile << std::endl;

    WritePlotFile(state,amrData,outfile,verb,nnames);
  }
  }
  amrex::Finalize();
  return 0;
}
