#include <string>
#include <iostream>
#include <set>

#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include <AMReX_DataServices.H>
#include <AMReX_BCRec.H>
#include <AMReX_Interpolater.H>
#include <WritePlotFile.H>

#include <AMReX_BLFort.H>
#include <mechanism.h>
#include <chemistry_file.H>
#include <EOS.H>
#include <util.H>

using namespace amrex;

static
void 
print_usage (int,
             char* argv[])
{
  std::cerr << "usage:\n";
  std::cerr << argv[0] << " infile infile=f1 [options] \n\tOptions:\n";
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
  Initialize(argc,argv);
  {
    if (argc < 2)
      print_usage(argc,argv);

    ParmParse pp;

    if (pp.contains("help"))
      print_usage(argc,argv);

    if (pp.contains("verbose"))
      AmrData::SetVerbose(true);

    std::string plotFileName; pp.get("infile",plotFileName);
    DataServices::SetBatchMode();
    Amrvis::FileType fileType(Amrvis::NEWPLT);

    DataServices dataServices(plotFileName, fileType);
    if( ! dataServices.AmrDataOk()) {
      DataServices::Dispatch(DataServices::ExitRequest, NULL);
      // ^^^ this calls ParallelDescriptor::EndParallel() and exit()
    }
    AmrData& amrData = dataServices.AmrDataRef();

    EOS::init();

    int finestLevel = amrData.FinestLevel();
    pp.query("finestLevel",finestLevel);
    int Nlev = finestLevel + 1;

    int idYin = -1;
    int idTin = -1;
    Vector<std::string> spec_names;
    EOS::speciesNames(spec_names);
    const Vector<std::string>& plotVarNames = amrData.PlotVarNames();
    const std::string spName= "Y(" + spec_names[0] + ")";
    const std::string TName = "temp";
    for (int i=0; i<plotVarNames.size(); ++i)
    {
      if (plotVarNames[i] == spName) idYin = i;
      if (plotVarNames[i] == TName) idTin = i;
    }
    if (idYin<0 || idTin<0)
      Print() << "Cannot find required data in pltfile" << std::endl;

    const int idTout = 0;
    const int iddTout = idTout + 1;
    const int nCompOut = iddTout + 1;
    const int nCompIn = NUM_SPECIES + 1;

    Vector<std::string> outNames(nCompOut);
    Vector<std::string> inNames(nCompIn);
    Vector<int> destFillComps(nCompIn);
    const int idYlocal = 0; // Ys start here
    const int idTlocal = NUM_SPECIES; // T start here
    for (int i=0; i<NUM_SPECIES; ++i)
    {
      destFillComps[i] = idYlocal + i;
      inNames[i] =  "Y(" + spec_names[i] + ")";
    }
    destFillComps[idTlocal] = idTlocal;
    inNames[idTlocal] = TName;
    outNames[idTout] = TName;
    outNames[iddTout] = "d" + TName;
    
    Vector<std::unique_ptr<MultiFab>> outdata(Nlev);
    const int nGrow = 0;

    for (int lev=0; lev<Nlev; ++lev)
    {
      const BoxArray ba = amrData.boxArray(lev);
      const DistributionMapping dm(ba);
      outdata[lev].reset(new MultiFab(ba,dm,nCompOut,nGrow));
      MultiFab indata(ba,dm,nCompIn,nGrow);

      Print() << "Reading data for level " << lev << std::endl;
      amrData.FillVar(indata,lev,inNames,destFillComps);
      Print() << "Data has been read for level " << lev << std::endl;

      for (MFIter mfi(indata,TilingIfNotGPU()); mfi.isValid(); ++mfi)
      {
        const Box& bx = mfi.tilebox();
        Array4<Real> const& Y = indata.array(mfi,idYlocal);
        Array4<Real> const& Tin = indata.array(mfi,idTlocal);
        Array4<Real> const& Tout = (*outdata[lev]).array(mfi,idTout);
        Array4<Real> const& dTout = (*outdata[lev]).array(mfi,iddTout);

        AMREX_PARALLEL_FOR_3D ( bx, i, j, k,
        {
          Real Yl[NUM_SPECIES];
          Real Xl[NUM_SPECIES];
          for (int n=0; n<NUM_SPECIES; ++n) {
            Yl[n] = Y(i,j,k,idYlocal+n);
          }
          Real H;
          EOS::TY2H(Tin(i,j,k),Yl,H);

          Real Tsolve = 300;
          EOS::HY2T(H,Yl,Tsolve);
          Tout(i,j,k) = Tsolve;
          dTout(i,j,k) = Tsolve - Tin(i,j,k);
        });
      }

      Print() << "Derive finished for level " << lev << std::endl;
    }

    std::string outfile(getFileRoot(plotFileName) + "_T");
    Print() << "Writing new data to " << outfile << std::endl;
    const bool verb = false;
    WritePlotFile(GetVecOfPtrs(outdata),amrData,outfile,verb,outNames);
  }
  Finalize();
  return 0;
}
