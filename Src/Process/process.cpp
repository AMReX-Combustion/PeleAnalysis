#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include <AMReX_DataServices.H>
#include <AMReX_WritePlotFile.H>

using namespace amrex;

static
void 
print_usage (int,
             char* argv[])
{
  std::cerr << "usage:\n";
  std::cerr << argv[0] << " infile=f1 [options] \n\tOptions:\n";
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

    // Open plotfile header and create an amrData object pointing into it
    std::string plotFileName; pp.get("infile",plotFileName);
    DataServices::SetBatchMode();
    Amrvis::FileType fileType(Amrvis::NEWPLT);
    DataServices dataServices(plotFileName, fileType);
    if( ! dataServices.AmrDataOk()) {
      DataServices::Dispatch(DataServices::ExitRequest, NULL);
    }
    AmrData& amrData = dataServices.AmrDataRef();

    // Set up input field data names, and destination components to load data upon read.
    Vector<std::string> inNames = {{AMREX_D_DECL("velx","vely","velz"),"density"}};
    //Vector<std::string> inNames = {{AMREX_D_DECL("x_velocity","y_velocity","z_velocity"),"density"}};
    Vector<int> destFillComps   = {{AMREX_D_DECL(0, 1, 2), AMREX_SPACEDIM}};
    const int idVlocal = destFillComps[0]; // Vs start here
    const int idRlocal = destFillComps[AMREX_SPACEDIM]; // R starts here

    int idVin = -1;
    int idRin = -1;
    const Vector<std::string>& plotVarNames = amrData.PlotVarNames();
    for (int i=0; i<plotVarNames.size(); ++i)
    {
      if (plotVarNames[i] == inNames[idVlocal]) idVin = i;
      if (plotVarNames[i] == inNames[idRlocal]) idRin = i;
    }
    if (idVin<0 || idRin<0) {
      Print() << "Cannot find required data in pltfile " << idVin << " " << idRin << std::endl;
      Abort();
    }

    // Loop over AMR levels in the plotfile, read the data and do work
    int finestLevel = amrData.FinestLevel();
    pp.query("finestLevel",finestLevel);
    int Nlev = finestLevel + 1;
    const int nGrow = 0;

    Vector<std::string> outNames = {{"momx","momy","momz"}};
    Vector<std::unique_ptr<MultiFab>> outdata(Nlev);

    for (int lev=0; lev<Nlev; ++lev)
    {
      const BoxArray ba = amrData.boxArray(lev);
      const DistributionMapping dm(ba);
      MultiFab indata(ba,dm,inNames.size(),nGrow); // Container to hold plotfile data
      outdata[lev].reset(new MultiFab(ba,dm,outNames.size(),nGrow)); // Container to hold derived data to be written out in new plotfile
      
      Print() << "Reading data for level " << lev << std::endl;
      amrData.FillVar(indata,lev,inNames,destFillComps);
      Print() << "Data has been read for level " << lev << std::endl;

      for (MFIter mfi(indata,TilingIfNotGPU()); mfi.isValid(); ++mfi)
      {
        const Box& bx = mfi.tilebox();
        Array4<Real> const& Vel = indata.array(mfi,idVlocal);
        Array4<Real> const& Rho = indata.array(mfi,idRlocal);
        Array4<Real> const& Mom = outdata[lev]->array(mfi);

        AMREX_PARALLEL_FOR_3D ( bx, i, j, k,
        {
          for (int n=0; n<AMREX_SPACEDIM; ++n) {
            Mom(i,j,k,n) = Rho(i,j,k) * Vel(i,j,k,n);
          }
        });
      }

      Print() << "Derive finished for level " << lev << std::endl;
    }

    std::string outfile(getFileRoot(plotFileName) + "_mom");
    Print() << "Writing new data to " << outfile << std::endl;
    const bool verb = false;
    WritePlotFile(GetVecOfPtrs(outdata),amrData,outfile,verb,outNames);
    
  }
  Finalize();
  return 0;
}
