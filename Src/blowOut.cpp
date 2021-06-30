#include <string>
#include <iostream>
#include <set>

#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include <AMReX_DataServices.H>
#include <WritePlotFile.H>

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

extern "C" {
    void process(const int* lo,  const int* hi,
                 const int* dlo, const int* dhi,
                 Real* U, const int* Ulo, const int* Uhi,
                 const int* nc, const amrex::Real* plo,  const amrex::Real* dx);
};

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

    int finestLevel = amrData.FinestLevel();
    pp.query("finestLevel",finestLevel);
    int Nlev = finestLevel + 1;

    Vector<std::unique_ptr<MultiFab>> outdata(Nlev);
    int ngrow = 0;
    int ncomp = amrData.NComp();

    const auto& domain = amrData.ProbDomain()[finestLevel];
    BoxArray ba(domain);

    int max_grid_size = 32;
    pp.query("max_grid_size",max_grid_size);
    ba.maxSize(max_grid_size);
    const DistributionMapping dm(ba);
    outdata[0].reset(new MultiFab(ba,dm,ncomp,ngrow));

    Vector<int> destFillComps(ncomp);
    for (int i=0; i<ncomp; ++i) destFillComps[i] = i;

    Print() << "Reading data" << std::endl;
    amrData.FillVar(*outdata[0],finestLevel,amrData.PlotVarNames(),destFillComps);
    Print() << "Data has been read" << std::endl;

    int coord = 0;
    RealBox rb(&(amrData.ProbLo()[0]),
               &(amrData.ProbHi()[0]));
    Vector<int> is_per(AMREX_SPACEDIM,0);
    Geometry geom(amrData.ProbDomain()[finestLevel],&rb,coord,&(is_per[0]));

    for (MFIter mfi(*outdata[0],TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
      const Box& box = mfi.tilebox();
      auto& fab = (*outdata[0])[mfi];
      process(BL_TO_FORTRAN_BOX(box),
              BL_TO_FORTRAN_BOX(domain),
              BL_TO_FORTRAN_ANYD(fab),
              &ncomp,geom.ProbLo(),geom.CellSize());
    }

    std::string outfile(getFileRoot(plotFileName) + "_b");
    Print() << "Writing new data to " << outfile << std::endl;
    bool verb = false;
    WritePlotFile(GetVecOfPtrs(outdata),amrData,outfile,verb,amrData.PlotVarNames());
  }
  Finalize();
  return 0;
}
