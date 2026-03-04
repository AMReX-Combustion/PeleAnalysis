#include <iostream>
#include <set>
#include <string>

#include <AMReX_DataServices.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_PlotFileUtilHDF5.H>

static void
print_usage(int, char* argv[])
{
  std::cerr << "Usage:\n"
            << "  " << argv[0] << " infile=FILE [OPTIONS]\n\n"

            << "Required arguments:\n"
            << "  infile=FILE              AMReX plotfile to convert\n\n"

            << "Options:\n"
            << "  -h, --help               Show this help message\n\n"

            << "Visit PeleAnalysis/Src/InputSamples for examples or refer to "
            << "the documentation.\n";

  std::exit(1);
}

std::string
getFileRoot(const std::string& infile)
{
  std::vector<std::string> tokens = amrex::Tokenize(infile, std::string("/"));
  return tokens[tokens.size() - 1];
}

int
main(int argc, char* argv[])
{
  amrex::Initialize(argc, argv);
  {
    if (argc < 2) {
      print_usage(argc, argv);
    } else if (
      (std::strcmp(argv[1], "-h") == 0) ||
      (std::strcmp(argv[1], "--help") == 0)) {
      print_usage(argc, argv);
    }
    amrex::Real dRunTime1 = amrex::ParallelDescriptor::second();

    int finestLevel = 1000;
    amrex::ParmParse pp;

    pp.query("finestLevel", finestLevel);
    std::string plotFileName;
    pp.get("infile", plotFileName);

    // Initialize DataService
    amrex::DataServices::SetBatchMode();
    amrex::Amrvis::FileType fileType(amrex::Amrvis::NEWPLT);
    amrex::DataServices dataServices(plotFileName, fileType);
    if (!dataServices.AmrDataOk()) {
      amrex::DataServices::Dispatch(amrex::DataServices::ExitRequest, NULL);
    }
    amrex::AmrData& amrData = dataServices.AmrDataRef();

    // Plotfile global infos
    finestLevel = std::min(finestLevel, amrData.FinestLevel());
    int Nlev = finestLevel + 1;
    const amrex::Vector<std::string>& plotVarNames = amrData.PlotVarNames();
    int nvars = plotVarNames.size();
    int id_comp_last = 0;
    amrex::RealBox rb(&(amrData.ProbLo()[0]), &(amrData.ProbHi()[0]));
    amrex::Vector<int> is_per(AMREX_SPACEDIM, 0);
    int coord = 0;

    // Copy data
    amrex::Vector<amrex::MultiFab*> fileData(Nlev);
    amrex::Vector<amrex::Geometry> geoms(Nlev);
    const int nGrow = 1;

    // Read data on all the levels
    for (int lev = 0; lev < Nlev; ++lev) {
      const amrex::DistributionMapping dm(amrData.boxArray(lev));
      geoms[lev] =
        amrex::Geometry(amrData.ProbDomain()[lev], &rb, coord, &(is_per[0]));
      fileData[lev] = new amrex::MultiFab(amrData.boxArray(lev), dm, nvars, 0);
    }

    amrex::Vector<int> idcomp;
    for (int i = 0; i < plotVarNames.size();
         ++i) {            // loop though current pltfile variable names
      idcomp.push_back(i); // sets index location of plotVarName that matches
    }

    for (int lev = 0; lev < Nlev; ++lev) {
      for (int i = 0; i < nvars; ++i) {
        fileData[lev]->ParallelCopy(
          amrData.GetGrids(lev, idcomp[i]), 0, id_comp_last + i, 1);
      }
    }

    // Write the results
    std::string hdf5_compression{"ZLIB@9"};
    pp.query("hdf5_compression", hdf5_compression);
    std::string outfile(getFileRoot(plotFileName) + "_hdf5");
    amrex::Print() << "Writing new data to " << outfile << std::endl;
    amrex::Vector<int> isteps(Nlev, 0);
    amrex::Vector<amrex::IntVect> refRatios(Nlev - 1, {AMREX_D_DECL(2, 2, 2)});
    amrex::WriteMultiLevelPlotfileHDF5SingleDset(
      outfile, Nlev, amrex::GetVecOfConstPtrs(fileData), plotVarNames, geoms,
      amrData.Time(), isteps, refRatios, hdf5_compression);

    amrex::Real dRunTime2 = amrex::ParallelDescriptor::second();
    amrex::Real runtime_total = dRunTime2 - dRunTime1;

    const int IOProc = amrex::ParallelDescriptor::IOProcessorNumber();
    amrex::ParallelDescriptor::ReduceRealMax(runtime_total, IOProc);

    if (amrex::ParallelDescriptor::IOProcessor()) {
      amrex::Print() << "Run time = " << runtime_total << std::endl;
    }
  }

  amrex::Finalize();
  return 0;
}
