#include <string>
#include <iostream>
#include <set>

#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include <AMReX_DataServices.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_VisMF.H>
#include <AMReX_BCRec.H>
#include <AMReX_Interpolater.H>
#include <AMReX_GpuLaunch.H>

#include <PelePhysics.H>

using namespace amrex;

static void
print_usage(int, char* argv[])
{
  std::cerr << "Usage:\n"
            << "  " << argv[0] << " infile=FILE [OPTIONS]\n\n"

            << "Required arguments:\n"
            << "  infile=FILE        AMReX plotfile\n\n"

            << "Options:\n"
            << "  -h, --help         Show this help message\n\n"
            << "Visit PeleAnalysis/Src/InputSamples for examples or refer to "
            << "the documentation.\n";

  std::exit(1);
}

std::string
getFileRoot(const std::string& infile)
{
  std::vector<std::string> tokens = Tokenize(infile, std::string("/"));
  return tokens[tokens.size() - 1];
}

int
main(int argc, char* argv[])
{
  Initialize(argc, argv);
  {
    if (argc < 2) {
      print_usage(argc, argv);
    } else if (
      (std::strcmp(argv[1], "-h") == 0) ||
      (std::strcmp(argv[1], "--help") == 0)) {
      print_usage(argc, argv);
    }

    ParmParse pp;

    if (pp.contains("verbose"))
      AmrData::SetVerbose(true);

    std::string plotFileName;
    pp.get("infile", plotFileName);
    DataServices::SetBatchMode();
    Amrvis::FileType fileType(Amrvis::NEWPLT);

    DataServices dataServices(plotFileName, fileType);
    if (!dataServices.AmrDataOk()) {
      DataServices::Dispatch(DataServices::ExitRequest, NULL);
      // ^^^ this calls ParallelDescriptor::EndParallel() and exit()
    }
    AmrData& amrData = dataServices.AmrDataRef();

    int finestLevel = amrData.FinestLevel();
    pp.query("finestLevel", finestLevel);
    int Nlev = finestLevel + 1;

    // Auxiliary variables: names copied unchanged from input to output plotfile
    int nAuxVar = pp.countval("Aux_Variables");
    Vector<std::string> auxVar(nAuxVar);
    for (int ivar = 0; ivar < nAuxVar; ++ivar) {
      pp.get("Aux_Variables", auxVar[ivar], ivar);
    }

    int idYin = -1;
    int idTin = -1;
    Vector<std::string> spec_names;
    pele::physics::eos::speciesNames<pele::physics::PhysicsType::eos_type>(
      spec_names);
    auto eos = pele::physics::PhysicsType::eos();
    const Vector<std::string>& plotVarNames = amrData.PlotVarNames();
    const std::string spName = "Y(" + spec_names[0] + ")";
    const std::string TName = "Temp";
    for (int i = 0; i < plotVarNames.size(); ++i) {
      if (plotVarNames[i] == spName)
        idYin = i;
      if (plotVarNames[i] == TName)
        idTin = i;
    }
    if (idYin < 0 || idTin < 0)
      Print() << "Cannot find required data in pltfile" << std::endl;

    const int idXout = 0;
    const int idTout = NUM_SPECIES;
    const int idAuxLocal = NUM_SPECIES + 1; // aux vars start here in input
    const int idAuxOut = idTout + 1;        // aux vars start here in output
    const int nCompIn = NUM_SPECIES + 1 + nAuxVar;
    const int nCompOut = idXout + NUM_SPECIES + 1 + nAuxVar;

    Vector<std::string> outNames(nCompOut);
    Vector<std::string> inNames(nCompIn);
    Vector<int> destFillComps(nCompIn);
    const int idYlocal = 0;           // Xs start here
    const int idTlocal = NUM_SPECIES; // T start here
    for (int i = 0; i < NUM_SPECIES; ++i) {
      destFillComps[i] = idYlocal + i;
      inNames[i] = "Y(" + spec_names[i] + ")";
      outNames[i] = "X(" + spec_names[i] + ")";
    }
    destFillComps[idTlocal] = idTlocal;
    inNames[idTlocal] = TName;
    outNames[idTout] = TName;

    // Auxiliary variables are read into the input MultiFab and appended,
    // unchanged, to the output plotfile after the computed fields
    for (int ivar = 0; ivar < nAuxVar; ++ivar) {
      if (amrData.StateNumber(auxVar[ivar]) < 0) {
        amrex::Abort("Unknown auxiliary variable name: " + auxVar[ivar]);
      }
      destFillComps[idAuxLocal + ivar] = idAuxLocal + ivar;
      inNames[idAuxLocal + ivar] = auxVar[ivar];
      outNames[idAuxOut + ivar] = auxVar[ivar];
    }

    Vector<std::unique_ptr<MultiFab>> outdata(Nlev);
    Vector<Geometry> geoms(Nlev);
    amrex::RealBox real_box(
      {AMREX_D_DECL(
        amrData.ProbLo()[0], amrData.ProbLo()[1], amrData.ProbLo()[2])},
      {AMREX_D_DECL(
        amrData.ProbHi()[0], amrData.ProbHi()[1], amrData.ProbHi()[2])});
    amrex::Array<int, AMREX_SPACEDIM> is_periodic{AMREX_D_DECL(1, 1, 1)};
    geoms[0] = amrex::Geometry(
      (amrData.ProbDomain())[0], real_box, amrData.CoordSys(), is_periodic);
    const int nGrow = 0;

    for (int lev = 0; lev < Nlev; ++lev) {
      const BoxArray ba = amrData.boxArray(lev);
      const DistributionMapping dm(ba);
      outdata[lev].reset(new MultiFab(ba, dm, nCompOut, nGrow));
      if (lev > 0) {
        geoms[lev] = amrex::refine(geoms[lev - 1], 2);
      }
      MultiFab indata(ba, dm, nCompIn, nGrow);

      Print() << "Reading data for level " << lev << std::endl;
      amrData.FillVar(indata, lev, inNames, destFillComps);
      Print() << "Data has been read for level " << lev << std::endl;

      for (MFIter mfi(indata, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();
        Array4<Real> const& Y = indata.array(mfi);
        Array4<Real> const& Tin = indata.array(mfi);
        Array4<Real> const& X = (*outdata[lev]).array(mfi);
        Array4<Real> const& Tout = (*outdata[lev]).array(mfi);

        AMREX_PARALLEL_FOR_3D(bx, i, j, k, {
          Real Yl[NUM_SPECIES];
          Real Xl[NUM_SPECIES];
          for (int n = 0; n < NUM_SPECIES; ++n) {
            Yl[n] = Y(i, j, k, idYlocal + n);
          }
          eos.Y2X(Yl, Xl);
          for (int n = 0; n < NUM_SPECIES; ++n) {
            X(i, j, k, idXout + n) = Xl[n];
          }
          Tout(i, j, k, idTout) = Tin(i, j, k, idTlocal);
        });
      }

      // Copy auxiliary variables unchanged from the input to the output
      for (int ivar = 0; ivar < nAuxVar; ++ivar) {
        MultiFab::Copy(
          *outdata[lev], indata, idAuxLocal + ivar, idAuxOut + ivar, 1, nGrow);
      }

      Print() << "Derive finished for level " << lev << std::endl;
    }

    std::string outfile(getFileRoot(plotFileName) + "_X");

    // Cap the number of plotfile data files via the n_files option (AMReX)
    int n_files = amrex::VisMF::GetNOutFiles();
    pp.query("n_files", n_files);
    amrex::VisMF::SetNOutFiles(n_files);

    Print() << "Writing new data to " << outfile << std::endl;
    Vector<int> isteps(Nlev, 0);
    Vector<IntVect> refRatios(Nlev - 1, {AMREX_D_DECL(2, 2, 2)});
    amrex::WriteMultiLevelPlotfile(
      outfile, Nlev, GetVecOfConstPtrs(outdata), outNames, geoms, 0.0, isteps,
      refRatios);
  }
  Finalize();
  return 0;
}
