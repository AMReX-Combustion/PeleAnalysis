#include <string>
#include <iostream>
#include <set>

#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include <AMReX_DataServices.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_VisMF.H>

using namespace amrex;

static void
print_usage(int, char* argv[])
{
  std::cerr
    << "Calculates a progress variable and its source term from a pltFile as C "
       "= (Sum(specNames) - unburntVal)/(burntVal - unburntVal)\n";
  std::cerr << "usage:\n";
  std::cerr << argv[0]
            << "inputs infile=<i> speciesNames=<s> unburntVal=<u> burntVal=<b> "
               "[options] \n\tOptions:\n";
  std::cerr << "\t     infile=<i> where <i> is a pltfile\n";
  std::cerr << "\t     speciesNames=<s> where <s> are all the species to be "
               "included in the progress variable\n";
  std::cerr << "\t     unburntVal=float This is the sum of the included "
               "species in the unburnt unburntVal[DEF->0]\n";
  std::cerr << "\t     burntVal=float This is the value of included species in "
               "the burnt comp[DEF->1]\n";
  std::cerr << "\t     outsuffix=string This is the ending of the output "
               "pltfile[DEF->_prog]\n";
  std::cerr << "\t     outname=string This is the name of the variable in the "
               "output [DEF->progVar]\n";
  std::cerr << "\t     printSource=int Set to 0 to skip computing and writing "
               "the source term I_R(<outname>) [DEF->1]\n";
  exit(1);
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
    if (argc < 2)
      print_usage(argc, argv);

    ParmParse pp;

    if (pp.contains("help"))
      print_usage(argc, argv);

    if (pp.contains("verbose"))
      AmrData::SetVerbose(true);

    std::string plotFileName;
    pp.get("infile", plotFileName);

    Vector<std::string> species;
    int nSpec = pp.countval("speciesNames");
    species.resize(nSpec);
    pp.getarr("speciesNames", species, 0, nSpec);

    Real unburnt = 0;
    pp.get("unburntVal", unburnt);
    Real burnt = 1;
    pp.get("burntVal", burnt);

    if (burnt == unburnt) {
      amrex::Abort(
        "burntVal must differ from unburntVal (denominator is zero)");
    }

    std::string outsuffix = "_prog";
    pp.query("outsuffix", outsuffix);
    std::string outname = "progVar";
    pp.query("outname", outname);

    // By default the progress-variable source term I_R(<outname>) is written.
    // Set printSource=0 to output only specSum and the progress variable.
    int printSource = 1;
    pp.query("printSource", printSource);

    // Auxiliary variables: names of variables copied unchanged from the input
    // plotfile to the output plotfile (see Aux_Variables in the documentation).
    int nAuxVar = pp.countval("Aux_Variables");
    Vector<std::string> auxVar(nAuxVar);
    for (int ivar = 0; ivar < nAuxVar; ++ivar) {
      pp.get("Aux_Variables", auxVar[ivar], ivar);
    }

    DataServices::SetBatchMode();
    Amrvis::FileType fileType(Amrvis::NEWPLT);

    DataServices dataServices(plotFileName, fileType);
    if (!dataServices.AmrDataOk()) {
      DataServices::Dispatch(DataServices::ExitRequest, NULL);
    }
    AmrData& amrData = dataServices.AmrDataRef();

    int finestLevel = amrData.FinestLevel();
    pp.query("finestLevel", finestLevel);
    finestLevel = std::max(0, std::min(finestLevel, amrData.FinestLevel()));
    int Nlev = finestLevel + 1;

    Vector<int> is_per(AMREX_SPACEDIM, 1);
    pp.queryarr("is_per", is_per, 0, AMREX_SPACEDIM);
    Print() << "Periodicity assumed for this case: ";
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
      Print() << is_per[idim] << " ";
    }
    Print() << std::endl;

    RealBox rb(&(amrData.ProbLo()[0]), &(amrData.ProbHi()[0]));

    int nCompIn = amrData.NComp();
    Vector<std::string> inNames = amrData.PlotVarNames();

    // Resolve the component index of each auxiliary variable in the input file
    Vector<int> auxIdx(nAuxVar, -1);
    for (int ivar = 0; ivar < nAuxVar; ++ivar) {
      auxIdx[ivar] = amrData.StateNumber(auxVar[ivar]);
      if (auxIdx[ivar] < 0) {
        amrex::Abort("Unknown auxiliary variable name: " + auxVar[ivar]);
      }
    }

    // Get ids of relevant species (mass fractions)
    Vector<int> idY(nSpec, -1);
    for (int j = 0; j < nSpec; ++j) {
      for (int i = 0; i < (int)inNames.size(); ++i) {
        if (inNames[i] == species[j])
          idY[j] = i;
      }
      if (idY[j] < 0) {
        amrex::Abort("Species " + species[j] + " not found in plotfile");
      }
    }

    // Get ids of production rates I_R(<species>) for each species (only needed
    // when the source term is requested)
    Vector<int> idIR(nSpec, -1);
    if (printSource) {
      for (int j = 0; j < nSpec; ++j) {
        // Strip "Y(" prefix and ")" suffix if present, e.g. "Y(H2)" -> "H2"
        std::string specBase = species[j];
        if (
          specBase.size() > 2 && specBase.substr(0, 2) == "Y(" &&
          specBase.back() == ')') {
          specBase = specBase.substr(2, specBase.size() - 3);
        }

        std::string irName = "I_R(" + specBase + ")";
        for (int i = 0; i < (int)inNames.size(); ++i) {
          if (inNames[i] == irName)
            idIR[j] = i;
        }
        if (idIR[j] < 0) {
          amrex::Abort("Production rate " + irName + " not found in plotfile");
        }
      }
    }

    Vector<MultiFab> outdata(Nlev);
    Vector<Geometry> geoms(Nlev);
    int nGrow = 0;
    // Output: specSum, progVar, and optionally I_R(progVar)
    Vector<std::string> outNames = {"specSum", outname};
    if (printSource)
      outNames.push_back("I_R(" + outname + ")");

    // Auxiliary variables are appended, unchanged, after the computed fields
    const int nBaseOut = outNames.size();
    for (int ivar = 0; ivar < nAuxVar; ++ivar) {
      outNames.push_back(auxVar[ivar]);
    }

    Vector<int> destFillComps(nCompIn);
    for (int i = 0; i < nCompIn; ++i)
      destFillComps[i] = i;

    for (int lev = 0; lev < Nlev; ++lev) {

      const BoxArray ba = amrData.boxArray(lev);
      const DistributionMapping dm(ba);

      outdata[lev] = MultiFab(ba, dm, outNames.size(), nGrow);
      MultiFab indata(ba, dm, nCompIn, nGrow);

      int coord = amrData.CoordSys();
      geoms[lev] =
        Geometry(amrData.ProbDomain()[lev], &rb, coord, &(is_per[0]));

      Print() << "Reading data for level " << lev << std::endl;
      amrData.FillVar(indata, lev, inNames, destFillComps);
      Print() << "Data has been read for level " << lev << std::endl;

      // Copy species index arrays to device
      amrex::Gpu::DeviceVector<int> d_idY(idY.size());
      amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, idY.begin(), idY.end(), d_idY.begin());
      int const* idY_d = d_idY.data();

      // Copy production rate index arrays to device (only when writing source)
      amrex::Gpu::DeviceVector<int> d_idIR;
      int const* idIR_d = nullptr;
      if (printSource) {
        d_idIR.resize(idIR.size());
        amrex::Gpu::copy(
          amrex::Gpu::hostToDevice, idIR.begin(), idIR.end(), d_idIR.begin());
        idIR_d = d_idIR.data();
      }

      Real denom = burnt - unburnt;

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (MFIter mfi(indata, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();
        auto const& out_a = outdata[lev].array(mfi);
        auto const& in_a = indata.array(mfi);
        amrex::ParallelFor(
          bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            // Sum of species mass fractions
            Real sum = 0.0_rt;
            for (int s = 0; s < nSpec; ++s) {
              sum += in_a(i, j, k, idY_d[s]);
            }
            out_a(i, j, k, 0) = sum;

            // Progress variable: C = (sum - unburnt) / (burnt - unburnt)
            out_a(i, j, k, 1) = (sum - unburnt) / denom;

            // Production rate of progress variable:
            // I_R(C) = Sum( I_R(species) ) / (burnt - unburnt)
            if (printSource) {
              Real irSum = 0.0_rt;
              for (int s = 0; s < nSpec; ++s) {
                irSum += in_a(i, j, k, idIR_d[s]);
              }
              out_a(i, j, k, 2) = irSum / denom;
            }
          });
      }

      // Copy auxiliary variables unchanged from the input to the output
      for (int ivar = 0; ivar < nAuxVar; ++ivar) {
        MultiFab::Copy(
          outdata[lev], indata, auxIdx[ivar], nBaseOut + ivar, 1, nGrow);
      }
    }

    std::string outfile(getFileRoot(plotFileName) + outsuffix);

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
