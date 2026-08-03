#include <string>
#include <iostream>
#include <set>
#include <cstring>

#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include <AMReX_DataServices.H>
#include <AMReX_VisMF.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_PlotFileUtil.H>
#include <PelePhysics.H>

using namespace amrex;

static void
print_usage(int, char* argv[])
{
  std::cerr
    << "Computes the Bilger mixture fraction Z from an AMReX plotfile, and\n"
       "optionally the mass fraction of every element in the mechanism.\n"
       "The fuel stream may be a single species or a blend.\n\n"

    << "Usage:\n"
    << "  " << argv[0] << " infile=FILE [OPTIONS]\n\n"

    << "Required arguments:\n"
    << "  infile=FILE             AMReX plotfile holding Y(<species>) for "
       "every\n"
    << "                          species in the compiled mechanism\n\n"

    << "Fuel stream (pick one; DEF: pure H2):\n"
    << "  fuelName=NAME           Single fuel species\n"
    << "  fuelNames=\"N1 N2 ...\"   Blend components, with one of:\n"
    << "  fuelMoleFracs=\"X1 ...\"    their mole fractions\n"
    << "  fuelMassFracs=\"Y1 ...\"    their mass fractions\n"
    << "                          Either is normalised to sum to one\n\n"

    << "Options:\n"
    << "  YO2ox=F                 O2 mass fraction of the oxidizer (DEF: "
       "0.233)\n"
    << "  YN2ox=F                 N2 mass fraction of the oxidizer (DEF: "
       "0.767)\n"
    << "  elementMassFracs=0|1    Also write Z_<elem> for each mechanism "
       "element (DEF: 0)\n"
    << "  outsuffix=STR           Suffix appended to the input name (DEF: "
       "_ZC)\n"
    << "  finestLevel=N           Finest AMR level used (DEF: finest in "
       "file)\n"
    << "  is_per=I J [K]          Periodicity flags per direction (DEF: 1 1 "
       "1)\n"
    << "  Aux_Variables=\"V1 ...\"  Variables copied unchanged to the output\n"
    << "  n_files=N               Max number of plotfile data files written\n"
    << "  -h, --help              Show this help message\n\n"

    << "For a reaction progress variable, use the progVar tool.\n"
    << "Visit PeleAnalysis/Src/InputSamples for examples or refer to "
       "the documentation.\n";

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

    if (pp.contains("help"))
      print_usage(argc, argv);

    if (pp.contains("verbose"))
      AmrData::SetVerbose(true);

    std::string plotFileName;
    pp.get("infile", plotFileName);

    // Fuel stream. Either a single species or a blend given as fuelNames plus
    // one of fuelMoleFracs / fuelMassFracs. The component list is settled first
    // so that reading, validating and normalising the fractions happens once,
    // on every path, rather than being duplicated per case.
    Vector<std::string> fuelNames;
    if (pp.countval("fuelNames") > 0) {
      pp.getarr("fuelNames", fuelNames);
    } else {
      std::string fuelName = "H2";
      pp.query("fuelName", fuelName);
      fuelNames = {fuelName};
    }
    const int nFuel = static_cast<int>(fuelNames.size());

    Vector<Real> fuelFracs;
    bool fuelFracsAreMole = true;
    const int nMole = pp.countval("fuelMoleFracs");
    const int nMass = pp.countval("fuelMassFracs");
    if (nMole > 0 && nMass > 0) {
      amrex::Abort("Give either fuelMoleFracs or fuelMassFracs, not both");
    }
    if (nMole > 0) {
      pp.getarr("fuelMoleFracs", fuelFracs);
    } else if (nMass > 0) {
      pp.getarr("fuelMassFracs", fuelFracs);
      fuelFracsAreMole = false;
    } else if (nFuel == 1) {
      // A lone component needs no fraction: any positive value normalises to 1
      // in the loop below.
      fuelFracs = {1.0};
    } else {
      amrex::Abort("fuelNames with more than one entry needs fuelMoleFracs or "
                   "fuelMassFracs");
    }

    if (static_cast<int>(fuelFracs.size()) != nFuel) {
      amrex::Abort("fuelNames and the fuel fractions must have equal length");
    }
    Real sumF = 0.0;
    for (const Real x : fuelFracs) {
      if (x < 0.0) {
        amrex::Abort("Fuel fractions must not be negative");
      }
      sumF += x;
    }
    if (sumF <= 0.0) {
      amrex::Abort("Fuel fractions must not sum to zero");
    }
    for (Real& x : fuelFracs) {
      x /= sumF;
    }

    // Oxidizer stream composition (mass fractions); defaults to air.
    Real YO2ox = 0.233;
    pp.query("YO2ox", YO2ox);
    Real YN2ox = 0.767;
    pp.query("YN2ox", YN2ox);
    std::string outsuffix = "_ZC";
    pp.query("outsuffix", outsuffix);
    int elementMassFracs = 0;
    pp.query("elementMassFracs", elementMassFracs);
    Vector<int> is_per(AMREX_SPACEDIM, 1);
    pp.queryarr("is_per", is_per, 0, AMREX_SPACEDIM);
    DataServices::SetBatchMode();
    Amrvis::FileType fileType(Amrvis::NEWPLT);

    DataServices dataServices(plotFileName, fileType);
    if (!dataServices.AmrDataOk()) {
      DataServices::Dispatch(DataServices::ExitRequest, NULL);
    }
    AmrData& amrData = dataServices.AmrDataRef();

    int finestLevel = amrData.FinestLevel();
    pp.query("finestLevel", finestLevel);
    int Nlev = finestLevel + 1;
    Vector<std::string> spec_names;
    pele::physics::eos::speciesNames<pele::physics::PhysicsType::eos_type>(
      spec_names);
    auto eos = pele::physics::PhysicsType::eos();

    // The CHON helpers this tool used to call raise an error for these two, and
    // the general mechanism calls below would instead return meaningless
    // numbers, so the check has to be made explicitly.
    const std::string eosName =
      pele::physics::PhysicsType::eos_type::identifier();
    if (eosName == "GammaLaw" || eosName == "Manifold") {
      amrex::Abort(
        "computeMixtureFraction needs a chemical mechanism, but the compiled "
        "EOS is " +
        eosName);
    }

    // The element set comes from the mechanism: its size, its names and its
    // order are all mechanism properties and none of them is CHON in general.
    // drm19's own order is O H C N Ar, alzeta carries F, IonizedAir carries a
    // free electron, and Ar or He appear in about half of the bundled
    // mechanisms. Assuming a CHON quartet silently drops every one of those.
    Vector<std::string> elemNames;
    CKSYME_STR(elemNames);
    Real atw[NUM_ELEMENTS];
    CKAWT(atw);
    int ecomp[NUM_SPECIES * NUM_ELEMENTS];
    CKNCF(ecomp);

    // Auxiliary variables: names copied unchanged from input to output plotfile
    int nAuxVar = pp.countval("Aux_Variables");
    Vector<std::string> auxVar(nAuxVar);
    for (int ivar = 0; ivar < nAuxVar; ++ivar) {
      pp.get("Aux_Variables", auxVar[ivar], ivar);
    }

    // Output layout: Z, then one elemental mass fraction per element in the
    // mechanism when they are requested, then the auxiliary variables.
    constexpr int idZlocal = 0; // Z out here
    const int idElemOut = 1;
    const int nElemOut = elementMassFracs ? NUM_ELEMENTS : 0;

    const int idAuxLocal = NUM_SPECIES; // aux vars start here in input
    const int idAuxOut = 1 + nElemOut;  // aux vars start here in output
    const int nCompIn = NUM_SPECIES + nAuxVar;
    const int nCompOut = 1 + nElemOut + nAuxVar;
    Vector<std::string> outNames(nCompOut);
    Vector<std::string> inNames(nCompIn);
    Vector<int> destFillComps(nCompIn);

    for (int i = 0; i < NUM_SPECIES; ++i) {
      destFillComps[i] = i;
      inNames[i] = "Y(" + spec_names[i] + ")";
    }
    // out
    outNames[idZlocal] = "Z";
    for (int e = 0; e < nElemOut; ++e) {
      outNames[idElemOut + e] = "Z_" + elemNames[e];
    }

    // Auxiliary variables are read into the input MultiFab and appended,
    // unchanged, to the output plotfile after the computed field
    for (int ivar = 0; ivar < nAuxVar; ++ivar) {
      if (amrData.StateNumber(auxVar[ivar]) < 0) {
        amrex::Abort("Unknown auxiliary variable name: " + auxVar[ivar]);
      }
      destFillComps[idAuxLocal + ivar] = idAuxLocal + ivar;
      inNames[idAuxLocal + ivar] = auxVar[ivar];
      outNames[idAuxOut + ivar] = auxVar[ivar];
    }

    Vector<MultiFab> outdata(Nlev);
    Vector<Geometry> geoms(Nlev);
    RealBox rb(&(amrData.ProbLo()[0]), &(amrData.ProbHi()[0]));
    constexpr int nGrow = 0;

    Real Zfu = -1.0;
    Real Zox = -1.0;

    Real YF[NUM_SPECIES], YO[NUM_SPECIES];
    for (int i = 0; i < NUM_SPECIES; ++i) {
      YF[i] = 0.0;
      YO[i] = 0.0;
      if (spec_names[i] == "O2") {
        YO[i] = YO2ox;
      }
      if (spec_names[i] == "N2") {
        YO[i] = YN2ox;
      }
    }

    // Locate each fuel component in the mechanism ordering used to fill indata
    // (destFillComps[i] = i, inNames[i] = Y(spec_names[i])).
    Vector<int> fuelIdx(nFuel, -1);
    for (int f = 0; f < nFuel; ++f) {
      for (int i = 0; i < NUM_SPECIES; ++i) {
        if (spec_names[i] == fuelNames[f]) {
          fuelIdx[f] = i;
        }
      }
      if (fuelIdx[f] < 0) {
        amrex::Abort(
          "Fuel species " + fuelNames[f] + " not found in mechanism");
      }
      for (int g = 0; g < f; ++g) {
        if (fuelIdx[g] == fuelIdx[f]) {
          amrex::Abort("Fuel species " + fuelNames[f] + " listed twice");
        }
      }
    }

    // Bilger's coupling function is defined on C, H and O:
    //   beta = 2 Z_C / W_C + Z_H / (2 W_H) - Z_O / W_O
    // Each of the three is located by name in the mechanism's element list
    // rather than assumed to occupy a fixed slot; an element the mechanism does
    // not carry simply contributes nothing. Every other element - N, Ar, He, a
    // free electron - has zero weight, which is correct: they are diluents that
    // the coupling function is meant to ignore.
    Array<amrex::Real, NUM_ELEMENTS> Beta_mix;
    for (int e = 0; e < NUM_ELEMENTS; ++e) {
      Beta_mix[e] = 0.0;
      if (atw[e] <= 0.0) {
        continue;
      }
      if (elemNames[e] == "C") {
        Beta_mix[e] = 2.0 / atw[e];
      } else if (elemNames[e] == "H") {
        Beta_mix[e] = 1.0 / (2.0 * atw[e]);
      } else if (elemNames[e] == "O") {
        Beta_mix[e] = -1.0 / atw[e];
      }
    }

    amrex::Array<amrex::Real, NUM_SPECIES> spec_Bilger_fact;
    amrex::Real mwt[NUM_SPECIES];
    eos.molecular_weight(mwt);

    // Fuel stream mass fractions. Mole fractions need converting with
    //   Y_k = X_k W_k / sum_j (X_j W_j)
    // whereas mass fractions were already normalised to sum to one above and
    // are therefore the answer as they stand.
    if (fuelFracsAreMole) {
      Real denom = 0.0;
      for (int f = 0; f < nFuel; ++f) {
        denom += fuelFracs[f] * mwt[fuelIdx[f]];
      }
      for (int f = 0; f < nFuel; ++f) {
        YF[fuelIdx[f]] = fuelFracs[f] * mwt[fuelIdx[f]] / denom;
      }
    } else {
      for (int f = 0; f < nFuel; ++f) {
        YF[fuelIdx[f]] = fuelFracs[f];
      }
    }
    Print() << "Fuel stream composition (mass fractions):\n";
    for (int f = 0; f < nFuel; ++f) {
      Print() << "  Y(" << fuelNames[f] << ") = " << YF[fuelIdx[f]] << "\n";
    }

    // Mass fraction of each element in each species, in mechanism element
    // order,
    //   w_elem(e,i) = n_{e,i} * A_e / W_i
    // so that Z_e = sum_i w_elem(e,i) Y_i. This is the same expression as
    // eos.Y2Z, but with the mechanism lookups hoisted out of the cell loop:
    // Y2Z rebuilds atomicWeight, get_imw and the whole NUM_SPECIES*NUM_ELEMENTS
    // composition matrix on every call, which is not something to do per cell.
    amrex::Array2D<amrex::Real, 0, NUM_ELEMENTS - 1, 0, NUM_SPECIES - 1> w_elem;
    for (int e = 0; e < NUM_ELEMENTS; ++e) {
      for (int i = 0; i < NUM_SPECIES; ++i) {
        w_elem(e, i) =
          (atw[e] > 0.0) ? ecomp[i * NUM_ELEMENTS + e] * atw[e] / mwt[i] : 0.0;
      }
    }

    Zfu = 0.0;
    Zox = 0.0;
    for (int i = 0; i < NUM_SPECIES; ++i) {
      spec_Bilger_fact[i] = 0.0;
      for (int e = 0; e < NUM_ELEMENTS; ++e) {
        // w_elem is exactly the n_{e,i} A_e / W_i factor this used to spell out
        // inline, so it is reused rather than recomputed.
        spec_Bilger_fact[i] += Beta_mix[e] * w_elem(e, i);
      }
      Zfu += spec_Bilger_fact[i] * YF[i];
      Zox += spec_Bilger_fact[i] * YO[i];
    }
    if (Zfu == Zox) {
      amrex::Abort("Fuel and oxidizer Bilger values are equal; check the fuel "
                   "stream and oxidizer composition (Zfu - Zox is zero)");
    }
    const Real denom_inv = 1.0 / (Zfu - Zox);
    Print() << "Zfu = " << Zfu << ", Zox = " << Zox << std::endl;

    // An element outside C, H and O is invisible to Bilger's coupling function.
    // That is the right answer while it stands alone - Ar, He and free
    // electrons are diluents - and it is also the right answer for N, which the
    // formulation excludes deliberately even in a mechanism with NOx chemistry.
    // It is not the right answer for anything else that is bonded to C, H or O,
    // because that element is then part of the combustion chemistry and Z
    // cannot see it. Warn only in that last case, so the message means
    // something when it appears: F in the alzeta mechanism triggers it, Ar and
    // NOx nitrogen do not.
    auto isBilgerElem = [&](int e) {
      return elemNames[e] == "C" || elemNames[e] == "H" || elemNames[e] == "O";
    };
    for (int e = 0; e < NUM_ELEMENTS; ++e) {
      if (isBilgerElem(e) || elemNames[e] == "N") {
        continue;
      }
      for (int i = 0; i < NUM_SPECIES; ++i) {
        if (ecomp[i * NUM_ELEMENTS + e] == 0) {
          continue;
        }
        bool bondedToCHO = false;
        for (int e2 = 0; e2 < NUM_ELEMENTS; ++e2) {
          if (isBilgerElem(e2) && ecomp[i * NUM_ELEMENTS + e2] > 0) {
            bondedToCHO = true;
            break;
          }
        }
        if (bondedToCHO) {
          Print() << "\n*** WARNING: element " << elemNames[e]
                  << " is bonded to C, H or O in " << spec_names[i] << ".\n"
                  << "    The Bilger coupling function is built on C, H and O "
                     "alone, so Z does not\n"
                  << "    account for " << elemNames[e] << " chemistry. Z_"
                  << elemNames[e] << " and the other elemental mass\n"
                  << "    fractions are unaffected.\n\n";
          break; // one warning per element is enough
        }
      }
    }

    for (int lev = 0; lev < Nlev; ++lev) {
      const BoxArray ba = amrData.boxArray(lev);
      const DistributionMapping dm(ba);
      MultiFab indata(ba, dm, nCompIn, nGrow);
      outdata[lev].define(ba, dm, nCompOut, nGrow);

      Print() << "Reading data for level " << lev << std::endl;
      amrData.FillVar(indata, lev, inNames, destFillComps);
      geoms[lev] = Geometry(
        amrData.ProbDomain()[lev], &rb, amrData.CoordSys(), &(is_per[0]));
      Print() << "Data has been read for level " << lev << std::endl;
      auto in_ma = indata.const_arrays();
      auto out_ma = outdata[lev].arrays();
      const int doElem = elementMassFracs;
      amrex::ParallelFor(
        outdata[lev],
        [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k) noexcept {
          Real Zloc = 0.0;
          for (int n = 0; n < NUM_SPECIES; ++n) {
            Zloc += in_ma[box_no](i, j, k, n) * spec_Bilger_fact[n];
          }
          out_ma[box_no](i, j, k, idZlocal) = (Zloc - Zox) * denom_inv;

          if (doElem) {
            for (int e = 0; e < NUM_ELEMENTS; ++e) {
              Real Ze = 0.0;
              for (int n = 0; n < NUM_SPECIES; ++n) {
                Ze += w_elem(e, n) * in_ma[box_no](i, j, k, n);
              }
              out_ma[box_no](i, j, k, idElemOut + e) = Ze;
            }
          }
        });

      // Copy auxiliary variables unchanged from the input to the output
      for (int ivar = 0; ivar < nAuxVar; ++ivar) {
        MultiFab::Copy(
          outdata[lev], indata, idAuxLocal + ivar, idAuxOut + ivar, 1, nGrow);
      }

      Print() << "Derive finished for level " << lev << std::endl;
    }

    std::string outfile(getFileRoot(plotFileName) + outsuffix);

    // Cap the number of plotfile data files via the n_files option (AMReX)
    int n_files = amrex::VisMF::GetNOutFiles();
    pp.query("n_files", n_files);
    amrex::VisMF::SetNOutFiles(n_files);

    Print() << "Writing new data to " << outfile << std::endl;
    Vector<int> isteps(Nlev, 0);
    Vector<IntVect> refRatios(Nlev - 1);
    for (int lev = 0; lev < Nlev - 1; ++lev) {
      const int rr = amrData.RefRatio()[lev];
      refRatios[lev] = IntVect{AMREX_D_DECL(rr, rr, rr)};
    }
    amrex::WriteMultiLevelPlotfile(
      outfile, Nlev, GetVecOfConstPtrs(outdata), outNames, geoms, 0.0, isteps,
      refRatios);
  }
  Finalize();
  return 0;
}
