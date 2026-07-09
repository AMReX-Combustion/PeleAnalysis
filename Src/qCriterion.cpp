#include <iostream>

#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include <AMReX_DataServices.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_VisMF.H>
#include <AMReX_MLMG.H>
#include <AMReX_MLPoisson.H>
#include <AMReX_MLABecLaplacian.H>

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
  amrex::Initialize(argc, argv);
  {
    static_assert(
      AMREX_SPACEDIM == 3, "This tool must be compiled with DIM = 3");

    if (argc < 2) {
      print_usage(argc, argv);
    } else if (
      (std::strcmp(argv[1], "-h") == 0) ||
      (std::strcmp(argv[1], "--help") == 0)) {
      print_usage(argc, argv);
    }

    // ---------------------------------------------------------------------
    // Set defaults input values
    // ---------------------------------------------------------------------
    std::string infile = "";
    int finestLevel = 1000;
    int nAuxVar = 0;

    // ---------------------------------------------------------------------
    // ParmParse
    // ---------------------------------------------------------------------
    ParmParse pp;

    pp.get("infile", infile);
    pp.query("finestLevel", finestLevel);

    // Initialize DataService
    DataServices::SetBatchMode();
    Amrvis::FileType fileType(Amrvis::NEWPLT);
    DataServices dataServices(infile, fileType);
    if (!dataServices.AmrDataOk()) {
      DataServices::Dispatch(DataServices::ExitRequest, NULL);
    }
    AmrData& amrData = dataServices.AmrDataRef();

    // Plotfile global infos
    finestLevel = std::min(finestLevel, amrData.FinestLevel());
    int Nlev = finestLevel + 1;
    const Vector<std::string>& plotVarNames = amrData.PlotVarNames();
    RealBox rb(&(amrData.ProbLo()[0]), &(amrData.ProbHi()[0]));

    int ID_VEL_X = -1, ID_VEL_Y = -1, ID_VEL_Z = -1;
    std::string gradVar_vel_x = "x_velocity";
    std::string gradVar_vel_y = "y_velocity";
    std::string gradVar_vel_z = "z_velocity";
    for (int i = 0; i < plotVarNames.size(); ++i) {
      if (plotVarNames[i] == gradVar_vel_x)
        ID_VEL_X = i;
      else if (plotVarNames[i] == gradVar_vel_y)
        ID_VEL_Y = i;
      else if (plotVarNames[i] == gradVar_vel_z)
        ID_VEL_Z = i;
    }
    if (ID_VEL_X == -1 || ID_VEL_Y == -1 || ID_VEL_Z == -1) {
      std::cerr
        << "At least one velocity component was not found in the plt-file!\n";
    }

    // Auxiliary variables
    nAuxVar = pp.countval("Aux_Variables");
    Vector<std::string> AuxVar(nAuxVar);
    for (int ivar = 0; ivar < nAuxVar; ++ivar) {
      pp.get("Aux_Variables", AuxVar[ivar], ivar);
    }

    // ---------------------------------------------------------------------
    // Variables index management
    // ---------------------------------------------------------------------
    int nCompIn = AMREX_SPACEDIM /*velocities*/;
    Vector<std::string> inVarNames(nCompIn);
    inVarNames[0] = plotVarNames[ID_VEL_X];
    inVarNames[1] = plotVarNames[ID_VEL_Y];
    inVarNames[2] = plotVarNames[ID_VEL_Z];

    if (nAuxVar > 0) {
      for (int ivar = 0; ivar < nAuxVar; ++ivar) {
        if (amrData.StateNumber(AuxVar[ivar]) < 0) {
          amrex::Abort("Unknown auxiliary variable name: " + AuxVar[ivar]);
        }
        inVarNames.push_back(AuxVar[ivar]);
      }
      nCompIn += nAuxVar;
    }

    Vector<int> destFillComps(nCompIn);
    for (int i = 0; i < nCompIn; ++i) {
      destFillComps[i] = i;
    }

    const int idGr_vel_x = nCompIn + 0 * AMREX_SPACEDIM;
    const int idGr_vel_y = nCompIn + 1 * AMREX_SPACEDIM;
    const int idGr_vel_z = nCompIn + 2 * AMREX_SPACEDIM;
    const int idQ = nCompIn + 3 * AMREX_SPACEDIM;
    const int idQNorm = nCompIn + 3 * AMREX_SPACEDIM + 1 /*Q*/;
    const int nCompOut = nCompIn + 3 * AMREX_SPACEDIM + 2 /*Q and Q_norm*/;

    // Check symmetry/periodicity in given coordinate direction
    Vector<int> sym_dir(AMREX_SPACEDIM, 0);
    pp.queryarr("sym_dir", sym_dir, 0, AMREX_SPACEDIM);

    Vector<int> is_per(AMREX_SPACEDIM, 1);
    pp.queryarr("is_per", is_per, 0, AMREX_SPACEDIM);
    Print() << "Periodicity assumed for this case: ";
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
      Print() << is_per[idim] << " ";
    }
    Print() << "\n";
    BCRec vel_x_BC;
    BCRec vel_y_BC;
    BCRec vel_z_BC;
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
      vel_x_BC.setLo(idim, BCType::foextrap);
      vel_y_BC.setLo(idim, BCType::foextrap);
      vel_z_BC.setLo(idim, BCType::foextrap);
      vel_x_BC.setHi(idim, BCType::foextrap);
      vel_y_BC.setHi(idim, BCType::foextrap);
      vel_z_BC.setHi(idim, BCType::foextrap);
      if (is_per[idim]) {
        vel_x_BC.setLo(idim, BCType::int_dir);
        vel_y_BC.setLo(idim, BCType::int_dir);
        vel_z_BC.setLo(idim, BCType::int_dir);
        vel_x_BC.setHi(idim, BCType::int_dir);
        vel_y_BC.setHi(idim, BCType::int_dir);
        vel_z_BC.setHi(idim, BCType::int_dir);
      }
    }

    int coord = 0;

    // ---------------------------------------------------------------------
    // Let's start the real work
    // ---------------------------------------------------------------------
    Vector<MultiFab> state(Nlev);
    Vector<Geometry> geoms(Nlev);
    Vector<BoxArray> grids(Nlev);
    Vector<DistributionMapping> dmap(Nlev);
    const int nGrow = 1;

    // Read data on all the levels
    for (int lev = 0; lev < Nlev; ++lev) {

      const BoxArray ba = amrData.boxArray(lev);
      grids[lev] = ba;
      dmap[lev] = DistributionMapping(ba);
      geoms[lev] =
        Geometry(amrData.ProbDomain()[lev], &rb, coord, &(is_per[0]));
      state[lev].define(grids[lev], dmap[lev], nCompOut, nGrow);

      Print() << "Reading data for level: " << lev << std::endl;
      amrData.FillVar(state[lev], lev, inVarNames, destFillComps);

      state[lev].FillBoundary(ID_VEL_X, 1, geoms[lev].periodicity());
      state[lev].FillBoundary(ID_VEL_Y, 1, geoms[lev].periodicity());
      state[lev].FillBoundary(ID_VEL_Z, 1, geoms[lev].periodicity());
    }

    // Get face-centered gradients from MLMG
    LPInfo info;
    info.setAgglomeration(1);
    info.setConsolidation(1);
    info.setMetricTerm(false);
    info.setMaxCoarseningLevel(0);
    MLPoisson poisson_vel_x({geoms}, {grids}, {dmap}, info);
    MLPoisson poisson_vel_y({geoms}, {grids}, {dmap}, info);
    MLPoisson poisson_vel_z({geoms}, {grids}, {dmap}, info);
    poisson_vel_x.setMaxOrder(4);
    poisson_vel_y.setMaxOrder(4);
    poisson_vel_z.setMaxOrder(4);
    std::array<LinOpBCType, AMREX_SPACEDIM> lo_bc;
    std::array<LinOpBCType, AMREX_SPACEDIM> hi_bc;
    for (int idim = 0; idim < AMREX_SPACEDIM; idim++) {
      if (is_per[idim] == 1) {
        lo_bc[idim] = hi_bc[idim] = LinOpBCType::Periodic;
      } else {
        if (sym_dir[idim] == 1) {
          lo_bc[idim] = hi_bc[idim] = LinOpBCType::reflect_odd;
        } else {
          lo_bc[idim] = hi_bc[idim] = LinOpBCType::Neumann;
        }
      }
    }
    poisson_vel_x.setDomainBC(lo_bc, hi_bc);
    poisson_vel_y.setDomainBC(lo_bc, hi_bc);
    poisson_vel_z.setDomainBC(lo_bc, hi_bc);

    // Need to apply the operator to ensure CF consistency with composite solve
    int nGrowGrad = 0; // No need for ghost face on gradient
    Vector<Array<MultiFab, AMREX_SPACEDIM>> grad_vel_x(Nlev);
    Vector<Array<MultiFab, AMREX_SPACEDIM>> grad_vel_y(Nlev);
    Vector<Array<MultiFab, AMREX_SPACEDIM>> grad_vel_z(Nlev);
    Vector<std::unique_ptr<MultiFab>> phi_vel_x;
    Vector<std::unique_ptr<MultiFab>> phi_vel_y;
    Vector<std::unique_ptr<MultiFab>> phi_vel_z;
    Vector<MultiFab> laps_vel_x;
    Vector<MultiFab> laps_vel_y;
    Vector<MultiFab> laps_vel_z;
    for (int lev = 0; lev < Nlev; ++lev) {
      for (int idim = 0; idim < AMREX_SPACEDIM; idim++) {
        const auto& ba = grids[lev];
        grad_vel_x[lev][idim].define(
          amrex::convert(ba, IntVect::TheDimensionVector(idim)), dmap[lev], 1,
          nGrowGrad);
        grad_vel_y[lev][idim].define(
          amrex::convert(ba, IntVect::TheDimensionVector(idim)), dmap[lev], 1,
          nGrowGrad);
        grad_vel_z[lev][idim].define(
          amrex::convert(ba, IntVect::TheDimensionVector(idim)), dmap[lev], 1,
          nGrowGrad);
      }
      phi_vel_x.push_back(
        std::make_unique<MultiFab>(state[lev], amrex::make_alias, ID_VEL_X, 1));
      phi_vel_y.push_back(
        std::make_unique<MultiFab>(state[lev], amrex::make_alias, ID_VEL_Y, 1));
      phi_vel_z.push_back(
        std::make_unique<MultiFab>(state[lev], amrex::make_alias, ID_VEL_Z, 1));
      poisson_vel_x.setLevelBC(lev, phi_vel_x[lev].get());
      poisson_vel_y.setLevelBC(lev, phi_vel_y[lev].get());
      poisson_vel_z.setLevelBC(lev, phi_vel_z[lev].get());
      laps_vel_x.emplace_back(grids[lev], dmap[lev], 1, 1);
      laps_vel_y.emplace_back(grids[lev], dmap[lev], 1, 1);
      laps_vel_z.emplace_back(grids[lev], dmap[lev], 1, 1);
    }

    MLMG mlmg_vel_x(poisson_vel_x);
    MLMG mlmg_vel_y(poisson_vel_y);
    MLMG mlmg_vel_z(poisson_vel_z);
    mlmg_vel_x.apply(GetVecOfPtrs(laps_vel_x), GetVecOfPtrs(phi_vel_x));
    mlmg_vel_y.apply(GetVecOfPtrs(laps_vel_y), GetVecOfPtrs(phi_vel_y));
    mlmg_vel_z.apply(GetVecOfPtrs(laps_vel_z), GetVecOfPtrs(phi_vel_z));
    mlmg_vel_x.getFluxes(
      GetVecOfArrOfPtrs(grad_vel_x), GetVecOfPtrs(phi_vel_x),
      MLMG::Location::FaceCenter);
    mlmg_vel_y.getFluxes(
      GetVecOfArrOfPtrs(grad_vel_y), GetVecOfPtrs(phi_vel_y),
      MLMG::Location::FaceCenter);
    mlmg_vel_z.getFluxes(
      GetVecOfArrOfPtrs(grad_vel_z), GetVecOfPtrs(phi_vel_z),
      MLMG::Location::FaceCenter);

    for (int lev = 0; lev < Nlev; ++lev) {
      // Convert to cell avg gradient
      MultiFab gradAlias_vel_x(
        state[lev], amrex::make_alias, idGr_vel_x, AMREX_SPACEDIM);
      MultiFab gradAlias_vel_y(
        state[lev], amrex::make_alias, idGr_vel_y, AMREX_SPACEDIM);
      MultiFab gradAlias_vel_z(
        state[lev], amrex::make_alias, idGr_vel_z, AMREX_SPACEDIM);
      average_face_to_cellcenter(
        gradAlias_vel_x, 0, GetArrOfConstPtrs(grad_vel_x[lev]));
      average_face_to_cellcenter(
        gradAlias_vel_y, 0, GetArrOfConstPtrs(grad_vel_y[lev]));
      average_face_to_cellcenter(
        gradAlias_vel_z, 0, GetArrOfConstPtrs(grad_vel_z[lev]));
      gradAlias_vel_x.mult(-1.0);
      gradAlias_vel_y.mult(-1.0);
      gradAlias_vel_z.mult(-1.0);
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (MFIter mfi(state[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();
        auto const& grad_vel_x_a = gradAlias_vel_x.const_array(mfi);
        auto const& grad_vel_y_a = gradAlias_vel_y.const_array(mfi);
        auto const& grad_vel_z_a = gradAlias_vel_z.const_array(mfi);

        // Intermediates:
        MultiFab S_mf(
          state[lev].boxArray(), state[lev].DistributionMap(), 9, 0);
        auto const& S_a = S_mf[mfi].array();
        MultiFab S_abs_mf(
          state[lev].boxArray(), state[lev].DistributionMap(), 1, 0);
        auto const& S_abs_a = S_abs_mf[mfi].array();
        MultiFab Omega_mf(
          state[lev].boxArray(), state[lev].DistributionMap(), 9, 0);
        auto const& Omega_a = Omega_mf[mfi].array();
        MultiFab Omega_abs_mf(
          state[lev].boxArray(), state[lev].DistributionMap(), 1, 0);
        auto const& Omega_abs_a = Omega_abs_mf[mfi].array();

        auto const& Q = state[lev].array(mfi, idQ);
        auto const& Q_norm = state[lev].array(mfi, idQNorm);

        amrex::ParallelFor(
          bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            // Compute Strain Rate Tensor S_ij:
            S_a(i, j, k, AMREX_SPACEDIM * 0 + 0) = grad_vel_x_a(i, j, k, 0);
            S_a(i, j, k, AMREX_SPACEDIM * 0 + 1) =
              0.5 * (grad_vel_x_a(i, j, k, 1) + grad_vel_y_a(i, j, k, 0));
            S_a(i, j, k, AMREX_SPACEDIM * 0 + 2) =
              0.5 * (grad_vel_x_a(i, j, k, 2) + grad_vel_z_a(i, j, k, 0));

            S_a(i, j, k, AMREX_SPACEDIM * 1 + 0) =
              0.5 * (grad_vel_y_a(i, j, k, 0) + grad_vel_x_a(i, j, k, 1));
            S_a(i, j, k, AMREX_SPACEDIM * 1 + 1) = grad_vel_y_a(i, j, k, 1);
            S_a(i, j, k, AMREX_SPACEDIM * 1 + 2) =
              0.5 * (grad_vel_y_a(i, j, k, 2) + grad_vel_z_a(i, j, k, 1));

            S_a(i, j, k, AMREX_SPACEDIM * 2 + 0) =
              0.5 * (grad_vel_z_a(i, j, k, 0) + grad_vel_x_a(i, j, k, 2));
            S_a(i, j, k, AMREX_SPACEDIM * 2 + 1) =
              0.5 * (grad_vel_z_a(i, j, k, 1) + grad_vel_y_a(i, j, k, 2));
            S_a(i, j, k, AMREX_SPACEDIM * 2 + 2) = grad_vel_z_a(i, j, k, 2);

            S_abs_a(i, j, k) = S_a(i, j, k, AMREX_SPACEDIM * 0 + 0) *
                                 S_a(i, j, k, AMREX_SPACEDIM * 0 + 0) +
                               S_a(i, j, k, AMREX_SPACEDIM * 0 + 1) *
                                 S_a(i, j, k, AMREX_SPACEDIM * 0 + 1) +
                               S_a(i, j, k, AMREX_SPACEDIM * 0 + 2) *
                                 S_a(i, j, k, AMREX_SPACEDIM * 0 + 2) +
                               S_a(i, j, k, AMREX_SPACEDIM * 1 + 0) *
                                 S_a(i, j, k, AMREX_SPACEDIM * 1 + 0) +
                               S_a(i, j, k, AMREX_SPACEDIM * 1 + 1) *
                                 S_a(i, j, k, AMREX_SPACEDIM * 1 + 1) +
                               S_a(i, j, k, AMREX_SPACEDIM * 1 + 2) *
                                 S_a(i, j, k, AMREX_SPACEDIM * 1 + 2) +
                               S_a(i, j, k, AMREX_SPACEDIM * 2 + 0) *
                                 S_a(i, j, k, AMREX_SPACEDIM * 2 + 0) +
                               S_a(i, j, k, AMREX_SPACEDIM * 2 + 1) *
                                 S_a(i, j, k, AMREX_SPACEDIM * 2 + 1) +
                               S_a(i, j, k, AMREX_SPACEDIM * 2 + 2) *
                                 S_a(i, j, k, AMREX_SPACEDIM * 2 + 2);

            // Compute Vorticity Tensor Omega_ij:
            Omega_a(i, j, k, AMREX_SPACEDIM * 0 + 0) = 0.0;
            Omega_a(i, j, k, AMREX_SPACEDIM * 0 + 1) =
              0.5 * (grad_vel_x_a(i, j, k, 1) - grad_vel_y_a(i, j, k, 0));
            Omega_a(i, j, k, AMREX_SPACEDIM * 0 + 2) =
              0.5 * (grad_vel_x_a(i, j, k, 2) - grad_vel_z_a(i, j, k, 0));

            Omega_a(i, j, k, AMREX_SPACEDIM * 1 + 0) =
              0.5 * (grad_vel_y_a(i, j, k, 0) - grad_vel_x_a(i, j, k, 1));
            Omega_a(i, j, k, AMREX_SPACEDIM * 1 + 1) = 0.0;
            Omega_a(i, j, k, AMREX_SPACEDIM * 1 + 2) =
              0.5 * (grad_vel_y_a(i, j, k, 2) - grad_vel_z_a(i, j, k, 1));

            Omega_a(i, j, k, AMREX_SPACEDIM * 2 + 0) =
              0.5 * (grad_vel_z_a(i, j, k, 0) - grad_vel_x_a(i, j, k, 2));
            Omega_a(i, j, k, AMREX_SPACEDIM * 2 + 1) =
              0.5 * (grad_vel_z_a(i, j, k, 1) - grad_vel_y_a(i, j, k, 2));
            Omega_a(i, j, k, AMREX_SPACEDIM * 2 + 2) = 0.0;

            Omega_abs_a(i, j, k) = Omega_a(i, j, k, AMREX_SPACEDIM * 0 + 0) *
                                     Omega_a(i, j, k, AMREX_SPACEDIM * 0 + 0) +
                                   Omega_a(i, j, k, AMREX_SPACEDIM * 0 + 1) *
                                     Omega_a(i, j, k, AMREX_SPACEDIM * 0 + 1) +
                                   Omega_a(i, j, k, AMREX_SPACEDIM * 0 + 2) *
                                     Omega_a(i, j, k, AMREX_SPACEDIM * 0 + 2) +
                                   Omega_a(i, j, k, AMREX_SPACEDIM * 1 + 0) *
                                     Omega_a(i, j, k, AMREX_SPACEDIM * 1 + 0) +
                                   Omega_a(i, j, k, AMREX_SPACEDIM * 1 + 1) *
                                     Omega_a(i, j, k, AMREX_SPACEDIM * 1 + 1) +
                                   Omega_a(i, j, k, AMREX_SPACEDIM * 1 + 2) *
                                     Omega_a(i, j, k, AMREX_SPACEDIM * 1 + 2) +
                                   Omega_a(i, j, k, AMREX_SPACEDIM * 2 + 0) *
                                     Omega_a(i, j, k, AMREX_SPACEDIM * 2 + 0) +
                                   Omega_a(i, j, k, AMREX_SPACEDIM * 2 + 1) *
                                     Omega_a(i, j, k, AMREX_SPACEDIM * 2 + 1) +
                                   Omega_a(i, j, k, AMREX_SPACEDIM * 2 + 2) *
                                     Omega_a(i, j, k, AMREX_SPACEDIM * 2 + 2);

            // Compute Q-criterion:
            Q(i, j, k) = 0.5 * (Omega_abs_a(i, j, k) - S_abs_a(i, j, k));
            Q_norm(i, j, k) = (Omega_abs_a(i, j, k) - S_abs_a(i, j, k)) /
                              (Omega_abs_a(i, j, k) + S_abs_a(i, j, k));
          });
      }
    }

    // ---------------------------------------------------------------------
    // Write the results
    // ---------------------------------------------------------------------
    Vector<std::string> nnames(nCompOut);
    for (int i = 0; i < nCompIn; ++i) {
      nnames[i] = inVarNames[i];
    }
    nnames[idGr_vel_x + 0] = gradVar_vel_x + "_gx";
    nnames[idGr_vel_x + 1] = gradVar_vel_x + "_gy";
    nnames[idGr_vel_x + 2] = gradVar_vel_x + "_gz";
    nnames[idGr_vel_y + 0] = gradVar_vel_y + "_gx";
    nnames[idGr_vel_y + 1] = gradVar_vel_y + "_gy";
    nnames[idGr_vel_y + 2] = gradVar_vel_y + "_gz";
    nnames[idGr_vel_z + 0] = gradVar_vel_z + "_gx";
    nnames[idGr_vel_z + 1] = gradVar_vel_z + "_gy";
    nnames[idGr_vel_z + 2] = gradVar_vel_z + "_gz";

    nnames[idQ] = "Q";
    nnames[idQNorm] = "Q_norm";

    std::string outfile(getFileRoot(infile) + "_qCriterion");
    pp.query("outfile", outfile);

    // Cap the number of plotfile data files via the n_files option (AMReX)
    int n_files = amrex::VisMF::GetNOutFiles();
    pp.query("n_files", n_files);
    amrex::VisMF::SetNOutFiles(n_files);

    Print() << "Writing new data to " << outfile << std::endl;
    Vector<int> isteps(Nlev, 0);
    Vector<IntVect> refRatios(Nlev - 1, {AMREX_D_DECL(2, 2, 2)});
    amrex::WriteMultiLevelPlotfile(
      outfile, Nlev, GetVecOfConstPtrs(state), nnames, geoms, 0.0, isteps,
      refRatios);
  }
  amrex::Finalize();
  return 0;
}
