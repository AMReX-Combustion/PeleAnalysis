#include <string>
#include <iostream>

#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include <AMReX_DataServices.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_MLMG.H>
#include <AMReX_MLPoisson.H>
#include <AMReX_MLABecLaplacian.H>

using namespace amrex;

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
  std::vector<std::string> tokens = Tokenize(infile,std::string("/"));
  return tokens[tokens.size()-1];
}

int
main (int   argc,
      char* argv[])
{
  amrex::Initialize(argc,argv);
  {
    if (argc < 2) {
      print_usage(argc,argv);
    }

    // ---------------------------------------------------------------------
    // Set defaults input values
    // ---------------------------------------------------------------------
    //std::string gradVar       = "temp";
    std::string infile        = "";  
    int finestLevel           = 1000;
    //int nAuxVar               = 0;

    // ---------------------------------------------------------------------
    // ParmParse
    // ---------------------------------------------------------------------
    ParmParse pp;

    if (pp.contains("help")) {
      print_usage(argc,argv);
    }

    pp.get("infile",infile);
    //pp.query("gradVar",gradVar);
    pp.query("finestLevel",finestLevel);

    // Initialize DataService
    DataServices::SetBatchMode();
    Amrvis::FileType fileType(Amrvis::NEWPLT);
    DataServices dataServices(infile, fileType);
    if( ! dataServices.AmrDataOk()) {
      DataServices::Dispatch(DataServices::ExitRequest, NULL);
    }
    AmrData& amrData = dataServices.AmrDataRef();

    // Plotfile global infos
    finestLevel = std::min(finestLevel,amrData.FinestLevel());
    int Nlev = finestLevel + 1;
    const Vector<std::string>& plotVarNames = amrData.PlotVarNames();
    RealBox rb(&(amrData.ProbLo()[0]), 
               &(amrData.ProbHi()[0]));

    // Gradient variable
    //int idC = -1;
    //for (int i=0; i<plotVarNames.size(); ++i)
    //{
    //  if (plotVarNames[i] == gradVar) idC = i;
    //}
    //if (idC<0) {
    //  Print() << "Cannot find " << gradVar << " data in pltfile \n";
    //}
    int ID_VEL_X, ID_VEL_Y, ID_VEL_Z;
    Vector<int> ID_VEL_VEC = {ID_VEL_X, ID_VEL_Y, ID_VEL_Z};
    for (int i=0; i<plotVarNames.size(); ++i) {
      if (plotVarNames[i] == "x_velocity") ID_VEL_X = i;
      else if (plotVarNames[i] == "y_velocity") ID_VEL_Y = i;
      else if (plotVarNames[i] == "z_velocity") ID_VEL_Z = i;
    }

    // Auxiliary variables
    //nAuxVar = pp.countval("Aux_Variables");
    //Vector<std::string> AuxVar(nAuxVar);
    //for(int ivar = 0; ivar < nAuxVar; ++ivar) { 
    //     pp.get("Aux_Variables", AuxVar[ivar],ivar);
    //}

    // ---------------------------------------------------------------------
    // Variables index management
    // ---------------------------------------------------------------------
    //const int idCst = 0;
    int nCompIn = AMREX_SPACEDIM;
    Vector<std::string> inVarNames(nCompIn);
    inVarNames[0] = plotVarNames[ID_VEL_X];
    inVarNames[1] = plotVarNames[ID_VEL_Y];
    inVarNames[2] = plotVarNames[ID_VEL_Z];

    //if (nAuxVar>0)
    //{
    //    inVarNames.resize(nCompIn+nAuxVar);
    //    for (int ivar=0; ivar<nAuxVar; ++ivar) {
    //        if ( amrData.StateNumber(AuxVar[ivar]) < 0 ) {
    //           amrex::Abort("Unknown auxiliary variable name: "+AuxVar[ivar]);
    //        }
    //        inVarNames[nCompIn] = AuxVar[ivar];
    //        nCompIn ++;
    //    } 
    //}

    Vector<int> destFillComps(nCompIn);
    for (int i=0; i<nCompIn; ++i) {
      destFillComps[i] = i;
    }

    Vector<int> idGr(AMREX_SPACEDIM);
    Print() << "!!! nCompIn " << nCompIn << "\n";
    idGr[0] = nCompIn;
    idGr[1] = nCompIn+3;
#if AMREX_SPACEDIM==3
    idGr[2] = nCompIn+6;
#endif
    //const int nCompOut = idGr + AMREX_SPACEDIM +1 ; // 1 component stores the ||gradT||
    const int nCompOut = nCompIn + AMREX_SPACEDIM*AMREX_SPACEDIM;

    // Check symmetry/periodicity in given coordinate direction
    Vector<int> sym_dir(AMREX_SPACEDIM,0);
    pp.queryarr("sym_dir",sym_dir,0,AMREX_SPACEDIM);  

    Vector<int> is_per(AMREX_SPACEDIM,1);
    pp.queryarr("is_per",is_per,0,AMREX_SPACEDIM);
    Print() << "Periodicity assumed for this case: ";
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        Print() << is_per[idim] << " ";
    }
    Print() << "\n";
    BCRec gradVarBC;
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        gradVarBC.setLo(idim,BCType::foextrap);
        gradVarBC.setHi(idim,BCType::foextrap);
        if ( is_per[idim] ) {
            gradVarBC.setLo(idim, BCType::int_dir);
            gradVarBC.setHi(idim, BCType::int_dir);
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
    for (int jdim=0; jdim<AMREX_SPACEDIM; jdim++) {
    for (int lev=0; lev<Nlev; ++lev) {

      const BoxArray ba = amrData.boxArray(lev);
      grids[lev] = ba;
      dmap[lev] = DistributionMapping(ba);
      geoms[lev] = Geometry(amrData.ProbDomain()[lev],&rb,coord,&(is_per[0]));
      state[lev].define(grids[lev], dmap[lev], nCompOut, nGrow);

      Print() << "Reading data for level: " << lev << std::endl;
      amrData.FillVar(state[lev], lev, inVarNames, destFillComps);

      state[lev].FillBoundary(ID_VEL_VEC[jdim],1,geoms[lev].periodicity());
    }
    }

    // Get face-centered gradients from MLMG
    LPInfo info;
    info.setAgglomeration(1);
    info.setConsolidation(1);
    info.setMetricTerm(false);
    info.setMaxCoarseningLevel(0);
    MLPoisson poisson({geoms}, {grids}, {dmap}, info);
    poisson.setMaxOrder(4);
    std::array<LinOpBCType, AMREX_SPACEDIM> lo_bc;
    std::array<LinOpBCType, AMREX_SPACEDIM> hi_bc;
    for (int idim = 0; idim< AMREX_SPACEDIM; idim++){
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
    poisson.setDomainBC(lo_bc, hi_bc);

// Need to apply the operator to ensure CF consistency with composite solve
    int nGrowGrad = 0; // No need for ghost face on gradient

    Vector<Array<MultiFab,AMREX_SPACEDIM>> grad_x_vel(Nlev);
    Vector<Array<MultiFab,AMREX_SPACEDIM>> grad_y_vel(Nlev);
    Vector<Array<MultiFab,AMREX_SPACEDIM>> grad_z_vel(Nlev);
    std::array<Vector<Array<MultiFab,AMREX_SPACEDIM>>*, 3> grad_vel_ptrs = {
        &grad_x_vel,
        &grad_y_vel,
        &grad_z_vel
    };

    for (int jdim = 0; jdim < AMREX_SPACEDIM; jdim++) {
      Vector<std::unique_ptr<MultiFab>> phi;
      Vector<MultiFab> laps;
      for (int lev = 0; lev < Nlev; ++lev) {
        for (int idim = 0; idim < AMREX_SPACEDIM; idim++) {
          const auto& ba = grids[lev];
          (*grad_vel_ptrs[jdim])[lev][idim].define(amrex::convert(ba,IntVect::TheDimensionVector(idim)),
                                  dmap[lev], 1, nGrowGrad);
        }
        phi.push_back(std::make_unique<MultiFab> (state[lev],amrex::make_alias,ID_VEL_VEC[jdim],1));
        poisson.setLevelBC(lev, phi[lev].get());
        laps.emplace_back(grids[lev], dmap[lev], 1, 1);
      }

      MLMG mlmg(poisson);
      mlmg.apply(GetVecOfPtrs(laps), GetVecOfPtrs(phi));
      mlmg.getFluxes(GetVecOfArrOfPtrs((*grad_vel_ptrs[jdim])), GetVecOfPtrs(phi), MLMG::Location::FaceCenter);

      for (int lev = 0; lev < Nlev; ++lev) {
        // Convert to cell avg gradient
        MultiFab gradAlias(state[lev], amrex::make_alias, idGr[jdim], AMREX_SPACEDIM);
        average_face_to_cellcenter(gradAlias, 0, GetArrOfConstPtrs((*grad_vel_ptrs[jdim])[lev]));
        gradAlias.mult(-1.0);
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
       // for (MFIter mfi(state[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
       // {    
       //   const Box& bx = mfi.tilebox();
       //   auto const& grad_a   = gradAlias.const_array(mfi);
       //   auto const& gradMag  = state[lev].array(mfi,idGr[0]+AMREX_SPACEDIM*AMREX_SPACEDIM);
       //   amrex::ParallelFor(bx, [=]
       //   AMREX_GPU_DEVICE (int i, int j, int k) noexcept
       //   {    
       //      gradMag(i,j,k) = std::sqrt(AMREX_D_TERM(  grad_a(i,j,k,0) * grad_a(i,j,k,0),
       //                                              + grad_a(i,j,k,1) * grad_a(i,j,k,1),
       //                                              + grad_a(i,j,k,2) * grad_a(i,j,k,2)));
       //   });  
       // }
      }
    }

    // ---------------------------------------------------------------------
    // Write the results
    // ---------------------------------------------------------------------
    Vector<std::string> nnames(nCompOut);
    Print() << "!!! nCompOut = " << nCompOut << "\n";
    for (int i=0; i<nCompIn; ++i) {
      nnames[i] = inVarNames[i];
    }
    for (int jdim=0; jdim<AMREX_SPACEDIM; jdim++) {
      nnames[idGr[jdim] + 0] = inVarNames[jdim] + "_gx";
      nnames[idGr[jdim] + 1] = inVarNames[jdim] + "_gy";
#if AMREX_SPACEDIM==3
      nnames[idGr[jdim] + 2] = inVarNames[jdim] + "_gz";
#endif
    }
    for (int i=0; i<nnames.size(); i++) {
      Print() << "!!!nnames(" << i << "): " << nnames[i] << "\n";
    }
//    nnames[idGr+0] = "x_velocity_gx";
//    nnames[idGr+1] = "y_velocity_gy";
//#if AMREX_SPACEDIM==3
//    nnames[idGr+2] = "z_velocity_gz";
//#endif
    //nnames[idGr[2]+AMREX_SPACEDIM] = "test";
    std::string outfile(getFileRoot(infile) + "_gt"); pp.query("outfile",outfile);

    Print() << "Writing new data to " << outfile << std::endl;
    Vector<int> isteps(Nlev, 0);
    Vector<IntVect> refRatios(Nlev-1,{AMREX_D_DECL(2, 2, 2)});
    amrex::WriteMultiLevelPlotfile(outfile, Nlev, GetVecOfConstPtrs(state), nnames,
                                   geoms, 0.0, isteps, refRatios);
  }
  amrex::Finalize();
  return 0;
}
