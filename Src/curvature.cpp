#include <string>
#include <iostream>
#include <set>

#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MultiFabUtil_C.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_DataServices.H>
#include <AMReX_BCRec.H>
#include <AMReX_Interpolater.H>
#include <AppendToPlotFile.H>
#include <WritePlotFile.H>
#include <AMReX_VisMF.H>

#include <AMReX_MLMG.H>
#include <AMReX_MLPoisson.H>
#include <AMReX_MLABecLaplacian.H>
#include <AMReX_BLFort.H>

using namespace amrex;

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
    // if (argc < 2)
    //    print_usage(argc,argv);
    
    
    // if (pp.contains("help"))
    //    print_usage(argc,argv);

    // ---------------------------------------------------------------------
    // Set defaults input values
    // ---------------------------------------------------------------------
    std::string progressName  = "temp";
    Real progMin              = 1.0e20;
    Real progMax              = -1.0e20;
    int finestLevel           = 1000;
    int verbose               = 0;
    int floorIt               = 0;
    int useFileMinMax         = 1;
    bool do_strain            = false;
    bool do_smooth            = false; 
    bool do_gaussCurv         = false; 
    bool getStrainTensor      = false;
    bool do_velnormal         = false; 
    bool getProgGrad          = false;
    Real smooth_time          = 1.0e-7;
    int nAuxVar               = 0;  
    bool appendPlotFile       = false;

    
    // ---------------------------------------------------------------------
    // ParmParse
    // ---------------------------------------------------------------------
    ParmParse pp;

    // IO
    pp.query("verbose",verbose);
    std::string plotFileName;
    pp.get("infile",plotFileName);
    std::string outfile(getFileRoot(plotFileName) + "_K");
    pp.query("outfile",outfile);
    pp.query("finestLevel",finestLevel);
    pp.query("appendPlotFile",appendPlotFile);
    pp.query("do_gaussCurv",do_gaussCurv);

    // Progress variable
    pp.query("progressName",progressName);
    pp.query("progMin",progMin);
    pp.query("progMax",progMax);
    pp.query("floorIt",floorIt);
    pp.query("useFileMinMax",useFileMinMax);
    pp.query("getProgGrad",getProgGrad);

    // Progress variable smoothing
    pp.query("do_smooth",do_smooth);
    pp.query("smoothing_time",smooth_time);
    
    // Pertaining to the velocity computation
    pp.query("do_strain",do_strain);
    if (do_strain) {
       pp.query("getStrainTensor",getStrainTensor);
    }
    pp.query("do_velnormal",do_velnormal);  

    // Auxiliary variables
    nAuxVar = pp.countval("Aux_Variables");
    Vector<string> AuxVar(nAuxVar);
    for(int ivar = 0; ivar < nAuxVar; ++ivar) { 
         pp.get("Aux_Variables", AuxVar[ivar],ivar);
    }

    if (verbose>1) AmrData::SetVerbose(true);
    
    Print() << "infile = " << plotFileName << "\n";
    Print() << "reading plt file = " << plotFileName << "\n";
    
    // Initialize DataService
    DataServices::SetBatchMode();
    Amrvis::FileType fileType(Amrvis::NEWPLT);
    DataServices dataServices(plotFileName, fileType);
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

    // ---------------------------------------------------------------------
    // Progress variable construction
    // ---------------------------------------------------------------------
    int idC = amrData.StateNumber(progressName);
    if ( idC < 0 ) {
         amrex::Abort("Wrong progress variable name: "+progressName);
    }
    Real progMinlvl = 1.0e20;
    Real progMaxlvl = -1.0e20;
    if (useFileMinMax || floorIt)
    {
        if (useFileMinMax) {
            for (int lev=0; lev<Nlev; ++lev) {
                amrData.MinMax(amrData.ProbDomain()[lev], progressName, lev, progMinlvl, progMaxlvl);
                progMin = std::min(progMin,progMinlvl);
                progMax = std::max(progMax,progMaxlvl);
            }
            ParallelDescriptor::ReduceRealMin(progMin);
            ParallelDescriptor::ReduceRealMax(progMax);
        }

        Print() << "progressName = " << progressName << " at index: " << amrData.StateNumber(progressName) << "\n";
        Print() << "useFileMinMax = " << useFileMinMax << "\n";
        Print() << "Min/Max = " << progMin << " / " << progMax << "\n";

        ParallelDescriptor::Barrier();

        if (progMin >= progMax) {
            amrex::Abort("progMin must be less than progMax");
        }
    }

    // ---------------------------------------------------------------------
    // Variables index management
    // ---------------------------------------------------------------------
    const int idCst = 0;
    const int idVst = idCst + 1;
    int nCompIn = idVst;

    Vector<std::string> inVarNames(nCompIn);
    inVarNames[idCst] = plotVarNames[idC];

    if (do_strain) {
      nCompIn += AMREX_SPACEDIM;
      inVarNames.resize(nCompIn);
      inVarNames[idVst+0] = "x_velocity";
      inVarNames[idVst+1] = "y_velocity";
#if AMREX_SPACEDIM == 3
      inVarNames[idVst+2] = "z_velocity";
#endif
    }

    if (nAuxVar>0)
    {
        inVarNames.resize(nCompIn+nAuxVar);
        for (int ivar=0; ivar<nAuxVar; ++ivar) {
            if ( amrData.StateNumber(AuxVar[ivar]) < 0 ) {
               amrex::Abort("Unknown auxiliary variable name: "+AuxVar[ivar]);
            }
            inVarNames[nCompIn] = AuxVar[ivar];
            nCompIn ++;
        } 
    }

    Vector<int> destFillComps(nCompIn);
    for (int i=0; i<nCompIn; ++i) {
        destFillComps[i] = i;    
    }

    // Start appending variables that we will compute
    const int idProg = nCompIn;              // progress variable
    const int idSmProg = idProg + 1;         // smoothed progress variable
    const int idKm = idSmProg + 1;           // Mean curvature
    const int idN= idKm + 1;                 // Flame normal components
    int idSR = -1;                           // Strain rate
    int nCompOut = 0;                        // Total number of output variables

#if AMREX_SPACEDIM > 2
    const int idKg = idN + AMREX_SPACEDIM;   // Gaussian curvature
    if (do_strain) {
      idSR = idKg + 1;
      nCompOut = idSR + 1;
    }
    else {
      nCompOut = idKg + 1;
    }
#else
    if (do_strain) {
      idSR = idN + AMREX_SPACEDIM;
      nCompOut = idSR + 1;
    }
    else {
      nCompOut = idN + AMREX_SPACEDIM;
    }
#endif

    int idROST=-1;
    if (getStrainTensor) {
        idROST = nCompOut;
        nCompOut = idROST + 9; // Rate-of-strain, always return 3D set
    }

    int idVelNormal=-1;
    if ( do_velnormal ) {  
        idVelNormal = nCompOut; 
        nCompOut += AMREX_SPACEDIM; 
    }  

    int idProgGrad=-1;
    if (getProgGrad)
    {
        idProgGrad = nCompOut;
        nCompOut = idProgGrad + 3; // Always return 3D set
    }

    if (verbose) {  
       Print() << "Will read the following states: ";
       for (int i=0; i<nCompIn; ++i) {
           Print() << " " << amrData.StateNumber(inVarNames[i]);
       }
       Print() << '\n';
       Print() << "States out will be those plus: " << '\n';
       Print() << "   idProg:   " << idProg << '\n';
       Print() << "   idSmProg: " << idSmProg << '\n';
       Print() << "   idKm:     " << idKm << '\n';
#if AMREX_SPACEDIM > 2
       Print() << "   idKg:     " << idKg << '\n';
#endif
       if (do_strain) {
          Print() << "   idSR:     " << idSR << '\n';
       }
    }
    
    // Check symmetry/periodicity in given coordinate direction
    Vector<int> sym_dir(AMREX_SPACEDIM,0);
    pp.queryarr("sym_dir",sym_dir,0,AMREX_SPACEDIM);  

    Vector<int> is_per(AMREX_SPACEDIM,1);
    pp.queryarr("is_per",is_per,0,AMREX_SPACEDIM);

    Print() << "Periodicity assumed for this case: ";
    for (int i=0; i<AMREX_SPACEDIM; ++i) {
        Print() << is_per[i] << " ";
    }
    Print() << "\n";

    int coord = 0;

    // ---------------------------------------------------------------------
    // Let's start the real work
    // ---------------------------------------------------------------------
    Vector<MultiFab*> state(Nlev);
    Vector<MultiFab*> flame_normal(Nlev);
    Vector<MultiFab*> cell_normal(Nlev);
    Vector<Geometry*> geoms(Nlev);
    Vector<Geometry>  geomsOP(Nlev);
    Vector<BoxArray>  grids(Nlev);
    Vector<DistributionMapping> dmap(Nlev);
    const int nGrow = 2 ;

    // Read the required data from pltfile and compute progress variable while we're at it.
    FArrayBox tmp;
    for (int lev=0; lev<Nlev; ++lev)
    {
        const BoxArray ba = amrData.boxArray(lev);
        grids[lev] = ba;
        DistributionMapping dm(ba);
        dmap[lev] = dm; 
        geoms[lev] = new Geometry(amrData.ProbDomain()[lev],&rb,coord,&(is_per[0]));
        geomsOP[lev] = *geoms[lev];

        state[lev] = new MultiFab(ba,dm,nCompOut,nGrow);
        flame_normal[lev] = new MultiFab(ba,dm,AMREX_SPACEDIM,nGrow);
        cell_normal[lev] = new MultiFab(ba,dm,AMREX_SPACEDIM,nGrow);
        const Vector<Real>& delta = amrData.DxLevel()[lev];
        Real dxn[3];
        for (int i=0; i<AMREX_SPACEDIM; ++i) {
           dxn[i] = delta[i];
        }

        // Get input state data
        if (verbose) Print() << "Reading data for level " << lev << "\n";
        amrData.FillVar(*state[lev],lev,inVarNames,destFillComps);

        // Build progress variable from state at idCst, put into idProg
        MultiFab StateVar(*state[lev], amrex::make_alias, idCst, 1);
        MultiFab ProgressVar(*state[lev], amrex::make_alias, idProg, 1);
        for (MFIter mfi(*state[lev]); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.validbox();
            const auto& StateVarFab = StateVar.array(mfi); 
            const auto& ProgVarFab  = ProgressVar.array(mfi); 
            Real invdenom = 1.0 / (progMax - progMin);
            AMREX_HOST_DEVICE_PARALLEL_FOR_3D(bx, i, j, k,
            {
                ProgVarFab(i,j,k) = ( StateVarFab(i,j,k) - progMin ) * invdenom;
            });
        }
        ProgressVar.FillBoundary(geoms[lev]->periodicity());

        if (verbose) Print() << "Progress variable computed for level " << lev << "\n";

    }  // End lev loop

    if ( do_smooth ) {
       // Set a composite ABec solve to smooth the progress variable 
       // Try solving c^{n+1} - ∆t \nabla \cdot b_i \nabla c^{n+1} = c^{n)
       // In ABec terms  (\alpha A - \beta \nabla \cdot B_i \nabla) \phi = rhs
       // \alpha = 1, A = I
       // \beta = ∆t, b_i = ?? let's start with 1, switch to a D_c later if need be. Adapt ∆t accordingly ...
       // rhs = c^{n}
       
       LPInfo info;
       info.setAgglomeration(1);
       info.setConsolidation(1);
       info.setMetricTerm(false);

       MLABecLaplacian mlabec(geomsOP, grids, dmap, info);
       mlabec.setMaxOrder(4);

       const Real tol_rel = 1.e-12;
       const Real tol_abs = 1.e-12;

       // Problem with Periodic or Neumann BC 
       std::array<LinOpBCType, AMREX_SPACEDIM> lo_bc;
       std::array<LinOpBCType, AMREX_SPACEDIM> hi_bc;
       for (int idim = 0; idim< AMREX_SPACEDIM; idim++){
          if (is_per[idim] == 1) {
              lo_bc[idim] = hi_bc[idim] = LinOpBCType::Periodic;
          } else {
              lo_bc[idim] = hi_bc[idim] = LinOpBCType::Neumann;
          }
       }
       mlabec.setDomainBC(lo_bc,hi_bc);

       for (int lev = 0; lev < Nlev; ++lev)
       {
          // for problem with homogeneous Neumann BC, we need to pass nullptr
          mlabec.setLevelBC(lev, nullptr);
       }

       Real alpha = 1.0;
       Real beta = smooth_time;
       mlabec.setScalars(alpha, beta);

       for (int lev = 0; lev < Nlev; ++lev) {
           // Set A = I  at each level
           MultiFab acoef(state[lev]->boxArray(), state[lev]->DistributionMap(), 1, 0);
           acoef.setVal(1.0);
           mlabec.setACoeffs(lev, acoef);

           // Set b_i = 1.0  at each level
           Array<MultiFab,AMREX_SPACEDIM> face_bcoef; 
           for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
           {   
               const BoxArray& ba = amrex::convert(state[lev]->boxArray(),
                                                   IntVect::TheDimensionVector(idim));
               face_bcoef[idim].define(ba, state[lev]->DistributionMap(), 1, 0); 
               face_bcoef[idim].setVal(1.0);   
           }   
           mlabec.setBCoeffs(lev, amrex::GetArrOfConstPtrs(face_bcoef));
       }

       Vector<MultiFab> solution(Nlev);
       Vector<MultiFab> rhs(Nlev);
       for (int lev = 0; lev < Nlev; ++lev) {
           rhs[lev].define(grids[lev], dmap[lev], 1, 0);
           MultiFab::Copy(rhs[lev],*state[lev], idProg, 0, 1, 0);
           solution[lev].define(grids[lev], dmap[lev], 1, 1);
           solution[lev].setVal(0.0);
       }

       MLMG mlmg(mlabec);
       mlmg.setMaxIter(100);
       mlmg.setVerbose(1);
       mlmg.solve(GetVecOfPtrs(solution), GetVecOfConstPtrs(rhs), tol_rel, tol_abs);

       for (int lev = 0; lev < Nlev; ++lev) {
          MultiFab::Copy(*state[lev], solution[lev],0,idSmProg, 1, nGrow);
          state[lev]->FillBoundary(idSmProg,1,geoms[lev]->periodicity()); 
       }  
       if (verbose) Print() << "Progress variable smoothed successfully \n";
    }

    int idprogvar = do_smooth ? idSmProg : idProg;  

//  Compute curvature using LinearOperators
    // Set-up Poisson Linear Solver
    LPInfo info;
    info.setAgglomeration(1);
    info.setConsolidation(1);
    info.setMetricTerm(false);
    info.setMaxCoarseningLevel(0);
        
    for (int lev=0; lev<Nlev; ++lev)
    {
        const BoxArray ba = amrData.boxArray(lev);
        DistributionMapping dm(ba);

        if (verbose) Print() << "Starting mean curvature on level " << lev << "\n";

        // Get face gradients of progress variable 
        MLPoisson poisson({*geoms[lev]}, {ba}, {dm}, info);
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
        if ( lev > 0 ) {
           MultiFab* ProgVarCoarse = new MultiFab(state[lev-1]->boxArray(), state[lev-1]->DistributionMap(), 1, state[lev-1]->nGrow()); 
           MultiFab::Copy(*ProgVarCoarse, *state[lev-1], idprogvar, 0, 1, nGrow);
           poisson.setCoarseFineBC(ProgVarCoarse,2);
        }
        MultiFab ProgVar(ba, dm, 1, nGrow); 
        MultiFab::Copy(ProgVar, *state[lev], idprogvar, 0, 1, nGrow);
        poisson.setLevelBC(0,&ProgVar);

        MLMG mlmg(poisson);

        std::array<MultiFab,AMREX_SPACEDIM> face_gradient;
        AMREX_D_TERM(face_gradient[0].define(convert(ba,IntVect::TheDimensionVector(0)), dm, 1, 0); ,
                     face_gradient[1].define(convert(ba,IntVect::TheDimensionVector(1)), dm, 1, 0); ,
                     face_gradient[2].define(convert(ba,IntVect::TheDimensionVector(2)), dm, 1, 0); );
        mlmg.getFluxes({amrex::GetArrOfPtrs(face_gradient)},{&ProgVar});

        // Convert to cell avg gradient
        MultiFab cellavg_gradient(ba, dm, AMREX_SPACEDIM, 0);
        average_face_to_cellcenter(cellavg_gradient, 0, amrex::GetArrOfConstPtrs(face_gradient));
        cellavg_gradient.mult(-1.0,0,AMREX_SPACEDIM);

        // Compute ||\nabla c||
        MultiFab cellnorm_gradient(ba, dm, 1, 1);
        cellnorm_gradient.setVal(0.0);
        for (MFIter mfi(cellnorm_gradient); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.validbox();
            AMREX_D_TERM(const auto& Cx = cellavg_gradient.array(mfi,0);,
                         const auto& Cy = cellavg_gradient.array(mfi,1);,
                         const auto& Cz = cellavg_gradient.array(mfi,2); );
            const auto& normgrad  = cellnorm_gradient.array(mfi); 
            const auto& progvar   = ProgVar.array(mfi); 
            AMREX_HOST_DEVICE_PARALLEL_FOR_3D(bx, i, j, k,
            {
                normgrad(i,j,k) = std::max(1e-6, std::sqrt( AMREX_D_TERM (   std::pow(Cx(i,j,k),2.0),
                                                                           + std::pow(Cy(i,j,k),2.0),
                                                                           + std::pow(Cz(i,j,k),2.0)) ) ) ;
                normgrad(i,j,k) = -normgrad(i,j,k);
            });
        }
        cellnorm_gradient.FillBoundary(geoms[lev]->periodicity());

        // Get the cell centered flame normal N_i (pointing toward fresh gases)
        MultiFab::Copy(*cell_normal[lev],cellavg_gradient,0,0,AMREX_SPACEDIM,0);
        cell_normal[lev]->FillBoundary(0, AMREX_SPACEDIM, geoms[lev]->periodicity());

        if (verbose) Print() << "Done with flame normal on level " << lev << "\n";

//      At this point, I got the flame normal from clean gradients provided by the MLMG
//      Now try to compute the divergence of the flame normal using another linear solver:
//      1 -> get the clean face gradient of flame normal components from MLMG: d N_i/d x_i
//      2 -> manually get the divergence from the gradients: ∑_i d N_i/d x_i 

//      Copy into level aware flame_normal MF and fill same level ghost cells on normal 
        MultiFab::Copy(*flame_normal[lev], *cell_normal[lev], 0, 0, AMREX_SPACEDIM, 0);
        for (int idim = 0; idim< AMREX_SPACEDIM; idim++){
            MultiFab::Divide(*flame_normal[lev],cellnorm_gradient,0,idim,1,0);
        }  
        flame_normal[lev]->FillBoundary(0, AMREX_SPACEDIM, geoms[lev]->periodicity());

//      Define curvature        
        MultiFab Curv(ba, dm, 1, 0);     
        Curv.setVal(0.0); 

        for (int idim = 0; idim< AMREX_SPACEDIM; idim++){

           MLPoisson poisson2({*geoms[lev]}, {ba}, {dm}, info);
           poisson2.setMaxOrder(4);
           poisson2.setDomainBC(lo_bc, hi_bc);
           if ( lev > 0 ) {
               MultiFab* FlameNormalIdimCoarse = new MultiFab(flame_normal[lev-1]->boxArray(),
                                                              flame_normal[lev-1]->DistributionMap(),
                                                              1, 0); 
               MultiFab::Copy(*FlameNormalIdimCoarse, *flame_normal[lev-1], idim, 0, 1, 0);
               poisson2.setCoarseFineBC(FlameNormalIdimCoarse,2);
           }
           MultiFab* FlameNormalIdim = new MultiFab(ba, dm, 1, 1); 
           MultiFab::Copy(*FlameNormalIdim, *flame_normal[lev], idim, 0, 1, 1);
           poisson2.setLevelBC(0,FlameNormalIdim);
            
           MLMG mlmg2(poisson2);

           // Get the fluxes : d N_i / d x_j   , i = idim, j = 0, 1 (,2)
           std::array<MultiFab,AMREX_SPACEDIM> faceg;
           AMREX_D_TERM(faceg[0].define(convert(ba,IntVect::TheDimensionVector(0)), dm, 1, 0); ,
                        faceg[1].define(convert(ba,IntVect::TheDimensionVector(1)), dm, 1, 0); ,
                        faceg[2].define(convert(ba,IntVect::TheDimensionVector(2)), dm, 1, 0); );
           mlmg2.getFluxes({amrex::GetArrOfPtrs(faceg)},{FlameNormalIdim});

           // Get cell centered d N_i / d x_y
           MultiFab cell_avgg(ba, dm, AMREX_SPACEDIM, 0);
           average_face_to_cellcenter(cell_avgg, 0, amrex::GetArrOfConstPtrs(faceg));
           cell_avgg.mult(-1.0,0,AMREX_SPACEDIM);

           // Add d N_i / d x_i to curvature
           MultiFab::Add(Curv,cell_avgg,idim,0,1,0);
        }

#if AMREX_SPACEDIM == 3
        // Mean curvature : 0.5 * \div \cdot n
        // TODO: I only need to do that in 3D ... right ?
        Curv.mult(0.5,0,1);  
#endif

        // Clip curvature & flame normal for c < 0.01
        for (MFIter mfi(Curv); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.validbox();
            const auto& CurvFab  = Curv.array(mfi); 
            AMREX_D_TERM(const auto& FnormXFab  = flame_normal[lev]->array(mfi,0); ,
                         const auto& FnormYFab  = flame_normal[lev]->array(mfi,1); ,
                         const auto& FnormZFab  = flame_normal[lev]->array(mfi,2); );
            const auto& progvar   = ProgVar.array(mfi); 
            AMREX_HOST_DEVICE_PARALLEL_FOR_3D(bx, i, j, k,
            {
                if ( progvar(i,j,k) < 0.01 || progvar(i,j,k) > 0.99 ) {
                   CurvFab(i,j,k) = 0.0;
                   AMREX_D_TERM( FnormXFab(i,j,k) = 0.0; ,
                                 FnormYFab(i,j,k) = 0.0; ,
                                 FnormZFab(i,j,k) = 0.0; );
                }
            });
        }

        MultiFab::Copy(*state[lev], Curv, 0, idKm, 1, 0);
        MultiFab::Copy(*state[lev], *flame_normal[lev], 0, idN, AMREX_SPACEDIM, 0);

        if (verbose) Print() << "Mean curvature has been computed on level " << lev << "\n";

        // Now work on the gaussian curvature: only if 3D and required
#if AMREX_SPACEDIM == 3
        if ( do_gaussCurv ) { 

           // Start by getting the Hessian of the progress variable
           MultiFab Hessian(ba, dm, 9, 0);     

           // Compute grad of grad in each dim and store in Hessian
           for (int idim = 0; idim< AMREX_SPACEDIM; idim++){

              MLPoisson poisson2({*geoms[lev]}, {ba}, {dm}, info);
              poisson2.setMaxOrder(4);
              poisson2.setDomainBC(lo_bc, hi_bc);
              if ( lev > 0 ) {
                  MultiFab* gradIdimCoarse = new MultiFab(cell_normal[lev-1]->boxArray(),
                                                          cell_normal[lev-1]->DistributionMap(),
                                                          1, 0); 
                  MultiFab::Copy(*gradIdimCoarse, *cell_normal[lev-1], idim, 0, 1, 0);
                  poisson2.setCoarseFineBC(gradIdimCoarse,2);
              }
              MultiFab* gradIdim = new MultiFab(ba, dm, 1, 1); 
              MultiFab::Copy(*gradIdim, *cell_normal[lev], idim, 0, 1, 1);
              poisson2.setLevelBC(0,gradIdim);

              MLMG mlmg2(poisson2);

              std::array<MultiFab,AMREX_SPACEDIM> faceg;
              AMREX_D_TERM(faceg[0].define(convert(ba,IntVect::TheDimensionVector(0)), dm, 1, 0); ,
                           faceg[1].define(convert(ba,IntVect::TheDimensionVector(1)), dm, 1, 0); ,
                           faceg[2].define(convert(ba,IntVect::TheDimensionVector(2)), dm, 1, 0); );
              mlmg2.getFluxes({amrex::GetArrOfPtrs(faceg)},{gradIdim});

              // Get cell centered d C / d idim _x_y_z
              MultiFab cell_avgg(ba, dm, AMREX_SPACEDIM, 0);
              average_face_to_cellcenter(cell_avgg, 0, amrex::GetArrOfConstPtrs(faceg));
              cell_avgg.mult(-1.0,0,AMREX_SPACEDIM);

              // Copy in Hessian
              MultiFab::Copy(Hessian,cell_avgg,0,(idim)*3,AMREX_SPACEDIM,0);
           } 

           // Get the adjugate of the Hessian
           MultiFab AdjHessian(ba, dm, 9, 0);     

           for (MFIter mfi(Hessian); mfi.isValid(); ++mfi)
           {
               const Box& bx = mfi.validbox();
               const auto& HxiFab  = Hessian.array(mfi,0);           // Cxx, Cxy, Cxz
               const auto& HyiFab  = Hessian.array(mfi,3);           // Cyx, Cyy, Cyz
               const auto& HziFab  = Hessian.array(mfi,6);           // Czx, Czy, Czz
               const auto& AdjHxiFab  = AdjHessian.array(mfi,0); 
               const auto& AdjHyiFab  = AdjHessian.array(mfi,3); 
               const auto& AdjHziFab  = AdjHessian.array(mfi,6); 
               AMREX_HOST_DEVICE_PARALLEL_FOR_3D(bx, i, j, k,
               {
                   AdjHxiFab(i,j,k,0) = HyiFab(i,j,k,1) * HziFab(i,j,k,2) - HziFab(i,j,k,1) * HyiFab(i,j,k,2); 
                   AdjHyiFab(i,j,k,0) = HyiFab(i,j,k,2) * HziFab(i,j,k,0) - HziFab(i,j,k,2) * HyiFab(i,j,k,0);
                   AdjHziFab(i,j,k,0) = HyiFab(i,j,k,0) * HziFab(i,j,k,1) - HziFab(i,j,k,0) * HyiFab(i,j,k,1);
                   AdjHxiFab(i,j,k,1) = HxiFab(i,j,k,2) * HziFab(i,j,k,1) - HziFab(i,j,k,2) * HxiFab(i,j,k,1);
                   AdjHyiFab(i,j,k,1) = HxiFab(i,j,k,0) * HziFab(i,j,k,2) - HziFab(i,j,k,0) * HxiFab(i,j,k,2);
                   AdjHziFab(i,j,k,1) = HxiFab(i,j,k,1) * HziFab(i,j,k,0) - HziFab(i,j,k,1) * HxiFab(i,j,k,0);
                   AdjHxiFab(i,j,k,2) = HxiFab(i,j,k,1) * HyiFab(i,j,k,2) - HyiFab(i,j,k,1) * HxiFab(i,j,k,2);
                   AdjHyiFab(i,j,k,2) = HxiFab(i,j,k,2) * HyiFab(i,j,k,0) - HyiFab(i,j,k,2) * HxiFab(i,j,k,0);
                   AdjHziFab(i,j,k,2) = HxiFab(i,j,k,0) * HyiFab(i,j,k,1) - HyiFab(i,j,k,0) * HxiFab(i,j,k,1); 
               });
           }

           // Now get the gausian curvature
           MultiFab gCurv(ba, dm, 1, 0);     
           for (MFIter mfi(gCurv); mfi.isValid(); ++mfi)
           {
               const Box& bx = mfi.validbox();
               const auto& progvar    = ProgVar.array(mfi);
               const auto& gCurvFab   = gCurv.array(mfi); 
               const auto& CgradNorm  = cellnorm_gradient.array(mfi); 
               const auto& CxFab      = cellavg_gradient.array(mfi,0); 
               const auto& CyFab      = cellavg_gradient.array(mfi,1); 
               const auto& CzFab      = cellavg_gradient.array(mfi,2); 
               const auto& AdjHxiFab  = AdjHessian.array(mfi,0); 
               const auto& AdjHyiFab  = AdjHessian.array(mfi,3); 
               const auto& AdjHziFab  = AdjHessian.array(mfi,6); 
               AMREX_HOST_DEVICE_PARALLEL_FOR_3D(bx, i, j, k,
               {
                  gCurvFab(i,j,k) = ( CxFab(i,j,k) * ( AdjHxiFab(i,j,k,0) * CxFab(i,j,k) +
                                                       AdjHxiFab(i,j,k,1) * CyFab(i,j,k) + 
                                                       AdjHxiFab(i,j,k,2) * CzFab(i,j,k) ) +
                                      CyFab(i,j,k) * ( AdjHyiFab(i,j,k,0) * CxFab(i,j,k) +
                                                       AdjHyiFab(i,j,k,1) * CyFab(i,j,k) +
                                                       AdjHyiFab(i,j,k,2) * CzFab(i,j,k) ) +
                                      CzFab(i,j,k) * ( AdjHziFab(i,j,k,0) * CxFab(i,j,k) +
                                                       AdjHziFab(i,j,k,1) * CyFab(i,j,k) +
                                                       AdjHziFab(i,j,k,2) * CzFab(i,j,k) ) 
                                    ) / std::pow(CgradNorm(i,j,k),4.0);
                   if ( progvar(i,j,k) < 0.01 || progvar(i,j,k) > 0.99 ) {
                      gCurvFab(i,j,k) = 0.0;
                   }
               });
           }
           MultiFab::Copy(*state[lev], gCurv, 0, idKg, 1, 0);
        }
        if (verbose) Print() << "Gaussian curvature has been computed on level " << lev << "\n";
#endif

        if (do_strain) { 
           // Strain rate -nn:\nabla u + \nabla \cdot u

           // Start by building the strain tensor
           MultiFab StrainT(ba, dm, AMREX_SPACEDIM * AMREX_SPACEDIM, 0);     

           // Compute strain tensor with a MLMG in each direction
           for (int idim = 0; idim< AMREX_SPACEDIM; idim++){

              MLPoisson poisson2({*geoms[lev]}, {ba}, {dm}, info);
              poisson2.setMaxOrder(4);
              poisson2.setDomainBC(lo_bc, hi_bc);
              if ( lev > 0 ) {
                  MultiFab* velIdimCoarse = new MultiFab(state[lev-1]->boxArray(),
                                                         state[lev-1]->DistributionMap(),
                                                         1, 0); 
                  MultiFab::Copy(*velIdimCoarse, *state[lev-1], idVst+idim, 0, 1, 0);
                  poisson2.setCoarseFineBC(velIdimCoarse,2);
              }
              MultiFab* velIdim = new MultiFab(ba, dm, 1, 1); 
              MultiFab::Copy(*velIdim, *state[lev], idVst+idim, 0, 1, 1);
              poisson2.setLevelBC(0,velIdim);

              MLMG mlmg2(poisson2);

              std::array<MultiFab,AMREX_SPACEDIM> faceg;
              AMREX_D_TERM(faceg[0].define(convert(ba,IntVect::TheDimensionVector(0)), dm, 1, 0); ,
                           faceg[1].define(convert(ba,IntVect::TheDimensionVector(1)), dm, 1, 0); ,
                           faceg[2].define(convert(ba,IntVect::TheDimensionVector(2)), dm, 1, 0); );
              mlmg2.getFluxes({amrex::GetArrOfPtrs(faceg)},{velIdim});

              // Get cell centered d u_idim / d _x_y(_z)
              MultiFab cell_avgg(ba, dm, AMREX_SPACEDIM, 0);
              average_face_to_cellcenter(cell_avgg, 0, amrex::GetArrOfConstPtrs(faceg));
              cell_avgg.mult(-1.0,0,AMREX_SPACEDIM);

              // Copy in strain tensor
              MultiFab::Copy(StrainT,cell_avgg,0,(idim)*AMREX_SPACEDIM,AMREX_SPACEDIM,0);
           } 

           // Gather the components of strain rate
           MultiFab strainrate(ba, dm, 1, 0);

           for (MFIter mfi(strainrate); mfi.isValid(); ++mfi)
           {
               const Box& bx = mfi.validbox();
               const auto& progvar    = ProgVar.array(mfi);
               const auto& srFab      = strainrate.array(mfi); 
               AMREX_D_TERM(const auto& NxFab  = flame_normal[lev]->array(mfi,0); ,
                            const auto& NyFab  = flame_normal[lev]->array(mfi,1); ,
                            const auto& NzFab  = flame_normal[lev]->array(mfi,2); );
               AMREX_D_TERM(const auto& gradUx = StrainT.array(mfi,0); ,
                            const auto& gradUy = StrainT.array(mfi,AMREX_SPACEDIM); ,
                            const auto& gradUz = StrainT.array(mfi,AMREX_SPACEDIM*2); ); 
               AMREX_HOST_DEVICE_PARALLEL_FOR_3D(bx, i, j, k,
               {
                  srFab(i,j,k) = AMREX_D_TERM ( AMREX_D_TERM( - gradUx(i,j,k,0) * NxFab(i,j,k) * NxFab(i,j,k) ,
                                                              - gradUx(i,j,k,1) * NxFab(i,j,k) * NyFab(i,j,k) ,
                                                              - gradUx(i,j,k,2) * NxFab(i,j,k) * NzFab(i,j,k) ),
                                                AMREX_D_TERM( - gradUy(i,j,k,0) * NyFab(i,j,k) * NxFab(i,j,k) ,
                                                              - gradUy(i,j,k,1) * NyFab(i,j,k) * NyFab(i,j,k) ,
                                                              - gradUy(i,j,k,2) * NyFab(i,j,k) * NzFab(i,j,k) ),
                                                AMREX_D_TERM( - gradUz(i,j,k,0) * NzFab(i,j,k) * NxFab(i,j,k) ,
                                                              - gradUz(i,j,k,1) * NzFab(i,j,k) * NyFab(i,j,k) ,
                                                              - gradUz(i,j,k,2) * NzFab(i,j,k) * NzFab(i,j,k) ) );
               });
           }
        }

    }

    // ---------------------------------------------------------------------
    // Set-up the output
    // ---------------------------------------------------------------------
    Vector<std::string> nnames(nCompOut);

    // Plot file variables
    for (int i=0; i<nCompIn; ++i) {
        nnames[i] = inVarNames[i];
    }

    // Computed variables
    nnames[idProg]   = "Progress";
    nnames[idSmProg] = "SmoothedProgress";
    nnames[idKm]     = "MeanCurvature_" + progressName;
    nnames[idN]      = "FlameNormalX_" + progressName;
    nnames[idN+1]    = "FlameNormalY_" + progressName;
#if AMREX_SPACEDIM == 3
    nnames[idN+2]    = "FlameNormalZ_" + progressName;
    nnames[idKg]     = "GaussianCurvature_" + progressName;
#endif
    if (do_strain) nnames[idSR] = "StrainRate_" + progressName;

    if (getStrainTensor)
    {
        for (int i=0; i<9; ++i)
        {
            int n=i%3;
            int m=i/3;
            char buf[40];
            sprintf(buf,"ROST_%1d%1d",n+1,m+1);
            nnames[idROST+i] = string(buf);
        }
    }
    
    if (do_velnormal)
    {
        nnames[idVelNormal+0] = "VelFlameNormal_X";
        nnames[idVelNormal+1] = "VelFlameNormal_Y";
#if AMREX_SPACEDIM>2
        nnames[idVelNormal+2] = "VelFlameNormal_Z";
#endif
    }

    if (getProgGrad)
    {
        nnames[idProgGrad+0] = progressName + "_g1";
        nnames[idProgGrad+1] = progressName + "_g2";
        nnames[idProgGrad+2] = progressName + "_g3";
    }

    // Write to plotfile
    bool verb=true;
    if (appendPlotFile)
    {

        int nStateOut = nCompOut - nCompIn;
        Vector<MultiFab*> ostate(Nlev);
        for (int lev=0; lev<Nlev; ++lev)
        {
            const BoxArray ba = state[lev]->boxArray();
            ostate[lev] = new MultiFab(ba,DistributionMapping(ba),nStateOut,0);
            MultiFab::Copy(*ostate[lev],*state[lev],nCompIn,0,nStateOut,0);
        }
        Vector<std::string> namesOut(nStateOut);
        for (int i=0; i<nStateOut; ++i)
            namesOut[i] = nnames[nCompIn+i];

        std::string newMFBaseName = "NEWDAT"; pp.query("newMFBaseName",newMFBaseName);
        std::string newHeaderName = "NewHeader"; pp.query("newHeaderName",newHeaderName);
        AppendToPlotFile(amrData,ostate,plotFileName,namesOut,newMFBaseName,newHeaderName,verb);
        Print() << "...finished.  Note: to see new data, you must rename NewHeader in the" << "\n";
        Print() << "              pltfile to Header (probably want to save the original first)" << "\n";
    }
    else
    {
        // Remove GC from outstate
        Vector<MultiFab*> ostate(Nlev);
        for (int lev=0; lev<Nlev; ++lev)
        {
            const BoxArray ba = state[lev]->boxArray();
            ostate[lev] = new MultiFab(ba,DistributionMapping(ba),nCompOut,0);
            MultiFab::Copy(*ostate[lev],*state[lev],0,0,nCompOut,0);
        }
        Print() << "Writing new data to " << outfile << "\n";
        WritePlotFile(ostate,amrData,outfile,verb,nnames);
    }

    }
    amrex::Finalize();
    return 0;
}
