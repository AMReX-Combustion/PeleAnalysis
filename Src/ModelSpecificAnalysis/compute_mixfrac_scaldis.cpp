#include <string>
#include <iostream>

#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include <AMReX_DataServices.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_BLFort.H>
#include <AMReX_MLMG.H>
#include <AMReX_MLPoisson.H>
#include <AMReX_MLABecLaplacian.H>

#include "mechanism.H"
#include <PelePhysics.H>
#include <ReactorBase.H>

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
    std::string gradVar       = "mixture_fraction";
    std::string infile        = "";  
    int finestLevel           = 1000;
    int nAuxVar               = 0;

    // ---------------------------------------------------------------------
    // ParmParse
    // ---------------------------------------------------------------------
    ParmParse pp;

    if (pp.contains("help")) {
      print_usage(argc,argv);
    }

    pp.get("infile",infile);
    pp.query("gradVar",gradVar);
    pp.query("finestLevel",finestLevel);
    int max_grid_size = 32;pp.query("max_grid_size",max_grid_size);
    int plt_src = 1; pp.query("plt_src",plt_src);


    // pp.get("chem_integrator", chem_integrator);
    // pele::physics::reactions::ReactorBase::create(chem_integrator);

    // Initialize transport module
    pele::physics::transport::TransportParams<
      pele::physics::PhysicsType::transport_type>
      trans_parms;
    trans_parms.allocate();

    // Initialize reactor module
    std::string chem_integrator = "ReactorCvode";
    amrex::sundials::Initialize(amrex::OpenMP::get_max_threads());

    int ode_ncells = 1;
    int ode_iE = 1;
    std::unique_ptr<pele::physics::reactions::ReactorBase> reactor =
      pele::physics::reactions::ReactorBase::create(chem_integrator);
    reactor->init(ode_iE, ode_ncells);

    Vector<std::string> spec_names;
    auto eos = pele::physics::PhysicsType::eos();
    // pele::physics::eos::speciesNames(spec_names);
    pele::physics::eos::speciesNames<pele::physics::PhysicsType::eos_type>(spec_names);
    for (int i=0; i<spec_names.size(); ++i) {
      Print() << spec_names[i] << " ";
    }
    Print() << std::endl;

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

    int ncomp = amrData.NComp();
    Print() << "Number of scalars in the original file: " << ncomp << std::endl;

    //=== Get index for relevant scalars present in the plt file ===
    int idRlocal  = -1; //density index
    int idTlocal  = -1; // T  index
    int idYlocal  = -1; // Ys index
    int idvFlocal = -1; // volume fraction index
    int idMixFrac = -1; // mixture fraction index
    int idScalDis = ncomp; // mixture fraction index
    const int nCompOutput = ncomp+1;

    const std::string spName  = "Y(" + spec_names[0] + ")";
    const std::string RHOName = "density";
    std::string TName;
    std::string vFracName;
    amrex::Real convert_rho;
    if(plt_src == 0){//PeleC
      Print() << "Assuming plt file was generated in PeleC" << std::endl;
      TName       = "Temp";
      vFracName   = "vfrac";
      convert_rho = 1.0;
    }
    else{ //PeleLM
      Print() << "Assuming plt file was generated in PeleLM" << std::endl;
      TName       = "temp";
      vFracName   = "volFrac";
      convert_rho = 1.0e-3;
    }
    std::string ZName = "mixture_fraction";

    // Gradient variable
    int idC = -1;
    for (int i=0; i<plotVarNames.size(); ++i)
    {
      if (plotVarNames[i] == gradVar) idC = i;
      if (plotVarNames[i] == spName)    idYlocal = i;
      if (plotVarNames[i] == TName)     idTlocal = i;
      if (plotVarNames[i] == RHOName)   idRlocal = i;
      if (plotVarNames[i] == ZName)     idMixFrac = i;
    }
    
    Print() << "Temperature index: " << idTlocal  << std::endl;
    Print() << "Density index    : " << idRlocal  << std::endl;
    Print() << "First spec index : " << idYlocal  << std::endl;

    if (idC<0) {
      Print() << "Cannot find " << gradVar << " data in pltfile \n";
    }
    // ===========================================

    // Auxiliary variables
    nAuxVar = pp.countval("Aux_Variables");
    Vector<std::string> AuxVar(nAuxVar);
    for(int ivar = 0; ivar < nAuxVar; ++ivar) { 
         pp.get("Aux_Variables", AuxVar[ivar],ivar);
    }

    // ---------------------------------------------------------------------
    // Variables index management
    // ---------------------------------------------------------------------
    const int idCst = 0;
    int nCompIn = idCst + 1;
    Vector<std::string> inVarNames(nCompIn);
    inVarNames[idCst] = plotVarNames[idC];

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
    Vector<int> destFillComps_indata(ncomp);
    for (int i=0; i<ncomp; ++i) {
      destFillComps_indata[i] = i;
    }

    const int idGr = nCompIn;
    const int nCompOut = idGr + AMREX_SPACEDIM +1 ; // 1 component stores the ||gradT||

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
    Vector<MultiFab> state(Nlev);  //At the moment this FAB stores the field which we will compute the gradient on
    Vector<MultiFab> indata(Nlev); //This FAB stores all the variables in the plt file
    Vector<Geometry> geoms(Nlev);
    Vector<BoxArray> grids(Nlev);
    Vector<DistributionMapping> dmap(Nlev);
    const int nGrow = 1;

    // Read data on all the levels
    for (int lev=0; lev<Nlev; ++lev) {

      const BoxArray ba = amrData.boxArray(lev);
      grids[lev] = ba;
      dmap[lev] = DistributionMapping(ba);
      geoms[lev] = Geometry(amrData.ProbDomain()[lev],&rb,coord,&(is_per[0]));
      state[lev] .define(grids[lev], dmap[lev], nCompOut, nGrow);
      indata[lev].define(grids[lev], dmap[lev], ncomp   , 0);

      Print() << "Reading data for level: " << lev << std::endl;
      amrData.FillVar(state[lev], lev, inVarNames, destFillComps);
      amrData.FillVar(indata[lev], lev, amrData.PlotVarNames(), destFillComps_indata);

      state[lev].FillBoundary(idCst,1,geoms[lev].periodicity());
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
    int nGrowGrad = 0;                   // No need for ghost face on gradient
    Vector<Array<MultiFab,AMREX_SPACEDIM> > grad(Nlev);
    Vector<std::unique_ptr<MultiFab>> phi;
    Vector<MultiFab> laps;
    for (int lev = 0; lev < Nlev; ++lev) {
      for (int idim = 0; idim <AMREX_SPACEDIM; idim++) {
         const auto& ba = grids[lev];
         grad[lev][idim].define(amrex::convert(ba,IntVect::TheDimensionVector(idim)),
                                dmap[lev], 1, nGrowGrad);
      }    
      phi.push_back(std::make_unique<MultiFab> (state[lev],amrex::make_alias,idCst,1));
      poisson.setLevelBC(lev, phi[lev].get());
      laps.emplace_back(grids[lev], dmap[lev], 1, 1);
    }

    MLMG mlmg(poisson);
    mlmg.apply(GetVecOfPtrs(laps), GetVecOfPtrs(phi));
    mlmg.getFluxes(GetVecOfArrOfPtrs(grad), GetVecOfPtrs(phi), MLMG::Location::FaceCenter);

    //declare MultiFabs for transport properties
    int num_grow = 0;
    Vector<MultiFab> D  (Nlev);
    Vector<MultiFab> mu (Nlev);
    Vector<MultiFab> xi (Nlev);
    Vector<MultiFab> lam(Nlev);

    // Data MFs
    Vector<MultiFab> mass_frac  (Nlev);
    Vector<MultiFab> temperature(Nlev);
    Vector<MultiFab> density    (Nlev);

    for (int lev = 0; lev < Nlev; ++lev) {
        // const auto& ba = grids[lev];
        // const auto& domain = amrData.ProbDomain()[lev];
        // BoxArray ba(domain);
        // ba.maxSize(max_grid_size);
        // const DistributionMapping dm(ba);

        //declare MultiFabs for transport properties
        int num_grow = 0;
        D  [lev].define(grids[lev], dmap[lev], NUM_SPECIES, num_grow);
        mu [lev].define(grids[lev], dmap[lev], 1, num_grow);
        xi [lev].define(grids[lev], dmap[lev], 1, num_grow);
        lam[lev].define(grids[lev], dmap[lev], 1, num_grow);

        // Data MFs
        mass_frac  [lev].define(grids[lev], dmap[lev], NUM_SPECIES, num_grow);
        temperature[lev].define(grids[lev], dmap[lev], 1, num_grow);
        density    [lev].define(grids[lev], dmap[lev], 1, num_grow);

        // Initialize Fabs with data from plt file
        for (MFIter mfi(indata[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
          const Box& box = mfi.tilebox();

          auto const& Y_a       = mass_frac[lev].array(mfi);
          auto const& T_a       = temperature[lev].array(mfi);
          auto const& rho_a     = density[lev].array(mfi);
          auto const& sfab_init = indata[lev].array(mfi);
          // amrex::Array4<const amrex::Real> sfab_init = indata[lev].array(mfi);

          // amrex::ParallelFor(
          //   box, [Y_a, T_a, rho_a,
          //        geom] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
          AMREX_PARALLEL_FOR_3D ( box, i, j, k,
          {
            
            // Print() << "Massfrac: " << sfab_init(i,j,k,0) << std::endl;               
                for (int n = 0; n < NUM_SPECIES; n++){
                    Y_a(i,j,k,n) = sfab_init(i,j,k,n+idYlocal);
                }

                T_a  (i,j,k) = sfab_init(i,j,k,idTlocal);
                rho_a(i,j,k) = sfab_init(i,j,k,idRlocal);
                
            });
        }

        // Get the transport data pointer
        auto const* ltransparm = trans_parms.device_trans_parm();

        // Convert to cell avg gradient
        MultiFab gradAlias(state[lev], amrex::make_alias, idGr, AMREX_SPACEDIM);
        average_face_to_cellcenter(gradAlias, 0, GetArrOfConstPtrs(grad[lev]));
        gradAlias.mult(-1.0);
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(state[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {

           const Box& box = mfi.tilebox();
           
           //Now that we have all the Fabs with data lets get the transp. coeff. and compute scaldis
           amrex::Array4<amrex::Real> const& Y_a   = mass_frac[lev].array(mfi);
           amrex::Array4<amrex::Real> const& T_a   = temperature[lev].array(mfi);
           amrex::Array4<amrex::Real> const& rho_a = density[lev].array(mfi);
           amrex::Array4<amrex::Real> const& D_a   = D[lev].array(mfi);
           amrex::Array4<amrex::Real> const& mu_a  = mu[lev].array(mfi);
           amrex::Array4<amrex::Real> const& xi_a  = xi[lev].array(mfi);
           amrex::Array4<amrex::Real> const& lam_a = lam[lev].array(mfi);

           amrex::launch(box, [=] AMREX_GPU_DEVICE(amrex::Box const& tbx) {
           auto trans = pele::physics::PhysicsType::transport();
          
           trans.get_transport_coeffs(
             tbx, Y_a, T_a, rho_a, D_a, mu_a, xi_a, lam_a, ltransparm);
           
             });

           const Box& bx = mfi.tilebox();
           auto const& grad_a   = gradAlias.const_array(mfi);
           auto const& gradMag  = state[lev].array(mfi,idGr+AMREX_SPACEDIM);
           amrex::ParallelFor(bx, [=]
           AMREX_GPU_DEVICE (int i, int j, int k) noexcept
           {    
              gradMag(i,j,k) = std::sqrt(AMREX_D_TERM(  grad_a(i,j,k,0) * grad_a(i,j,k,0),
                                                      + grad_a(i,j,k,1) * grad_a(i,j,k,1),
                                                      + grad_a(i,j,k,2) * grad_a(i,j,k,2)));

              gradMag(i,j,k) = 2.*D_a(i,j,k,N2_ID)*pow(gradMag(i,j,k),2);
           });  
        } 
    }

    // ---------------------------------------------------------------------
    // Write the results
    // ---------------------------------------------------------------------
    
    // Combine the data from state FAB containing gradient info and data from indata FAB with all other scalars
    Vector<MultiFab> data_output(Nlev);
    int ncomp_all = ncomp + nCompIn;
    int ngrow = 0;
    for (int lev = 0; lev < Nlev; ++lev) {
        // Box box_output = amrData.boxArray(lev);
        const BoxArray box_output = amrData.boxArray(lev);

        BoxArray ba_sub(box_output);
        ba_sub.maxSize(max_grid_size);
        DistributionMapping dmap_sub(ba_sub);

        data_output[lev].define(grids[lev], dmap[lev],ncomp_all,ngrow);
        // data_output[lev].define(grids[lev], dmap[lev], NUM_SPECIES, num_grow);
        data_output[lev].ParallelCopy(indata[lev],0,0,ncomp);
        data_output[lev].ParallelCopy(state[lev] ,idGr+AMREX_SPACEDIM,ncomp,1);
    }

    Vector<std::string> nnames(ncomp_all);
    for (int i=0; i<ncomp; ++i) {
      nnames[i] = amrData.PlotVarNames()[i];
    }
    nnames[ncomp] = "scaldis";
//     nnames[idGr+0] = gradVar + "_gx";
//     nnames[idGr+1] = gradVar + "_gy";
// #if AMREX_SPACEDIM==3
//     nnames[idGr+2] = gradVar + "_gz";
// #endif
//     nnames[idGr+AMREX_SPACEDIM] = "||grad"+ gradVar+ "||";
    std::string output_folder; pp.get("outFolder",output_folder);
    std::string outfile(output_folder+getFileRoot(infile) + "_scaldis"); pp.query("outfile",outfile);

    Print() << "Writing new data to " << outfile << std::endl;
    Vector<int> isteps(Nlev, 0);
    Vector<IntVect> refRatios(Nlev-1,{AMREX_D_DECL(2, 2, 2)});
    amrex::WriteMultiLevelPlotfile(outfile, Nlev, GetVecOfConstPtrs(data_output), nnames,
                                   geoms, amrData.Time(), isteps, refRatios);
  }


//   }
  amrex::Finalize();
  return 0;
}