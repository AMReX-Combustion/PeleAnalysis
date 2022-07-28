#include <string>
#include <iostream>
#include <set>

#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include <AMReX_DataServices.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_WritePlotFile.H>
#include <AMReX_Sundials.H>

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
  std::cerr << argv[0] << " infile infile=f1 [options] \n\tOptions:\n";
  exit(1);
}

std::string
getFileRoot(const std::string& infile)
{
  // vector<std::string> tokens = Tokenize(infile,std::string("/"));
  std::vector<std::string> tokens = Tokenize(infile,std::string("/"));
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

    std::string chem_integrator = "ReactorCvode";
    std::string solve_type_str = "none";

    int ode_ncells = 1;
    int ode_iE = 1;
    // pp.get("chem_integrator", chem_integrator);
    // pele::physics::reactions::ReactorBase::create(chem_integrator);
    pele::physics::transport::TransportParams<
      pele::physics::PhysicsType::transport_type>
      trans_parms;
    trans_parms.allocate();

    amrex::sundials::Initialize(amrex::OpenMP::get_max_threads());

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

    std::string output_folder; pp.get("outFolder",output_folder);
    std::string plotFileName; pp.get("infile",plotFileName);
    DataServices::SetBatchMode();
    Amrvis::FileType fileType(Amrvis::NEWPLT);
    Print() << "Processing " << plotFileName << std::endl;

    DataServices dataServices(plotFileName, fileType);
    if( ! dataServices.AmrDataOk()) {
      DataServices::Dispatch(DataServices::ExitRequest, NULL);
      // ^^^ this calls ParallelDescriptor::EndParallel() and exit()
    }
    AmrData& amrData = dataServices.AmrDataRef();

    int finestLevel = amrData.FinestLevel();
    pp.query("finestLevel",finestLevel);
    int Nlev = finestLevel + 1;
    Print() << "Working with " << Nlev << " levels" << std::endl;

    Vector<std::unique_ptr<MultiFab>> indata(Nlev);
    MultiFab outdata;
    // Vector<std::unique_ptr<MultiFab>> outdata(Nlev);

    int ngrow = 1;
    int ncomp = amrData.NComp();
    Print() << "Number of scalars in the original file: " << ncomp << std::endl;

    // const auto& domain = amrData.ProbDomain()[finestLevel];
    // BoxArray ba(domain);
    // int max_grid_size = 32;
    // pp.query("max_grid_size",max_grid_size);
    // ba.maxSize(max_grid_size);
    // const DistributionMapping dm(ba);

    int plt_src; pp.get("plt_src",plt_src);
    
    int idRlocal  = -1; //density index
    int idTlocal  = -1; // T  index
    int idYlocal  = -1; // Ys index
    int idvFlocal = -1; // volume fraction index
    int idMixFrac = -1; // mixture fraction index
    int idScalDis = ncomp; // mixture fraction index
    const int nCompOut = ncomp+1;

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
    for (int i=0; i<amrData.PlotVarNames().size(); ++i)
    {
      if (amrData.PlotVarNames()[i] == spName)    idYlocal = i;
      if (amrData.PlotVarNames()[i] == TName)     idTlocal = i;
      if (amrData.PlotVarNames()[i] == RHOName)   idRlocal = i;
      if (amrData.PlotVarNames()[i] == vFracName) idvFlocal = i;
      if (amrData.PlotVarNames()[i] == ZName)     idMixFrac = i;
    }

    Print() << "Temperature index: " << idTlocal  << std::endl;
    Print() << "Density index    : " << idRlocal  << std::endl;
    Print() << "Vfrac index      : " << idvFlocal << std::endl;
    Print() << "First spec index : " << idYlocal  << std::endl;
    Print() << "Mix Frac index   : " << idMixFrac  << std::endl;
    Print() << "Scaldis  index   : " << idScalDis  << std::endl;

    Vector<std::string> outNames(nCompOut);
    for (int i=0; i<amrData.PlotVarNames().size(); ++i)
    {
      outNames[i] = amrData.PlotVarNames()[i];
    }
    outNames[ncomp] = "scaldis";

    RealBox rb(&(amrData.ProbLo()[0]),
               &(amrData.ProbHi()[0]));
    Vector<int> is_per(AMREX_SPACEDIM,0);

    int coord = 0;
    int max_grid_size = 32;
    pp.query("max_grid_size",max_grid_size);

    for (int lev=0; lev<Nlev; ++lev)
    {        
      Print() << "Working on level " << lev << std::endl; 
      //define domain stuff
      const auto& domain = amrData.ProbDomain()[lev];
      BoxArray ba(domain);
      ba.maxSize(max_grid_size);
      const DistributionMapping dm(ba);

      Print() << "Number of boxes: " << ba.size() << std::endl; 

      //declare MultiFabs for transport properties
      int num_grow = 0;
      amrex::MultiFab D(ba, dm, NUM_SPECIES, num_grow);
      amrex::MultiFab mu (ba, dm, 1, num_grow);
      amrex::MultiFab xi (ba, dm, 1, num_grow);
      amrex::MultiFab lam(ba, dm, 1, num_grow);

      // Data MFs
      amrex::MultiFab mass_frac(ba, dm, NUM_SPECIES, num_grow);
      amrex::MultiFab temperature(ba, dm, 1, num_grow);
      amrex::MultiFab density(ba, dm, 1, num_grow);

      //vector containing the number of components to be read
      Vector<int> destFillComps(ncomp);
      for (int i=0; i<ncomp; ++i) destFillComps[i] = i;

      //reset MultiFab pointer
      indata[lev].reset(new MultiFab(ba,dm,ncomp,ngrow));

      Print() << "Reading data..." << std::endl; 
      amrData.FillVar(*indata[lev],lev,amrData.PlotVarNames(),destFillComps);
      Print() << "Done!" << std::endl;
      if (lev == Nlev - 1) {
        int ngrow_out = 0;
        // outdata[lev].reset(new MultiFab(ba,dm,ncomp,ngrow));
        outdata.define(ba,dm,nCompOut,ngrow_out);

        Geometry geom(amrData.ProbDomain()[lev],&rb,coord,&(is_per[0]));
        
        // Initialize Fabs with data from plt file
        for (MFIter mfi(*indata[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
          const Box& box = mfi.tilebox();

          auto const& Y_a = mass_frac.array(mfi);
          auto const& T_a = temperature.array(mfi);
          auto const& rho_a = density.array(mfi);
          amrex::Array4<const amrex::Real> sfab_init = (*indata[lev]).array(mfi);

          // amrex::ParallelFor(
          //   box, [Y_a, T_a, rho_a,
          //        geom] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
          AMREX_PARALLEL_FOR_3D ( box, i, j, k,
          {
                  
                  for (int n = 0; n < NUM_SPECIES; n++){
                    Y_a(i,j,k,n) = sfab_init(i,j,k,n+idYlocal);
                  }

                  T_a  (i,j,k) = sfab_init(i,j,k,idTlocal);
                  rho_a(i,j,k) = sfab_init(i,j,k,idRlocal);
                ;
            });
        }

        // Get the transport data pointer
        auto const* ltransparm = trans_parms.device_trans_parm();

        //Now that we have all the Fabs with data lets get the transp. coeff. and compute scaldis
        for (MFIter mfi(*indata[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
          const Box& box = mfi.tilebox();

          amrex::Array4<amrex::Real> const& Y_a = mass_frac.array(mfi);
          amrex::Array4<amrex::Real> const& T_a = temperature.array(mfi);
          amrex::Array4<amrex::Real> const& rho_a = density.array(mfi);
          amrex::Array4<amrex::Real> const& D_a = D.array(mfi);
          amrex::Array4<amrex::Real> const& mu_a = mu.array(mfi);
          amrex::Array4<amrex::Real> const& xi_a = xi.array(mfi);
          amrex::Array4<amrex::Real> const& lam_a = lam.array(mfi);

          amrex::launch(box, [=] AMREX_GPU_DEVICE(amrex::Box const& tbx) {
          auto trans = pele::physics::PhysicsType::transport();
          
          trans.get_transport_coeffs(
            tbx, Y_a, T_a, rho_a, D_a, mu_a, xi_a, lam_a, ltransparm);
          
          });

          amrex::Array4<const amrex::Real> sfab = (*indata[lev]).array(mfi);
          amrex::Array4<Real> const sfab_out = outdata.array(mfi);
          // auto& sfab_out = (*outdata[0])[mfi];
          const Real* dx = geom.CellSize();
          
          amrex::Real dxInv = 0.5 / dx[0];

          // amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
          AMREX_PARALLEL_FOR_3D ( box, i, j, k,
          {

            sfab_out(i,j,k,ncomp) = 2.*D_a(i,j,k,N2_ID)*pow(dxInv*(sfab(i+1,j,k,idMixFrac) - sfab(i-1,j,k,idMixFrac)),2);

          });

        }
      }
    }

    std::string outfile(output_folder+getFileRoot(plotFileName) + "_scaldis");
    Print() << "Writing new data to " << outfile << std::endl;
    bool verb = false;

    // This was working in 3D but is not in 2D
    // WritePlotFile(GetVecOfPtrs(outdata),amrData,outfile,verb,outNames);
    
    // Output the finest level only    
    int ref_ratio = 2;
//     Vector<Geometry> geoms(1);
//     Vector<int> levelSteps(1);
    Geometry geoms;
    int levelSteps;
    {
      geoms = Geometry(amrData.ProbDomain()[Nlev - 1],&rb,coord,&(is_per[0]));
      levelSteps = 0;
    }

    WriteSingleLevelPlotfile(outfile, outdata, outNames, geoms, amrData.Time(), levelSteps);

  }
  Finalize();
  return 0;
}
