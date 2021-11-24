#include <string>
#include <iostream>
#include <set>

#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include <AMReX_DataServices.H>
#include <WritePlotFile.H>
#include <AMReX_PlotFileUtil.H>

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
  vector<std::string> tokens = Tokenize(infile,std::string("/"));
  return tokens[tokens.size()-1];
}

extern "C" {
    void process(const int* lo,  const int* hi,
                 const int* dlo, const int* dhi,
                 Real* U, const int* Ulo, const int* Uhi,
                 const int* nc, const amrex::Real* plo,  const amrex::Real* dx);
};

void
init_composition(amrex::Real massfrac[NUM_SPECIES])
{
  amrex::Real vt, ek, a, yl, yr, sumY, T, rho, e;

  ParmParse pp;
  amrex::Real phi; pp.get("phi",phi);
  amrex::Real egr; pp.get("egr",egr);
  
  amrex::Real molefrac[NUM_SPECIES];

  for (int n = 0; n < NUM_SPECIES; n++){
    massfrac[n] = 0.0;
    molefrac[n] = 0.0;
  }

  // number of carbon and hydrogen atoms in each fuel component
  //order: IC8H18; NC7H16_ID; C6H5CH3; C2H5OH
  amrex::Real C_atoms[4] = {8.0 ,7.0 ,7.0,2.0}; 
  amrex::Real H_atoms[4] = {18.0,16.0,4.0,6.0}; 
  amrex::Real O_atoms[4] = {0.0 ,0.0 ,0.0,1.0}; 

  amrex::Real volfrac_IC8H18  = 0.284;
  amrex::Real volfrac_NC7H16  = 0.172;
  amrex::Real volfrac_C6H5CH3 = 0.339;
  amrex::Real volfrac_C2H5OH  = 0.205;

  amrex::Real volfrac_CO2 =  C_atoms[0]*volfrac_IC8H18 + C_atoms[1]*volfrac_NC7H16 + C_atoms[2]*volfrac_C6H5CH3 + C_atoms[3]*volfrac_C2H5OH;
  amrex::Real volfrac_H2O = (H_atoms[0]*volfrac_IC8H18 + H_atoms[1]*volfrac_NC7H16 + H_atoms[2]*volfrac_C6H5CH3 + H_atoms[3]*volfrac_C2H5OH)/2.0;

  amrex::Real volfrac_O2  = (2.0*volfrac_CO2+volfrac_H2O-(O_atoms[0]*volfrac_IC8H18  + 
                                                         O_atoms[1]*volfrac_NC7H16  + 
                                                         O_atoms[2]*volfrac_C6H5CH3 + 
                                                         O_atoms[3]*volfrac_C2H5OH))/2.0/phi;
  amrex::Real volfrac_N2  = volfrac_O2*(79./21.);

  amrex::Real sum_oxi  = volfrac_IC8H18 + volfrac_NC7H16 + volfrac_C6H5CH3 + volfrac_C2H5OH + volfrac_O2  + volfrac_N2;
  amrex::Real sum_prod = volfrac_CO2    + volfrac_H2O    + volfrac_N2;

  amrex::Real X_IC8H18  = volfrac_IC8H18   / sum_oxi;
  amrex::Real X_NC7H16  = volfrac_NC7H16   / sum_oxi;
  amrex::Real X_C6H5CH3 = volfrac_C6H5CH3  / sum_oxi;
  amrex::Real X_C2H5OH  = volfrac_C2H5OH   / sum_oxi;
  amrex::Real X_O2      = volfrac_O2       / sum_oxi;
  amrex::Real X_N2      = volfrac_N2       / sum_oxi;

  amrex::Real X_co2    = volfrac_CO2 / sum_prod;
  amrex::Real X_h2o    = volfrac_H2O / sum_prod;

  molefrac[IC8_ID]     = X_IC8H18  * (1.0-egr);
  molefrac[NC7H16_ID]  = X_NC7H16  * (1.0-egr);
  molefrac[C6H5CH3_ID] = X_C6H5CH3 * (1.0-egr);
  molefrac[C2H5OH_ID]  = X_C2H5OH  * (1.0-egr);
  molefrac[O2_ID]      = X_O2      * (1.0-egr);
  molefrac[N2_ID]      = X_N2      * (1.0-egr) + (egr*volfrac_N2 / sum_prod);
  molefrac[CO2_ID]     = X_co2*egr;
  molefrac[H2O_ID]     = X_h2o*egr;

  auto eos = pele::physics::PhysicsType::eos();
  eos.X2Y(molefrac,massfrac);
}

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
    // "denseAJ_direct"
    std::string solve_type_str = "none";

    int ode_ncells = 1;
    int ode_iE = 1;
    // pp.get("chem_integrator", chem_integrator);
    // pele::physics::reactions::ReactorBase::create(chem_integrator);
    std::unique_ptr<pele::physics::reactions::ReactorBase> reactor =
      pele::physics::reactions::ReactorBase::create(chem_integrator);
    {
      reactor->init(ode_iE, ode_ncells);
    }

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
    //Vector<std::unique_ptr<MultiFab>> outdata(1);
    int ngrow = 0;
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
    const int nCompOut = NUM_SPECIES;

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
    for (int i=0; i<amrData.PlotVarNames().size(); ++i)
    {
      if (amrData.PlotVarNames()[i] == spName)    idYlocal = i;
      if (amrData.PlotVarNames()[i] == TName)     idTlocal = i;
      if (amrData.PlotVarNames()[i] == RHOName)   idRlocal = i;
      if (amrData.PlotVarNames()[i] == vFracName) idvFlocal = i;
    }

    Print() << "Temperature index: " << idTlocal  << std::endl;
    Print() << "Density index    : " << idRlocal  << std::endl;
    Print() << "Vfrac index      : " << idvFlocal << std::endl;
    Print() << "First spec index : " << idYlocal  << std::endl;

    Vector<std::string> outNames(NUM_SPECIES);
    for (int i=0; i<NUM_SPECIES; ++i)
    {
      outNames[i] = "rho_omega_" + spec_names[i];
    }

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

      //vector containing the number of components to be read
      Vector<int> destFillComps(ncomp);
      for (int i=0; i<ncomp; ++i) destFillComps[i] = i;

      //reset MultiFab pointer
      indata[lev].reset(new MultiFab(ba,dm,ncomp,ngrow));
      Print() << "Reading data..." << std::endl; 
      amrData.FillVar(*indata[lev],lev,amrData.PlotVarNames(),destFillComps);
      Print() << "Done!" << std::endl;
      if (lev == Nlev - 1) {
        outdata.define(ba,dm,ncomp,ngrow);

        Geometry geom(amrData.ProbDomain()[lev],&rb,coord,&(is_per[0]));
        
        for (MFIter mfi(*indata[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
          const Box& box = mfi.tilebox();
          amrex::Array4<const amrex::Real> sfab = (*indata[lev]).array(mfi);
          amrex::Array4<Real> const output = outdata.array(mfi);
          // auto& fab = (*indata[0])[mfi];

          // amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
          AMREX_PARALLEL_FOR_3D ( box, i, j, k,
          {

          for (int n=0; n<amrData.PlotVarNames().size(); ++n){
            output(i,j,k,n) = sfab(i,j,k,n);
          }

          });

          // process(BL_TO_FORTRAN_BOX(box),
          //         BL_TO_FORTRAN_BOX(domain),
          //         BL_TO_FORTRAN_ANYD(fab),
          //         &ncomp,geom.ProbLo(),geom.CellSize());
        }
      }
    }

    std::string outfile(output_folder+getFileRoot(plotFileName) + "_finest");
    Print() << "Writing new data to " << outfile << std::endl;
    bool verb = false;

    // This was working in 3D but is not in 2D
    // WritePlotFile(GetVecOfPtrs(outdata),amrData,outfile,verb,outNames_test);
    
    //this version is working in 2D    
    int ref_ratio = 2;
//     Vector<Geometry> geoms(1);
//     Vector<int> levelSteps(1);
    Geometry geoms;
    int levelSteps;
    {
      geoms = Geometry(amrData.ProbDomain()[Nlev - 1],&rb,coord,&(is_per[0]));
      levelSteps = 0;
    }

    WriteSingleLevelPlotfile(outfile, outdata, amrData.PlotVarNames(), geoms, amrData.Time(), levelSteps);

  }
  Finalize();
  return 0;
}
