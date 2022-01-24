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
const bool verbose_DEF   = false;

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
clip_normalize_Y(amrex::Real* Y)
{

  // Clip
  amrex::Real sum = 0.0;
  for (int n = 0; n < NUM_SPECIES; n++) {
    Y[n] = amrex::max(amrex::min(Y[n], 1.0), 0.0);
    sum += Y[n];
  }

  // Normalize
  for (int n = 0; n < NUM_SPECIES; n++) {
    Y[n] /= sum;
  }
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
    int ncomp_all = amrData.NComp()+NUM_SPECIES;
    Print() << "Number of scalars in the original file: " << ncomp << std::endl;
    Print() << "Number of scalars in output file: " << ncomp_all << std::endl;

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
      // Print() << "domain " << domain <`< std::endl; 
      // Print() << "ba " << ba << std::endl; 
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
        outdata.define(ba,dm,ncomp_all,ngrow);

        Geometry geom(amrData.ProbDomain()[lev],&rb,coord,&(is_per[0]));
        // Print() << "amrData.ProbDomain()[lev] " << amrData.ProbDomain()[lev] << std::endl;
        
        for (MFIter mfi(*indata[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
          const Box& box = mfi.tilebox();
          amrex::Array4<const amrex::Real> sfab = (*indata[lev]).array(mfi);
          amrex::Array4<Real> const output = outdata.array(mfi);
          // amrex::Array4<Real> const WDOT = outdata.array(mfi);

          AMREX_PARALLEL_FOR_3D ( box, i, j, k,
          {

            //output PeleC scalars without computing anything
            for (int n=0; n<amrData.PlotVarNames().size(); ++n){
              output(i,j,k,n) = sfab(i,j,k,n);
            }

            //PeleLM processing requires to calculate reaction rates
            if(plt_src == 1){//PeleLM
              if(idvFlocal> -1){ //Dealing with EB cases
                if(sfab(i,j,k,idvFlocal) > 1.e-10){ //only compute reaction rate for volume fraction greater than zero
                  amrex::Real Yl[NUM_SPECIES];
                  for (int n=0; n<NUM_SPECIES; ++n) {
                    Yl[n] = sfab(i,j,k,idYlocal+n);
                  }
                  
                  amrex::Real reaction_rate[NUM_SPECIES];
                  amrex::Real rho = sfab(i,j,k,idRlocal)*convert_rho;
                  // printf("rho = %e  T = %e  \n", rho,sfab(i,j,k,idTlocal));
                  eos.RTY2WDOT(rho, sfab(i,j,k,idTlocal), Yl, reaction_rate);
                  for (int n=0; n<NUM_SPECIES; ++n) {
                    output(i,j,k,n+ncomp) = reaction_rate[n];
                    // if(n == O2_ID) printf("%e\n", WDOT(i,j,k,n));
                  }
                }
                else{//setting reaction rate to zero at vfrac = 0
                  for (int n=0; n<NUM_SPECIES; ++n) {
                    output(i,j,k,n) = 0.0;
                  }          
                }
              }

              else{ //The code below is used for cases without EB
                  amrex::Real Yl[NUM_SPECIES];
                  amrex::Real rho_from_Y=0.0;
                  amrex::Real rho = sfab(i,j,k,idRlocal)*convert_rho;
                  for (int n=0; n<NUM_SPECIES; ++n) {
                    Yl[n] = sfab(i,j,k,idYlocal+n);
                    rho_from_Y += Yl[n]*rho;
                  }
                  clip_normalize_Y(Yl);
                  
                  amrex::Real reaction_rate[NUM_SPECIES];
                  // printf("rho = %e  T = %e  \n", rho,sfab(i,j,k,idTlocal));
                  eos.RTY2WDOT(rho_from_Y, sfab(i,j,k,idTlocal), Yl, reaction_rate);
                  for (int n=0; n<NUM_SPECIES; ++n) {
                    output(i,j,k,n+ncomp) = reaction_rate[n];
                  }    
              }
            }
          });

          // process(BL_TO_FORTRAN_BOX(box),
          //         BL_TO_FORTRAN_BOX(domain),
          //         BL_TO_FORTRAN_ANYD(fab),
          //         &ncomp,geom.ProbLo(),geom.CellSize());
        }
      }
    }

    int ref_ratio = 2;
    Geometry geoms;
    int levelSteps;

    {
      geoms = Geometry(amrData.ProbDomain()[Nlev - 1],&rb,coord,&(is_per[0]));
      // geoms = Geometry(subbox,&rb,coord,&(is_per[0]));
      levelSteps = 0;
    }

    Box small_box = amrData.boxArray(Nlev - 1).minimalBox();

    BoxArray ba_sub(small_box);
    ba_sub.maxSize(max_grid_size);
    DistributionMapping dmap_sub(ba_sub);

    Print() << "Box small_box: " << small_box << std::endl; 
    Print() << "Number of boxes in subset: " << ba_sub.size() << std::endl; 

    // Vector<MultiFab> data_sub(1);
    MultiFab data_sub;
    data_sub.define(ba_sub,dmap_sub,ncomp_all,ngrow);
    data_sub.ParallelCopy(outdata,0,0,ncomp_all);

    Vector<std::string> outNames(ncomp_all);
    for (int i=0; i<amrData.PlotVarNames().size(); ++i)
    {
      outNames[i] = amrData.PlotVarNames()[i];
    }
    if(plt_src == 1){ //Only PeleLM requires to compute rho*omega_i
      for (int i=0; i<NUM_SPECIES; ++i)
      {
        outNames[i+ncomp] = "rho_omega_" + spec_names[i];
      }
    }

    std::string outfile(output_folder+getFileRoot(plotFileName) + "_finest_small");
    Print() << "Writing new data to " << outfile << std::endl;
    WriteSingleLevelPlotfile(outfile, data_sub, outNames, geoms, amrData.Time(), levelSteps);

  }
  Finalize();
  return 0;
}