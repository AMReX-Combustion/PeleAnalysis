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
 //   std::unique_ptr<pele::physics::reactions::ReactorBase> reactor =
 //     pele::physics::reactions::ReactorBase::create(chem_integrator);
 //   {
 //     reactor->init(ode_iE, ode_ncells);
 //   }

    amrex::sundials::Initialize(amrex::OpenMP::get_max_threads());

    std::unique_ptr<pele::physics::reactions::ReactorBase> reactor =
      pele::physics::reactions::ReactorBase::create(chem_integrator);
    reactor->init(ode_iE, ode_ncells);

    
    Vector<std::unique_ptr<MultiFab>> indata_ensemble(1);
    Vector<std::unique_ptr<MultiFab>> indata_snapshot(1);

    Vector<std::string> spec_names;
    auto eos = pele::physics::PhysicsType::eos();

    std::string ensemble_file; pp.get("ensemble_file",ensemble_file);
    std::string snapshot_file; pp.get("snapshot_file",snapshot_file);
    std::string output_folder; pp.get("outFolder",output_folder);

    // std::string plotFileName;

    DataServices::SetBatchMode();
    Amrvis::FileType fileType(Amrvis::NEWPLT);
    
    Print() << "Processing ensemble from " << ensemble_file << std::endl;
    Print() << "Processing snapshot from " << snapshot_file << std::endl;

    DataServices dataServices_en(ensemble_file, fileType);
    if( ! dataServices_en.AmrDataOk()) {
      DataServices::Dispatch(DataServices::ExitRequest, NULL);
      // ^^^ this calls ParallelDescriptor::EndParallel() and exit()
    }
    AmrData& amrData_en = dataServices_en.AmrDataRef();

    DataServices dataServices_sn(snapshot_file, fileType);
    if( ! dataServices_sn.AmrDataOk()) {
      DataServices::Dispatch(DataServices::ExitRequest, NULL);
      // ^^^ this calls ParallelDescriptor::EndParallel() and exit()
    }
    AmrData& amrData_sn = dataServices_sn.AmrDataRef();

    int finestLevel = amrData_en.FinestLevel();
    pp.query("finestLevel",finestLevel);
    int Nlev = finestLevel + 1;
    Print() << "Working with " << Nlev << " levels" << std::endl;

    //Vector<std::unique_ptr<MultiFab>> indata_ensemble(Nlev);
    MultiFab outdata;
    // MultiFab output_storage;

    //Vector<std::unique_ptr<MultiFab>> outdata(1);
    int ngrow = 0;
    int ncomp = amrData_en.NComp();
    int ncomp_sn = amrData_sn.NComp();
    Print() << "Number of scalars in the ensemble file: " << ncomp << std::endl;
    Print() << "Number of scalars in the snapshot file: " << ncomp_sn << std::endl;

    // const auto& domain = amrData_en.ProbDomain()[finestLevel];
    // BoxArray ba(domain);
    // int max_grid_size = 32;
    // pp.query("max_grid_size",max_grid_size);
    // ba.maxSize(max_grid_size);
    // const DistributionMapping dm(ba);
    
    int idUlocal_en  = -1; // Velocity index
    int idUlocal_sn  = -1; // Velocity index

    const std::string UName = "x_velocity";

    for (int i=0; i<amrData_en.PlotVarNames().size(); ++i)
    {
      if (amrData_en.PlotVarNames()[i] == UName)    idUlocal_en = i;
    }

    for (int i=0; i<amrData_sn.PlotVarNames().size(); ++i)
    {
      if (amrData_sn.PlotVarNames()[i] == UName)    idUlocal_sn = i;
    }

    Print() << "x_velocity index in ensemble file: " << idUlocal_en  << std::endl;
    Print() << "x_velocity index in snapshot file: " << idUlocal_sn  << std::endl;

    int nCompOut = 3;
    Vector<std::string> outNames(nCompOut);
    outNames[0] = "ux_prime";
    outNames[1] = "uy_prime";
    outNames[2] = "uz_prime";


    RealBox rb_en(&(amrData_en.ProbLo()[0]),
               &(amrData_en.ProbHi()[0]));
    RealBox rb_sn(&(amrData_sn.ProbLo()[0]),
               &(amrData_sn.ProbHi()[0]));
    Vector<int> is_per(AMREX_SPACEDIM,0);

    int coord = 0;
    int max_grid_size = 32;
    pp.query("max_grid_size",max_grid_size);

    for (int lev=0; lev<Nlev; ++lev)
    {        
      Print() << "Working on level " << lev << std::endl; 
      //define domain stuff
      const auto& domain_en = amrData_en.ProbDomain()[lev];
      BoxArray ba_en(domain_en);
      ba_en.maxSize(max_grid_size);
      const DistributionMapping dm_en(ba_en);
      const Vector<Real>& delta = amrData_en.DxLevel()[lev];

      Print() << "Number of boxes in ensemble: " << ba_en.size() << std::endl; 

      //vector containing the number of components to be read
      Vector<int> destFillComps(ncomp);
      for (int i=0; i<ncomp; ++i) destFillComps[i] = i;

      //reset MultiFab pointer
      indata_ensemble[lev].reset(new MultiFab(ba_en,dm_en,ncomp,ngrow));
      Print() << "Reading ensemble data..." << std::endl; 
      amrData_en.FillVar(*indata_ensemble[lev],lev,amrData_en.PlotVarNames(),destFillComps);
      Print() << "Done!" << std::endl;

      const auto& domain_sn = amrData_sn.ProbDomain()[lev];
      BoxArray ba_sn(domain_sn);
      ba_sn.maxSize(max_grid_size);
      const DistributionMapping dm_sn(ba_sn);

      Print() << "Number of boxes in snapshot: " << ba_sn.size() << std::endl; 

      Vector<int> destFill(ncomp_sn);
      for (int i=0; i<ncomp_sn; ++i) destFill[i] = i;

      indata_snapshot[lev].reset(new MultiFab(ba_sn,dm_sn,ncomp_sn,ngrow));
      Print() << "Reading snapshot data..." << std::endl; 

      amrData_sn.FillVar(*indata_snapshot[lev],lev,amrData_sn.PlotVarNames(),destFill);
      Print() << "Done!" << std::endl;
      if (lev == Nlev - 1) {

        BoxArray ba_sub(domain_sn);
        ba_sub.maxSize(max_grid_size);
        DistributionMapping dmap_sub(ba_sub);

        //copy into a new MultiFab common to both ensemble and snapshot data
        MultiFab data_ensemble;
        int ensemble_size = amrData_en.PlotVarNames().size();
        data_ensemble.define(ba_sub,dmap_sub,ensemble_size,ngrow);
        data_ensemble.ParallelCopy(*indata_ensemble[lev],0,0,ensemble_size);

        MultiFab data_snapshot;
        int snapshot_size = amrData_sn.PlotVarNames().size();
        data_snapshot.define(ba_sub,dmap_sub,snapshot_size,ngrow);
        data_snapshot.ParallelCopy(*indata_snapshot[lev],0,0,snapshot_size);

        outdata.define(ba_sub,dmap_sub,nCompOut,ngrow);


        Geometry geom(amrData_sn.ProbDomain()[lev],&rb_sn,coord,&(is_per[0]));
        // Print() << "amrData.ProbDomain()[lev] " << amrData.ProbDomain()[lev] << std::endl;
        
        for (MFIter mfi(data_ensemble,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
          const Box& box = mfi.tilebox();
          // amrex::Array4<const amrex::Real> sfab_ensemble = (*indata_ensemble[lev]).array(mfi);
          amrex::Array4<const amrex::Real> sfab_snapshot = (*indata_snapshot[lev]).array(mfi);

          amrex::Array4<Real> const sfab_ensemble = data_ensemble.array(mfi);
          amrex::Array4<Real> const output = outdata.array(mfi);
          // amrex::Array4<Real> const WDOT = outdata.array(mfi);

          AMREX_PARALLEL_FOR_3D ( box, i, j, k,
          {

            //output PeleC scalars without computing anything
            for (int n=0; n<nCompOut; ++n){
              output(i,j,k,n) = sfab_snapshot(i,j,k,n+idUlocal_sn) - sfab_ensemble(i,j,k,n+idUlocal_en);
              // output(i,j,k,n) = sfab_ensemble(i,j,k,n+idUlocal_sn);
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
      geoms = Geometry(amrData_en.ProbDomain()[Nlev - 1],&rb_sn,coord,&(is_per[0]));
      // geoms = Geometry(subbox,&rb,coord,&(is_per[0]));
      levelSteps = 0;
    }

    // Box small_box = amrData_en.boxArray(Nlev - 1).minimalBox();

    // BoxArray ba_sub(small_box);
    // ba_sub.maxSize(max_grid_size);
    // DistributionMapping dmap_sub(ba_sub);

    // Print() << "Box small_box: " << small_box << std::endl; 
    // Print() << "Number of boxes in subset: " << ba_sub.size() << std::endl; 

    // // Vector<MultiFab> data_sub(1);
    // MultiFab data_sub;
    // data_sub.define(ba_sub,dmap_sub,ncomp_all,ngrow);
    // data_sub.ParallelCopy(outdata,0,0,ncomp_all);


    std::string outfile(output_folder+getFileRoot(snapshot_file) + "_uprime");
    Print() << "Writing new data to " << outfile << std::endl;
    WriteSingleLevelPlotfile(outfile, outdata, outNames, geoms, amrData_sn.Time(), levelSteps);

  }
  Finalize();
  return 0;
}
