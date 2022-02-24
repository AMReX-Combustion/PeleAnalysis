#include <string>
#include <iostream>
#include <set>

#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include <AMReX_DataServices.H>
#include <AMReX_WritePlotFile.H>
#include <AMReX_PlotFileUtil.H>

#include "mechanism.H"
#include <PelePhysics.H>
#include <ReactorBase.H>

using namespace amrex;

extern "C" {
  void pushvtog(const int* lo,  const int* hi,
                const int* dlo, const int* dhi,
                Real* U, const int* Ulo, const int* Uhi,
                const int* nc);

  void compute_Kolmogorov(const int* lo,  const int* hi,
                const Real* U , const int* Ulo  , const int* Uhi,
                const Real* Um, const int* Um_lo, const int* Um_hi,
                Real* eta, const int* eta_lo, const int* eta_hi,
                const Real* dx);
}

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
    // vector<std::string> tokens = Tokenize(infile,std::string("/"))
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

    std::string chem_integrator = "ReactorNull";
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

    Vector<std::unique_ptr<MultiFab>> indata_ensemble(1);
    Vector<std::unique_ptr<MultiFab>> indata_snapshot(1);
    Vector<std::unique_ptr<MultiFab>> outdata(1);
    Vector<Geometry*> geoms(1);
    //MultiFab outdata;

    auto eos = pele::physics::PhysicsType::eos();
    // pele::physics::eos::speciesNames(spec_names);


    //std::string output_folder; pp.get("outFolder",output_folder);
    // const int number_of_plts = pp.countval("plt_files");
    std::string ensemble_file; pp.get("ensemble_file",ensemble_file);
    std::string snapshot_file; pp.get("snapshot_file",snapshot_file);

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
    //MultiFab outdata;
    MultiFab output_storage;

    //Vector<std::unique_ptr<MultiFab>> outdata(1);
    int ngrow = 1;
    int ncomp = amrData_en.NComp();
    Print() << "Number of scalars in the original file: " << ncomp << std::endl;

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
    outNames[0] = "eta";
    outNames[1] = "epsilon";
    outNames[2] = "gradx_dx";

    RealBox rb(&(amrData_en.ProbLo()[0]),
               &(amrData_en.ProbHi()[0]));
    Vector<int> is_per(AMREX_SPACEDIM,0);

    int coord = 0;
    int max_grid_size = 32;
    pp.query("max_grid_size",max_grid_size);

    for (int lev=0; lev<Nlev; ++lev)
    {        
      Print() << "Working on level " << lev << std::endl; 
      //define domain stuff
      const auto& domain = amrData_en.ProbDomain()[lev];
      BoxArray ba(domain);
      ba.maxSize(max_grid_size);
      const DistributionMapping dm(ba);
      const Vector<Real>& delta = amrData_en.DxLevel()[lev];

      Print() << "Number of boxes: " << ba.size() << std::endl; 

      //vector containing the number of components to be read
      Vector<int> destFillComps(ncomp);
      for (int i=0; i<ncomp; ++i) destFillComps[i] = i;

      //reset MultiFab pointer
      indata_ensemble[lev].reset(new MultiFab(ba,dm,ncomp,ngrow));
      Print() << "Reading ensemble data..." << std::endl; 
      amrData_en.FillVar(*indata_ensemble[lev],lev,amrData_en.PlotVarNames(),destFillComps);
      Print() << "Done!" << std::endl;

      int ncomp_sn = amrData_sn.NComp();
      Vector<int> destFill(ncomp_sn);
      for (int i=0; i<ncomp_sn; ++i) destFill[i] = i;

      indata_snapshot[lev].reset(new MultiFab(ba,dm,ncomp_sn,ngrow));
      Print() << "Reading snapshot data..." << std::endl; 
      amrData_en.FillVar(*indata_snapshot[lev],lev,amrData_en.PlotVarNames(),destFill);
      Print() << "Done!" << std::endl;

      if (lev == Nlev - 1) {
        //outdata.define(ba,dm,nCompOut,ngrow);
        output_storage.define(ba,dm,nCompOut,ngrow);
	int ngrow_out=0;
        outdata[lev].reset(new MultiFab(ba,dm,nCompOut,ngrow_out));

         // Fix up grow cells.  Use extrap for guess
         const Box& dbox = amrData_en.ProbDomain()[lev];
         for (MFIter mfi(*indata_ensemble[lev]); mfi.isValid(); ++mfi)
         {
           FArrayBox& fab = (*indata_ensemble[lev])[mfi];
           const Box& box = mfi.validbox();
           pushvtog(BL_TO_FORTRAN_BOX(box),
                    BL_TO_FORTRAN_BOX(dbox),
                    BL_TO_FORTRAN_ANYD(fab),
                    &ncomp);
         }

         // Fix up grow cells.  Use extrap for guess
         // const Box& dbox = amrData_en.ProbDomain()[lev];
         for (MFIter mfi(*indata_snapshot[lev]); mfi.isValid(); ++mfi)
         {
           FArrayBox& fab = (*indata_snapshot[lev])[mfi];
           const Box& box = mfi.validbox();
           pushvtog(BL_TO_FORTRAN_BOX(box),
                    BL_TO_FORTRAN_BOX(dbox),
                    BL_TO_FORTRAN_ANYD(fab),
                    &ncomp);
         }    

        Geometry geom(amrData_en.ProbDomain()[lev],&rb,coord,&(is_per[0]));
        geoms[0] = new Geometry(amrData_en.ProbDomain()[lev],&rb,coord,&(is_per[0]));

        //// Fix up fine-fine and periodic
        indata_ensemble[lev]->FillBoundary(geoms[0]->periodicity());
        indata_snapshot[lev]->FillBoundary(geoms[0]->periodicity());

        Print() << "Starting MFIter..." << std::endl;       
        for (MFIter mfi(*indata_ensemble[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
          const Box& box = mfi.tilebox();
          // amrex::Array4<const amrex::Real> sfab = (*indata_ensemble[lev]).array(mfi);

          FArrayBox& fab_en = (*indata_ensemble[lev])[mfi];
          FArrayBox& fab_sn = (*indata_snapshot[lev])[mfi];
          FArrayBox& output = (*outdata[lev])[mfi];

          // amrex::Array4<Real> const output = outdata.array(mfi);
          // auto& fab = (*indata_ensemble[0])[mfi];

          //FArrayBox& fab = (*state[lev])[mfi];
          //const Box& box = mfi.validbox();
          compute_Kolmogorov(BL_TO_FORTRAN_BOX(box),
                   BL_TO_FORTRAN_N_ANYD(fab_en,idUlocal_en),
                   BL_TO_FORTRAN_N_ANYD(fab_sn,idUlocal_sn),
                   BL_TO_FORTRAN_N_ANYD(output,0),
                   // BL_TO_FORTRAN_ANYD(output),
                   &(delta[0]));


          // process(BL_TO_FORTRAN_BOX(box),
          //         BL_TO_FORTRAN_BOX(domain),
          //         BL_TO_FORTRAN_ANYD(fab),
          //         &ncomp,geom.ProbLo(),geom.CellSize());
        }
      }
    }


      std::string outfile("plt_turbulence");
      Print() << "Writing new data to " << outfile << std::endl;
      bool verb = false;

      // This was working in 3D but is not in 2D
      WritePlotFile(GetVecOfPtrs(outdata),amrData_sn,outfile,verb,outNames);
      
//      //this version is working in 2D    
//      int ref_ratio = 2;
//  //     Vector<Geometry> geoms(1);
//  //     Vector<int> levelSteps(1);
//      Geometry geoms;
//      int levelSteps;
//      {
//        geoms = Geometry(amrData_en.ProbDomain()[Nlev - 1],&rb,coord,&(is_per[0]));
//        levelSteps = 0;
//      }
//
//      output_storage.copy(outdata,0,0,nCompOut);
//      WriteSingleLevelPlotfile(outfile, output_storage, outNames, geoms, amrData_sn.Time(), levelSteps);
//
  }
  Finalize();
  return 0;
}
