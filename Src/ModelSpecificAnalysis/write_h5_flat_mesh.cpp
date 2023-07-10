#include <string>
#include <iostream>
#include <set>

#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include <AMReX_DataServices.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_WritePlotFile.H>
#include <AMReX_Sundials.H>

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


int
main (int   argc,
      char* argv[])
{
  Initialize(argc,argv);
  {
    if (argc < 2)
      print_usage(argc,argv);
    
    // ---------------------------------------------------------------------
    // Set defaults input values
    // ---------------------------------------------------------------------
    Vector<std::string> VarsOut;
    Vector<int> scalar_index;
    // ---------------------------------------------------------------------


    ParmParse pp;

    if (pp.contains("help"))
      print_usage(argc,argv);

    if (pp.contains("verbose"))
      AmrData::SetVerbose(true);

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
    const Vector<std::string>& plotVarNames = amrData.PlotVarNames();
    Print() << "Number of scalars in the original file: " << ncomp << std::endl;

    int const nvars_user =  pp.countval("VarsOut");

    int nvars_out = 0;
    if(nvars_user == 0){
      Print() << "VarsOut not defined in input. Including all data from source plt file..." << std::endl;
      nvars_out = amrData.PlotVarNames().size();
      VarsOut.resize(nvars_out);
      scalar_index.resize(nvars_out);
      for (int i=0; i<plotVarNames.size(); ++i){
        VarsOut[i] = amrData.PlotVarNames()[i];
        scalar_index[i] = i;
      }
    }
    else{
        pp.queryarr("VarsOut",VarsOut,0,nvars_user);
        nvars_out = nvars_user;
        scalar_index.resize(nvars_out);
        // Check if variable exists
        int idC = -1;
        for (int iscal=0; iscal<nvars_out; ++iscal){
        {
          for (int i=0; i<plotVarNames.size(); ++i)
            if (plotVarNames[i] == VarsOut[iscal]) {
              idC = i;
              scalar_index[iscal] = i;
            }
          }
          if (idC<0) {
            Print() << "Cannot find " << VarsOut[iscal] << " data in pltfile \n";
            amrex::Abort();
          }
          idC = -1;
        }
    }

    Print() << "Number of scalars to be written in output: " << nvars_out << std::endl;

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

      if (lev == Nlev - 1) {
        //reset MultiFab pointer
        indata[lev].reset(new MultiFab(ba,dm,ncomp,ngrow));
        Print() << "Reading data..." << std::endl; 
        amrData.FillVar(*indata[lev],lev,amrData.PlotVarNames(),destFillComps);
        Print() << "Done!" << std::endl;
        outdata.define(ba,dm,nvars_out,ngrow);

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

          // for (int n=0; n<amrData.PlotVarNames().size(); ++n){
          //   output(i,j,k,n) = sfab(i,j,k,n);
          // }
            for (int n=0; n<nvars_out; ++n){
              output(i,j,k,n) = sfab(i,j,k,scalar_index[n]);
            }

          });

        }
      }
    }


    //Find minimal box size that encompass the finest level
    //This will be the new domain size
    Box small_box = amrData.boxArray(Nlev - 1).minimalBox();

    BoxArray ba_sub(small_box);
    ba_sub.maxSize(max_grid_size);
    DistributionMapping dmap_sub(ba_sub);
   
    //New MultiFab with finest level data in a domain encompassing only the finest level
    MultiFab outdata_sub;
    outdata_sub.define(ba_sub,dmap_sub,nvars_out,ngrow);
    outdata_sub.ParallelCopy(outdata,0,0,nvars_out);


    // std::string outfile(output_folder+getFileRoot(plotFileName) + "_finest");
    std::string outfile(getFileRoot(plotFileName) + "_finest");
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

 
    std::string compression = "None@0";
    int level_step = 0;
    WriteSingleLevelPlotfileHDF5 (outfile,
                                  outdata_sub,
                                  VarsOut,
                                  geoms,
                                  amrData.Time(),
                                  level_step,
                                  compression);

  }
  Finalize();
  return 0;
}
