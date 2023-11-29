/** \file unit-tests-main.cpp
 *  Entry point for unit tests
 */

#include <string>
#include <iostream>

#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include <AMReX_DataServices.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_BLFort.H>
#include <Filter.H>

std::string
getFileRoot(const std::string& infile)
{
  std::vector<std::string> tokens = amrex::Tokenize(infile,std::string("/"));
  return tokens[tokens.size()-1];
}


void write_plotfile (const std::string &plotfilename,
                              int nlevels,
                              const amrex::Vector<const amrex::MultiFab*> &mf,
                              const amrex::Vector<std::string> &varnames,
                              const amrex::Vector<amrex::Geometry> &geom,
                              amrex::Real time,
                              const amrex::Vector<int> &level_steps)
{

  amrex::Vector<amrex::IntVect> ref_ratio(nlevels-1,{AMREX_D_DECL(2, 2, 2)});
  amrex::WriteMultiLevelPlotfile(plotfilename, nlevels, mf, varnames,
                            geom, time, level_steps,ref_ratio);
}

int
main(int argc, char** argv)
{

    amrex::Initialize(argc,argv);
    // ---------------------------------------------------------------------
    // ParmParse
    // ---------------------------------------------------------------------
    amrex::ParmParse pp;

    std::string infile        = "";
    int finestLevel           = 1000;
    amrex::Vector<Filter> les_filter;
    int les_filter_type = 1;
    int les_filter_fgr  = 2;
    bool same_fgr_all_levels = false;
    int coord = 0;

    pp.get("infile",infile);
    pp.query("finestLevel",finestLevel);
    pp.query("filter_type",les_filter_type);
    pp.query("base_fgr",les_filter_fgr); // filter to grid ratio
    pp.query("same_fgr_all_levels", same_fgr_all_levels);
    int max_grid_size = 32;
    pp.query("max_grid_size",max_grid_size);

    // Initialize DataService
    amrex::DataServices::SetBatchMode();
    amrex::Amrvis::FileType fileType(amrex::Amrvis::NEWPLT);
    amrex::DataServices dataServices(infile, fileType);
    if( ! dataServices.AmrDataOk()) {
      amrex::DataServices::Dispatch(amrex::DataServices::ExitRequest, NULL);
    }
    amrex::AmrData& amrData = dataServices.AmrDataRef();

    // Plotfile global infos
    finestLevel = std::min(finestLevel,amrData.FinestLevel());
    int Nlev = finestLevel + 1;

    const amrex::Vector<std::string>& plotVarNames = amrData.PlotVarNames();
    amrex::RealBox rb(&(amrData.ProbLo()[0]),
               &(amrData.ProbHi()[0]));
    amrex::Vector<int> is_per(AMREX_SPACEDIM,0);

    int ncomp = amrData.NComp();
    amrex::Print() << "Number of scalars in the original file: " << ncomp << std::endl;


    std::string VxName  = "x_velocity";
    std::string VyName  = "y_velocity";
    std::string VzName  = "z_velocity";
    std::string rhoName = "density";
    amrex::Vector<std::string> compNames;
    int ncomp_filter = 4;
    compNames.resize(ncomp_filter);

    compNames[0] = "x_velocity";
    compNames[1] = "y_velocity";
    compNames[2] = "z_velocity";
    compNames[3] = "density";

    int idVx  = -1;
    int idVy  = -1;
    int idVz  = -1;
    int idRho = -1;

    for (int i=0; i<amrData.PlotVarNames().size(); ++i)
    {
      if (amrData.PlotVarNames()[i] == VxName)  idVx = i;
      if (amrData.PlotVarNames()[i] == VyName)  idVy = i;
      if (amrData.PlotVarNames()[i] == VzName)  idVz = i;
      if (amrData.PlotVarNames()[i] == rhoName) idRho = i;
    }

    amrex::Vector<std::unique_ptr<amrex::MultiFab>> indata(Nlev);
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> outdata(Nlev);

    amrex::Vector<amrex::Geometry> geoms(Nlev);

    amrex::Vector<int> destFillComps(ncomp_filter);
    for (int i=0; i<ncomp_filter; ++i) destFillComps[i] = i;

    int les_filter_fgr_lev = les_filter_fgr;
    for (int lev=0; lev<Nlev; ++lev)
    {
      // Initialize filter stuff
      if ((!same_fgr_all_levels) && (lev > 0)) {
        les_filter_fgr_lev *= amrData.RefRatio()[lev - 1];
      }

      les_filter[lev] = Filter(les_filter_type, les_filter_fgr_lev);
      int nGrowF = les_filter[lev].get_filter_ngrow();
      int nGrow  = 0;

      const auto& domain = amrData.ProbDomain()[lev];
      amrex::BoxArray ba(domain);
      ba.maxSize(max_grid_size);
      const amrex::DistributionMapping dm(ba);

      geoms[lev].define(amrData.ProbDomain()[lev], &rb, coord, &(is_per[0]));

      indata [lev].reset(new amrex::MultiFab(ba,dm,ncomp_filter,nGrowF));
      outdata[lev].reset(new amrex::MultiFab(ba,dm,ncomp_filter,nGrow));

      amrex::Print() << "Reading data..." << std::endl;
      amrData.FillVar(*indata[lev],lev,compNames,destFillComps);
      amrex::Print() << "Done!" << std::endl;
    }

    amrex::Print() << "Filtering data..." << std::endl;
    for (int lev=0; lev<Nlev; ++lev)
    {
      amrex::Print() << "on level "<< lev << std::endl;

      for (amrex::MFIter mfi(*indata[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {

        amrex::FArrayBox& fab_in  = (*indata[lev])[mfi];
        amrex::FArrayBox& fab_out = (*outdata[lev])[mfi];
        const amrex::Box& box = mfi.validbox();

        les_filter[lev].apply_filter(box, fab_in, fab_out);
      }
    }
    amrex::Print() << "Done!" << std::endl;

    amrex::Print() << "Saving filtered data..." << std::endl;
    std::string outfile(getFileRoot(infile) + "_filtered");
    //bool verb = false;

    //    int levelSteps = 0;

    write_plotfile(outfile,Nlev,amrex::GetVecOfConstPtrs(outdata),compNames,geoms,amrData.Time(),amrex::Vector<int>(Nlev, 0));
    amrex::Print() << "Done!" << std::endl;

}
