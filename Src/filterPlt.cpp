#include <string>
#include <iostream>

#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include <AMReX_DataServices.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_BLFort.H>
#include <Filter.H>
#include <PltFileManager.H>

static
void
print_usage (int,
             char* argv[])
{
  std::cerr << "Utility to average pltfiles on same domain but with non-matching AMR";
  std::cerr << "usage:\n";
  std::cerr << argv[0] << " infile=<s> [options] \n\tOptions:\n";
  std::cerr << "\t     infile=<s> where s is a pltfile \n";
  std::cerr << "\t     variables=<s1 s2 s3> where <s1> <s2> and <s3> are variable names to select for combined pltfile [DEF-> all possible]\n";
  std::cerr << "\t     max_filter_level=<int> where <int> is the max refinement level to filter, zero-indexed [DEF->1000]\n";
  std::cerr << "\t     filter_type=<int> where <int> is the filter type as defined in PeleC (1->box, 2->Gaussian, etc) [DEF->1000]\n";
  std::cerr << "\t     base_fgr=<int> where <int> is the desired filter to grid ratio on the base level, must be even [DEF->2]\n";
  std::cerr << "\t     same_fgr_all_levels=<bool> where if true the same filter to grid ratio is kept on all levels (rather than absolute filter width) [DEF->false]\n";
  std::cerr << "\t     max_grid_size=<int> where <int> is the AMReX max_grid_size for the output [DEF->32]\n";
exit(1);
}

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

    if (argc < 2 || pp.contains("help")) {
       print_usage(argc,argv);
     }

    std::string infile        = "";
    int finestLevel           = 1000;
    amrex::Vector<Filter> les_filter;
    int les_filter_type = 1;
    int les_filter_fgr  = 2;
    bool same_fgr_all_levels = false;

    pp.get("infile",infile);
    pp.query("max_filter_level",finestLevel);
    pp.query("filter_type",les_filter_type);
    pp.query("base_fgr",les_filter_fgr); // filter to grid ratio
    pp.query("same_fgr_all_levels", same_fgr_all_levels);
    int max_grid_size = 32;
    pp.query("max_grid_size",max_grid_size);

    // use PltFileManager to load data
    amrex::Vector<pele::physics::pltfilemanager::PltFileManager*> plt_file_data(1);
    plt_file_data[0] = new pele::physics::pltfilemanager::PltFileManager(infile);
    // Plotfile global infos
    int Nlev = std::min(finestLevel + 1, plt_file_data[0]->getNlev());

    // Variable names
    int ncomp_filter;
    amrex::Vector<std::string> variableNames;
    amrex::Vector<int> var_idxs;
    int nvar = pp.countval("variables");
    bool all_vars = nvar > 0 ? false : true;
    if (!all_vars) {
      ncomp_filter = nvar;
      pp.getarr("variables",variableNames);
      const amrex::Vector<std::string>& plotVarNames = plt_file_data[0]->getVariableList();
      for (int var = 0; var < nvar; ++var) {
        int pvar;
        for (pvar = 0; pvar < plotVarNames.size(); ++pvar) {
          if (variableNames[var] == plotVarNames[pvar]) {
            var_idxs.push_back(pvar);
            break;
          }
        }
        if (pvar == plotVarNames.size()) {
          amrex::Abort("Variable '" + variableNames[var] + "' not found in file");
        }
      }
    } else {
      variableNames = plt_file_data[0]->getVariableList();
      ncomp_filter = variableNames.size();
    }

    amrex::Vector<amrex::MultiFab> indata(Nlev);
    amrex::Vector<amrex::MultiFab> outdata(Nlev);
    amrex::Vector<amrex::Geometry> level_geometries;

    amrex::Print() << "Reading data..." << std::endl;
    int les_filter_fgr_lev = les_filter_fgr;
    for (int lev=0; lev<Nlev; ++lev)
    {
      amrex::Print() << "on level "<< lev << std::endl;

      // Initialize filter stuff
      if ((!same_fgr_all_levels) && (lev > 0)) {
        les_filter_fgr_lev *= plt_file_data[0]->getRefRatio(lev - 1);
      }

      les_filter.push_back(Filter(les_filter_type, les_filter_fgr_lev));
      int nGrowF = les_filter[lev].get_filter_ngrow();
      int nGrow  = 0;

      amrex::BoxArray ba(plt_file_data[0]->getGrid(lev));
      ba.maxSize(max_grid_size);
      const amrex::DistributionMapping dm = amrex::DistributionMapping(ba);
      level_geometries.push_back(plt_file_data[0]->getGeom(lev));
      indata[lev].define(ba,dm,ncomp_filter,nGrowF);
      outdata[lev].define(ba,dm,ncomp_filter,nGrow);

      if (all_vars) {
        plt_file_data[0]->fillPatchFromPlt(lev, level_geometries[lev], 0, 0, ncomp_filter, indata[lev], indata[lev].nGrowVect());
      } else {
        for (int var = 0; var < nvar; ++var) {
          plt_file_data[0]->fillPatchFromPlt(lev, level_geometries[lev], var_idxs[var], var, 1, indata[lev], indata[lev].nGrowVect());
        }
      }
    }
    amrex::Print() << "Done!" << std::endl;

    amrex::Print() << "Filtering data..." << std::endl;
    for (int lev=0; lev<Nlev; ++lev)
    {
      amrex::Print() << "on level " << lev << std::endl;

      for (amrex::MFIter mfi(indata[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {

        amrex::FArrayBox& fab_in  = (indata[lev])[mfi];
        amrex::FArrayBox& fab_out = (outdata[lev])[mfi];
        const amrex::Box& box = mfi.validbox();

        les_filter[lev].apply_filter(box, fab_in, fab_out);
      }
    }
    amrex::Print() << "Done!" << std::endl;

    amrex::Print() << "Saving filtered data..." << std::endl;
    std::string outfile(getFileRoot(infile) + "_filtered");

    write_plotfile(outfile,Nlev,amrex::GetVecOfConstPtrs(outdata),variableNames,level_geometries,plt_file_data[0]->getTime(),amrex::Vector<int>(Nlev, 0));
    amrex::Print() << "Done!" << std::endl;

}
