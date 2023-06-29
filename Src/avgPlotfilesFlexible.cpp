#include <AMReX_ParmParse.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_PlotFileUtil.H>
#include <PltFileManager.H>

using namespace amrex;

static
void
print_usage (int,
             char* argv[])
{
  std::cerr << "Utility to average pltfiles on same domain but with non-matching AMR";
  std::cerr << "usage:\n";
  std::cerr << argv[0] << "infiles=<s1 s2 s3> [options] \n\tOptions:\n";
  std::cerr << "\t     infiles=<s1 s2 s3> where <s1> <s2> amnd <s3> are pltfiles\n";
  std::cerr << "\t     outfile=<s> where <s> is the output pltfile\n";
  std::cerr << "\t     variables=<s1 s2 s3> where <s1> <s2> and <s3> are variable names to select for combined pltfile [DEF-> all possible]\n";
  std::cerr << "\t     output_max_level=<s> where <s> is the max refinement level to combine, zero-indexed [DEF->1000]\n";
  std::cerr << "\t     output_max_grid_size=<s> where <s> is the output max_grid_size. If all BoxArrays are the same, this is ignored. [DEF->32]\n";
exit(1);
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
     // ParmParse
     // ---------------------------------------------------------------------
     ParmParse pp;

     if (pp.contains("help")) {
       print_usage(argc,argv);
     }

     // Arbitrary number of input files will be combined and averaged
     int nf = pp.countval("infiles");
     AMREX_ALWAYS_ASSERT(nf>0);
     Vector<std::string> plotFileNames; pp.getarr("infiles",plotFileNames,0,nf);

     // into a single output file
     std::string outfile("plt_averaged");
     pp.query("outfile",outfile);

     // Vairables to keep - will keep all if not specified
     int nvar = pp.countval("variables");
     Vector<std::string> variableNames;
     pp.queryarr("variables",variableNames,0,nvar);
     bool all_vars = nvar > 0 ? false : true;
     Vector<Vector<int>> var_idxs;

     // Maximum number of levels to keep
     int output_max_level = 1000;
     pp.query("output_max_level", output_max_level);
     output_max_level +=1; // account for base level

     // Max grid size in output data
     int output_max_grid_size = 32;
     pp.query("output_max_grid_size", output_max_grid_size);


     // ---------------------------------------------------------------------
     // Execute
     // ---------------------------------------------------------------------

     // First load the metadata of each input plt file
     Print() << "Loading plt file metadata..." << std::endl;
     Vector<pele::physics::pltfilemanager::PltFileManager*> plt_file_data(nf);
     int nlevels = 0;
     for (int i = 0; i < plt_file_data.size(); ++i) {
       plt_file_data[i] = new pele::physics::pltfilemanager::PltFileManager(plotFileNames[i]);
       nlevels = max(nlevels, plt_file_data[i]->getNlev());
       // Verify we have the right vairables
       if (all_vars) {
         if (i==0) {
           variableNames = plt_file_data[i]->getVariableList();
           nvar = variableNames.size();
         } else {
           Vector<std::string> variableNamesTest = plt_file_data[i]->getVariableList();
           if (variableNamesTest.size() != nvar) {
             amrex::Abort("All plt files must have same number of variables unless variable list is specified. File: " + plotFileNames[i]);
           }
           for (int var = 0; var < nvar; ++var) {
             if (variableNames[var] != variableNamesTest[var]) {
               amrex::Abort("All plt files must have same variables unless variable list is specified. File: " + plotFileNames[i]);
             }
           }
         }
       } else {
         Vector<int> var_idx_loc;
         Vector<std::string> variableNamesPlt = plt_file_data[i]->getVariableList();
         for (int var = 0; var < nvar; ++var) {
           int pvar;
           for (pvar = 0; pvar < variableNamesPlt.size(); ++pvar) {
             if (variableNames[var] == variableNamesPlt[pvar]) {
               var_idx_loc.push_back(pvar);
               break;
             }
           }
           if (pvar == variableNamesPlt.size()) {
             amrex::Abort("Variable '" + variableNames[var] + "' not found in file: " + plotFileNames[i]);
           }
         }
         var_idxs.push_back(var_idx_loc);
       }
     }
     nlevels = min(nlevels, output_max_level);
     Print() << " -> Combining " << nf << " files across " << nlevels << " levels" << std::endl;

     // Loop over each input file to get union of boxes on each level
     // On any level with all same BoxArray, we use that without modification
     Print() << "Finding the combined grids..." << std::endl;
     Vector<BoxArray> combined_boxes;
     Vector<Geometry> level_geometries;
     Vector<int> boxarray_all_same(nlevels, 1);
     for (int i = 0; i < plt_file_data.size(); ++i) {
       int nlevels_file = min(nlevels, plt_file_data[i]->getNlev());
       for (int lev = 0; lev < nlevels_file; ++lev) {
         // Verify we have the same geometry
         if (level_geometries.size() <= lev) {
           level_geometries.push_back(plt_file_data[i]->getGeom(lev));
         } else {
           bool same_domain = AlmostEqual(plt_file_data[i]->getGeom(lev).ProbDomain(), level_geometries[lev].ProbDomain());
           bool same_box = plt_file_data[i]->getGeom(lev).Domain() == level_geometries[lev].Domain();
           if (!(same_domain and same_box)) {
             amrex::Abort("All plt files must have the same geometry");
           }
         }
         // Now combine the boxes - if not the same
         if (combined_boxes.size() <= lev) {
           combined_boxes.push_back(BoxArray(plt_file_data[i]->getGrid(lev)));
         } else {
           BoxList boxlist_file{plt_file_data[i]->getGrid(lev)};
           BoxList boxlist_combined{combined_boxes[lev]};
           if (boxlist_combined != boxlist_file) {
             boxlist_combined.catenate(boxlist_file);
             combined_boxes[lev] = BoxArray(boxlist_combined);
             combined_boxes[lev].removeOverlap();
             boxarray_all_same[lev] = 0;
           }
         }
       }
     }

     // Create the data structures to read in the data and keep running sums
     Vector<MultiFab> running_data(nlevels);
     Vector<MultiFab> tmp_data(nlevels);
     Vector<IntVect> refRatios(nlevels-1);
     for (int lev = 0; lev < nlevels; ++lev) {
       if (!boxarray_all_same[lev]) {
         combined_boxes[lev].maxSize(output_max_grid_size);
       }
       DistributionMapping dmap = DistributionMapping(combined_boxes[lev]);
       tmp_data[lev].define(combined_boxes[lev], dmap, nvar, 0);
       running_data[lev].define(combined_boxes[lev], dmap, nvar, 0);
       running_data[lev].setVal(0.0);
       if (lev > 0) {
         int rr = int(level_geometries[lev-1].CellSize(0) / level_geometries[lev].CellSize(0));
         refRatios[lev-1] = {AMREX_D_DECL(rr,rr,rr)};
       }
     }

     // Fillpatch tmp_data from each pltfile and add to running data
     Print() << "Fillpatching and combining..." << std::endl;
     for (int i = 0; i < plt_file_data.size(); ++i) {
       for (int lev = 0; lev < nlevels; ++lev) {
         if (all_vars) {
           plt_file_data[i]->fillPatchFromPlt(lev, level_geometries[lev], 0, 0, nvar, tmp_data[lev]);
         } else {
           for (int var = 0; var < nvar; ++var) {
             plt_file_data[i]->fillPatchFromPlt(lev, level_geometries[lev], var_idxs[i][var], var, 1, tmp_data[lev]);
           }
         }
         MultiFab::Add(running_data[lev], tmp_data[lev], 0, 0, nvar, 0);
       }
       delete plt_file_data[i];
     }

     // Divide by number of files to get average
     Real factor = 1.0 / Real(nf);
     for (int lev = 0; lev < nlevels; ++lev) {
       running_data[lev].mult(factor);
     }

     // Save the final plt file
     Print() << "Saving final plt file..." << std::endl;
     Vector<int> stepidx(nlevels,0);
     WriteMultiLevelPlotfile(outfile,nlevels, GetVecOfConstPtrs(running_data),variableNames,level_geometries,0.0,stepidx,refRatios);
     Print() << "Done." << std::endl;
   }
   amrex::Finalize();
   return 0;
}
