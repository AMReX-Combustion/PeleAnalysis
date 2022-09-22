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
  std::cerr << "usage:\n";
  std::cerr << argv[0] << " infile=<s> [options] \n\tOptions:\n";
  std::cerr << "\t     infile=<s> where <s> is a pltfile\n";
  std::cerr << "\t     output_file=<s> where <s> is the flatten pltfile\n";
  std::cerr << "\t     output_level=<s> where <s> is the level of infile kept in output_file [DEF->0]\n";
  std::cerr << "\t     output_max_grid_size=<s> where <s> is the flatten max_grid_size [DEF->64]\n";
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
     // ParmParse
     // ---------------------------------------------------------------------
     ParmParse pp;

     if (pp.contains("help")) {
       print_usage(argc,argv);
     }

     int verbose = 0;
     pp.query("verbose",verbose);
     std::string plotFileName;
     pp.get("infile",plotFileName);
     std::string outfile(getFileRoot(plotFileName) + "_flatten");
     pp.query("output_file",outfile);
     int output_level = 0;
     pp.query("output_level",output_level);
     int output_mgs = 64;
     pp.query("output_max_grid_size",output_mgs);

     // ---------------------------------------------------------------------
     // Initiate PelePhysics PltFileManager
     // ---------------------------------------------------------------------
     if (verbose) {
       Print() << " Reading data from pltfile " << plotFileName << "\n";
     }
     pele::physics::pltfilemanager::PltFileManager pltData(plotFileName);
     Vector<std::string> plt_vars = pltData.getVariableList();
     int nvars = plt_vars.size();

     // Check
     AMREX_ALWAYS_ASSERT(output_level < pltData.getNlev());

     // ---------------------------------------------------------------------
     // Set the outgoing geom/grid/dmap/MF
     // ---------------------------------------------------------------------
     if (verbose) {
       Print() << " Setting outgoing container \n";
     }
     auto geom = pltData.getGeom(output_level);
     BoxArray grid{geom.Domain()};
     grid.maxSize(output_mgs);
     const DistributionMapping dmap{grid};
     MultiFab output(grid,dmap,nvars,0);
         
     // ---------------------------------------------------------------------
     // Interpolate data
     // ---------------------------------------------------------------------
     if (verbose) {
       Print() << " Interpolate data \n";
     }
     pltData.fillPatchFromPlt(0, geom, 0, 0, nvars, output);
     
     // ---------------------------------------------------------------------
     // Write standard pltfile
     // ---------------------------------------------------------------------
     if (verbose) {
       Print() << " Writting data \n";
     }
     WriteSingleLevelPlotfile(outfile, output, plt_vars, geom, pltData.getTime(), 0);

   }
   amrex::Finalize();
   return 0;
}
