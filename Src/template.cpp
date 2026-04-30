#include <string>
#include <iostream>

#include <AMReX_AmrData.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>

#include <analysis_util.H>

using namespace amrex;

static void
print_usage(int, char* argv[])
{
  std::cerr
    << "Usage:\n"
    << "  " << argv[0] << " infile=FILE vars=\"VAR1 VAR2 ...\" [OPTIONS]\n\n"

    << "Required arguments:\n"
    << "  infile=FILE        AMReX plotfile to read\n"
    << "  vars=\"...\"         Space-separated list of variable names to "
       "process\n\n"

    << "Options:\n"
    << "  finestLevel=N      Cap AMR levels (default: all levels in file)\n"
    << "  is_per=\"1 1 1\"     Periodicity flags per dimension (default: all "
       "1)\n"
    << "  verbose            Enable verbose plotfile I/O\n"
    << "  -h, --help         Show this help message\n\n"

    << "Visit PeleAnalysis/Src/InputSamples for examples or refer to "
    << "the documentation.\n";

  std::exit(1);
}

int
main(int argc, char* argv[])
{
  Initialize(argc, argv);
  {
    if (
      argc < 2 || std::strcmp(argv[1], "-h") == 0 ||
      std::strcmp(argv[1], "--help") == 0) {
      print_usage(argc, argv);
    }
    ParmParse pp;

    if (pp.contains("verbose"))
      AmrData::SetVerbose(true);

    std::string infile;
    pp.get("infile", infile);

    Vector<std::string> vars;
    const int nvars = pp.countval("vars");
    if (nvars == 0) {
      amrex::Abort("template: 'vars' is required — specify the variable names "
                   "to process.\n"
                   "  Example: vars=\"density velocity_x\"");
    }
    pp.getarr("vars", vars, 0, nvars);

    int p_finestLevel = 1000;
    pp.query("finestLevel", p_finestLevel);

    Vector<int> is_per(AMREX_SPACEDIM, 1);
    pp.queryarr("is_per", is_per, 0, AMREX_SPACEDIM);

    // -------------------------------------------------------------------------
    // Read plotfile
    // -------------------------------------------------------------------------
    auto data = analysis_util::read_plotfile(
      infile, vars, p_finestLevel,
      /*n_grow=*/0, is_per);

    Print() << "Read " << data.n_lev << " level(s) from " << infile << "\n";

    // -------------------------------------------------------------------------
    // Process data (modify data.mf in-place)
    // -------------------------------------------------------------------------
    const int nComp = static_cast<int>(vars.size());
    for (int lev = 0; lev < data.n_lev; ++lev) {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (MFIter mfi(data.mf[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();
        auto const& a = data.mf[lev].array(mfi);
        amrex::ParallelFor(
          bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            for (int n = 0; n < nComp; ++n) {
              a(i, j, k, n) = a(i, j, k, n); // replace with actual computation
            }
          });
      }
    }

    // -------------------------------------------------------------------------
    // Write output
    // -------------------------------------------------------------------------
    const std::string outfile = analysis_util::get_file_root(infile) + "_temp";
    Print() << "Writing new data to " << outfile << "\n";
    analysis_util::write_plotfile(
      outfile, data.mf, vars, data.geoms, data.time);
  }
  Finalize();
  return 0;
}
