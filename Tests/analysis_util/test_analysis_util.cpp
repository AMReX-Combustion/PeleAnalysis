#include <sstream>
#include <string>

#include <AMReX.H>
#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_Print.H>
#include <AMReX_RealBox.H>
#include <AMReX_Vector.H>

#include <analysis_util.H>

using namespace amrex;

// ---------------------------------------------------------------------------
// Pass/fail tracking
// ---------------------------------------------------------------------------
static int s_pass = 0;
static int s_fail = 0;

static void do_pass(const std::string& label)
{
    amrex::Print() << "  PASS  " << label << "\n";
    ++s_pass;
}
static void do_fail(const std::string& label, const std::string& detail = "")
{
    amrex::Print() << "  FAIL  " << label;
    if (!detail.empty()) amrex::Print() << "  [" << detail << "]";
    amrex::Print() << "\n";
    ++s_fail;
}

#define ASSERT_TRUE(label, expr) \
    do { if (expr) do_pass(label); else do_fail(label, #expr " was false"); } while(0)

#define ASSERT_EQ(label, a, b) \
    do { if ((a) == (b)) do_pass(label); \
         else do_fail(label, "got " + std::to_string(a) + ", expected " + std::to_string(b)); } while(0)

#define ASSERT_STR_EQ(label, a, b) \
    do { if ((a) == (b)) do_pass(label); \
         else do_fail(label, "got '" + (a) + "', expected '" + (b) + "'"); } while(0)

#define ASSERT_NEAR(label, a, b, tol) \
    do { if (std::abs(static_cast<double>(a) - static_cast<double>(b)) <= (tol)) \
             do_pass(label); \
         else do_fail(label, "got " + std::to_string(static_cast<double>(a)) + \
                      ", expected " + std::to_string(static_cast<double>(b))); } while(0)

// ---------------------------------------------------------------------------
// Group 1: String utilities
// ---------------------------------------------------------------------------
static void group_string_utils()
{
    amrex::Print() << "\n=== Group 1: String Utilities ===\n";

    // T-STR-1
    ASSERT_STR_EQ("T-STR-1 get_file_root path",
                  analysis_util::get_file_root("path/to/plt00001"), "plt00001");

    // T-STR-2
    ASSERT_STR_EQ("T-STR-2 get_file_root bare filename",
                  analysis_util::get_file_root("plt00001"), "plt00001");

    // T-STR-3
    ASSERT_STR_EQ("T-STR-3 get_file_root root-relative",
                  analysis_util::get_file_root("/plt00001"), "plt00001");

    amrex::Vector<std::string> vars = {"a", "b", "c"};

    // T-STR-4
    ASSERT_EQ("T-STR-4 find_var_index middle",
              analysis_util::find_var_index(vars, "b", false), 1);

    // T-STR-5
    ASSERT_EQ("T-STR-5 find_var_index first element",
              analysis_util::find_var_index(vars, "a", false), 0);

    // T-STR-6
    ASSERT_EQ("T-STR-6 find_var_index last element",
              analysis_util::find_var_index(vars, "c", false), 2);

    // T-STR-7
    ASSERT_EQ("T-STR-7 find_var_index not found returns -1",
              analysis_util::find_var_index(vars, "d", false), -1);

    // T-STR-8: exact match — "var" must not match "var_extra"
    amrex::Vector<std::string> vars2 = {"var", "var_extra"};
    ASSERT_EQ("T-STR-8 find_var_index exact match",
              analysis_util::find_var_index(vars2, "var", false), 0);

    // T-STR-9: parse_title and parse_var_names space-separated
    {
        std::istringstream ss("My Title\nvar1 var2 var3\n");
        std::string title = analysis_util::parse_title(ss);
        ASSERT_STR_EQ("T-STR-9 parse_title", title, "My Title");
        auto names = analysis_util::parse_var_names(ss);
        ASSERT_EQ("T-STR-9 parse_var_names count", static_cast<int>(names.size()), 3);
        if (static_cast<int>(names.size()) == 3) {
            ASSERT_STR_EQ("T-STR-9 parse_var_names[0]", names[0], "var1");
            ASSERT_STR_EQ("T-STR-9 parse_var_names[1]", names[1], "var2");
            ASSERT_STR_EQ("T-STR-9 parse_var_names[2]", names[2], "var3");
        }
    }

    // T-STR-10: parse_var_names comma+space delimited
    {
        std::istringstream ss("var1, var2, var3\n");
        auto names = analysis_util::parse_var_names(ss);
        ASSERT_EQ("T-STR-10 parse_var_names comma count", static_cast<int>(names.size()), 3);
        if (static_cast<int>(names.size()) >= 1) {
            ASSERT_STR_EQ("T-STR-10 parse_var_names comma[0]", names[0], "var1");
        }
    }
}

// ---------------------------------------------------------------------------
// Group 2: AMReX plotfile round-trip
// ---------------------------------------------------------------------------
static void group_plotfile_io(const std::string& run_dir)
{
    amrex::Print() << "\n=== Group 2: Plotfile I/O ===\n";

    const IntVect lo(AMREX_D_DECL(0, 0, 0));
    const IntVect hi(AMREX_D_DECL(7, 7, 7));  // 8^3 cells
    Box domain(lo, hi);
    BoxArray ba(domain);
    DistributionMapping dm(ba);
    RealBox rb({AMREX_D_DECL(0.0, 0.0, 0.0)}, {AMREX_D_DECL(1.0, 1.0, 1.0)});
    const int coord = 0;
    const int is_per[AMREX_SPACEDIM] = {0};
    Geometry geom(domain, &rb, coord, is_per);
    Vector<Geometry> geoms = {geom};

    // T-PLT-1: single variable, constant value, round-trip
    {
        const Real val = 3.14;
        Vector<MultiFab> mf(1);
        mf[0].define(ba, dm, 1, 0);
        mf[0].setVal(val);

        const std::string plt = run_dir + "/plt_roundtrip_1var";
        analysis_util::write_plotfile(plt, mf, {"myvar"}, geoms);

        auto data = analysis_util::read_plotfile(plt, {"myvar"});
        ASSERT_EQ("T-PLT-1 n_lev", data.n_lev, 1);
        ASSERT_NEAR("T-PLT-1 value round-trip",
                    data.mf[0].min(0), val, 1e-12);
        ASSERT_NEAR("T-PLT-1 value max",
                    data.mf[0].max(0), val, 1e-12);
    }

    // T-PLT-2: two-variable file, load only one
    {
        Vector<MultiFab> mf(1);
        mf[0].define(ba, dm, 2, 0);
        mf[0].setVal(1.0, 0, 1);
        mf[0].setVal(2.0, 1, 1);

        const std::string plt = run_dir + "/plt_roundtrip_2var";
        analysis_util::write_plotfile(plt, mf, {"comp0", "comp1"}, geoms);

        auto data = analysis_util::read_plotfile(plt, {"comp1"});
        ASSERT_EQ("T-PLT-2 loaded ncomp", data.mf[0].nComp(), 1);
        ASSERT_NEAR("T-PLT-2 loaded value", data.mf[0].min(0), 2.0, 1e-12);
    }

    // T-PLT-3: PlotfileData.var_names matches the requested variables (not all file vars)
    //          so that write_plotfile(outfile, data) has consistent mf / var_names sizes.
    {
        Vector<MultiFab> mf(1);
        mf[0].define(ba, dm, 2, 0);
        mf[0].setVal(0.0);
        const std::string plt = run_dir + "/plt_allvars";
        analysis_util::write_plotfile(plt, mf, {"alpha", "beta"}, geoms);
        auto data = analysis_util::read_plotfile(plt, {"alpha"});
        ASSERT_EQ("T-PLT-3 var_names size", static_cast<int>(data.var_names.size()), 1);
        if (static_cast<int>(data.var_names.size()) == 1) {
            ASSERT_STR_EQ("T-PLT-3 var_names[0]", data.var_names[0], "alpha");
        }
    }

    // T-PLT-4: finest_level cap
    {
        Vector<MultiFab> mf(1);
        mf[0].define(ba, dm, 1, 0);
        mf[0].setVal(0.0);
        const std::string plt = run_dir + "/plt_levelcap";
        analysis_util::write_plotfile(plt, mf, {"f"}, geoms);
        auto data = analysis_util::read_plotfile(plt, {"f"}, /*finest_level=*/0);
        ASSERT_EQ("T-PLT-4 n_lev after cap", data.n_lev, 1);
    }

    // T-PLT-5: n_grow ghost cells
    {
        Vector<MultiFab> mf(1);
        mf[0].define(ba, dm, 1, 0);
        mf[0].setVal(0.0);
        const std::string plt = run_dir + "/plt_ngrow";
        analysis_util::write_plotfile(plt, mf, {"f"}, geoms);
        auto data = analysis_util::read_plotfile(plt, {"f"}, 1000, /*n_grow=*/1);
        ASSERT_TRUE("T-PLT-5 ghost cells",
                    data.mf[0].nGrowVect() == IntVect(1));
    }

    // T-PLT-6: is_per propagated to geometry
    {
        Vector<MultiFab> mf(1);
        mf[0].define(ba, dm, 1, 0);
        mf[0].setVal(0.0);
        const std::string plt = run_dir + "/plt_isper";
        analysis_util::write_plotfile(plt, mf, {"f"}, geoms);

        Vector<int> per = {1, 0, 0};
        auto data = analysis_util::read_plotfile(plt, {"f"}, 1000, 0, per);
        ASSERT_TRUE("T-PLT-6 isPeriodic(0)", data.geoms[0].isPeriodic(0));
        ASSERT_TRUE("T-PLT-6 not isPeriodic(1)", !data.geoms[0].isPeriodic(1));
    }

    // T-PLT-7: write_plotfile with default ref_ratios, single level
    {
        Vector<MultiFab> mf(1);
        mf[0].define(ba, dm, 1, 0);
        mf[0].setVal(0.0);
        const std::string plt = run_dir + "/plt_defratios";
        analysis_util::write_plotfile(plt, mf, {"f"}, geoms, 0.0, {});
        auto data = analysis_util::read_plotfile(plt, {"f"});
        ASSERT_EQ("T-PLT-7 single-level write succeeds", data.n_lev, 1);
    }

    // T-PLT-8: Header nComp matches variable count (checked via shell: sed in run_tests.sh)
    {
        Vector<MultiFab> mf(1);
        mf[0].define(ba, dm, 2, 0);
        mf[0].setVal(0.0);
        const std::string plt = run_dir + "/plt_ncomp_check";
        analysis_util::write_plotfile(plt, mf, {"v0", "v1"}, geoms);
        do_pass("T-PLT-8 write two-var plotfile for header check");
    }

    // T-PLT-GUARD: write_plotfile redistributes when #boxes < #MPI ranks.
    // Verifies that the utility-level MPI guard in write_plotfile prevents
    // collective deadlocks on single-box domains with many ranks.
    {
        const Real val = 7.77;
        // Intentionally single-box — the guard inside write_plotfile must
        // redistribute before calling WriteMultiLevelPlotfile.
        Box domain2(IntVect(AMREX_D_DECL(0,0,0)), IntVect(AMREX_D_DECL(7,7,7)));
        BoxArray ba2(domain2);               // 1 box
        DistributionMapping dm2(ba2);
        Geometry geom2(domain2, &rb, coord, is_per);
        Vector<MultiFab> mf2(1);
        mf2[0].define(ba2, dm2, 1, 0);
        mf2[0].setVal(val);
        const std::string plt2 = run_dir + "/plt_guard_1box";
        analysis_util::write_plotfile(plt2, mf2, {"q"}, {geom2});
        auto data2 = analysis_util::read_plotfile(plt2, {"q"});
        ASSERT_NEAR("T-PLT-GUARD value preserved", data2.mf[0].min(0), val, 1e-12);
        const int nprocs = ParallelDescriptor::NProcs();
        if (nprocs > 1) {
            ASSERT_TRUE("T-PLT-GUARD boxes redistributed (≥ nprocs)",
                        static_cast<int>(data2.mf[0].boxArray().size()) >= nprocs);
        }
    }
}

// ---------------------------------------------------------------------------
// Group 3: MEF file round-trip
// ---------------------------------------------------------------------------
static void group_mef_io(const std::string& run_dir)
{
    amrex::Print() << "\n=== Group 3: MEF File I/O ===\n";

    // Build a simple MEFData: 4 nodes (x,y,z,field), 2 triangles
    const int n_nodes = 4;
    const int n_comp  = 4;
    const int n_elts  = 2;
    const int npe     = 3;

    // 1D box along x: (0,0,0) to (n_nodes-1, 0, 0)
    IntVect lo(AMREX_D_DECL(0, 0, 0));
    IntVect hi(AMREX_D_DECL(n_nodes - 1, 0, 0));
    Box nodeBox(lo, hi);

    auto make_data = [&](const std::string& title) {
        analysis_util::MEFData d;
        d.title         = title;
        d.var_names     = {"x", "y", "z", "field"};
        d.n_elts        = n_elts;
        d.nodes_per_elt = npe;
        d.nodes.resize(nodeBox, n_comp);
        for (int comp = 0; comp < n_comp; ++comp)
            for (int i = 0; i < n_nodes; ++i)
                d.nodes(IntVect(AMREX_D_DECL(i, 0, 0)), comp) =
                    static_cast<Real>(comp * 10 + i);
        d.connectivity.resize(n_elts * npe);
        d.connectivity = {1, 2, 3, 1, 3, 4};
        return d;
    };

    // T-MEF-1/2/3/4/5: basic round-trip
    {
        auto data_in = make_data("Test MEF");
        const std::string mef_path = run_dir + "/test_basic.mef";
        analysis_util::write_mef(mef_path, data_in);
        auto data_out = analysis_util::read_mef(mef_path);

        ASSERT_STR_EQ("T-MEF-1 title", data_out.title, data_in.title);
        ASSERT_EQ("T-MEF-2 var_names size",
                  static_cast<int>(data_out.var_names.size()),
                  static_cast<int>(data_in.var_names.size()));
        if (data_out.var_names.size() == data_in.var_names.size()) {
            for (int i = 0; i < static_cast<int>(data_in.var_names.size()); ++i)
                ASSERT_STR_EQ("T-MEF-2 var_names[" + std::to_string(i) + "]",
                               data_out.var_names[i], data_in.var_names[i]);
        }
        ASSERT_EQ("T-MEF-3 n_elts",        data_out.n_elts,        data_in.n_elts);
        ASSERT_EQ("T-MEF-3 nodes_per_elt",  data_out.nodes_per_elt, data_in.nodes_per_elt);

        // T-MEF-4: node values
        bool nodes_ok = true;
        for (int comp = 0; comp < n_comp && nodes_ok; ++comp)
            for (int i = 0; i < n_nodes && nodes_ok; ++i) {
                Real expected = static_cast<Real>(comp * 10 + i);
                Real got = data_out.nodes(IntVect(AMREX_D_DECL(i, 0, 0)), comp);
                if (std::abs(got - expected) > 1e-12) nodes_ok = false;
            }
        ASSERT_TRUE("T-MEF-4 node values", nodes_ok);

        // T-MEF-5: connectivity
        bool conn_ok = data_out.connectivity.size() == data_in.connectivity.size();
        for (int i = 0; i < static_cast<int>(data_in.connectivity.size()) && conn_ok; ++i)
            conn_ok = (data_out.connectivity[i] == data_in.connectivity[i]);
        ASSERT_TRUE("T-MEF-5 connectivity", conn_ok);
    }

    // T-MEF-6: single triangle
    {
        analysis_util::MEFData d;
        d.title         = "Single triangle";
        d.var_names     = {"x", "y", "z"};
        d.n_elts        = 1;
        d.nodes_per_elt = 3;
        IntVect h(AMREX_D_DECL(2, 0, 0));
        d.nodes.resize(Box(IntVect(0), h), 3);
        d.nodes.setVal(0.0);
        d.connectivity = {0, 1, 2};

        const std::string p = run_dir + "/test_single_tri.mef";
        analysis_util::write_mef(p, d);
        auto out = analysis_util::read_mef(p);
        ASSERT_EQ("T-MEF-6 single tri n_elts", out.n_elts, 1);
        ASSERT_EQ("T-MEF-6 single tri conn size",
                  static_cast<int>(out.connectivity.size()), 3);
    }

    // T-MEF-7: multi-component (nComp > 3)
    {
        auto d = make_data("Multi-comp");
        d.var_names = {"x", "y", "z", "T", "rho", "u"};
        IntVect h(AMREX_D_DECL(n_nodes - 1, 0, 0));
        d.nodes.resize(Box(IntVect(0), h), 6);
        d.nodes.setVal(7.0);
        const std::string p = run_dir + "/test_multicomp.mef";
        analysis_util::write_mef(p, d);
        auto out = analysis_util::read_mef(p);
        ASSERT_EQ("T-MEF-7 nComp preserved", out.nodes.nComp(), 6);
        ASSERT_NEAR("T-MEF-7 multi-comp value", out.nodes.min(0), 7.0, 1e-12);
    }

    // T-MEF-8: title with spaces
    {
        auto d = make_data("My surface title with spaces");
        const std::string p = run_dir + "/test_title_spaces.mef";
        analysis_util::write_mef(p, d);
        auto out = analysis_util::read_mef(p);
        ASSERT_STR_EQ("T-MEF-8 title with spaces",
                      out.title, "My surface title with spaces");
    }
}

// ---------------------------------------------------------------------------
// Group 4: get_covered_mf
// ---------------------------------------------------------------------------
static void group_get_covered_mf()
{
    amrex::Print() << "\n=== Group 4: get_covered_mf ===\n";

    // T-COV-1: single level — all cells uncovered
    {
        IntVect lo(AMREX_D_DECL(0, 0, 0)), hi(AMREX_D_DECL(7, 7, 7));
        Box domain(lo, hi);
        BoxArray ba(domain);
        DistributionMapping dm(ba);
        RealBox rb({AMREX_D_DECL(0.0, 0.0, 0.0)}, {AMREX_D_DECL(1.0, 1.0, 1.0)});
        const int is_per[AMREX_SPACEDIM] = {0};
        Geometry geom(domain, &rb, 0, is_per);

        Vector<MultiFab> mf(1);
        mf[0].define(ba, dm, 1, 0);
        Vector<int> no_ratios;
        auto mask = analysis_util::get_covered_mf(mf, no_ratios);
        const long total_cells = AMREX_D_TERM(8, * 8, * 8);
        ASSERT_EQ("T-COV-1 all uncovered", mask[0].sum(0),
                  static_cast<long>(total_cells));
    }

    // T-COV-2: two-level — fine-covered coarse cells = 0
    {
        // Coarse: 8^3, Fine: 4^3 covering coarse cells [2..5]^3, ref=2
        Box coarse_domain(IntVect(AMREX_D_DECL(0,0,0)), IntVect(AMREX_D_DECL(7,7,7)));
        Box fine_domain  (IntVect(AMREX_D_DECL(4,4,4)), IntVect(AMREX_D_DECL(11,11,11)));

        BoxArray ba0(coarse_domain), ba1(fine_domain);
        DistributionMapping dm0(ba0), dm1(ba1);

        RealBox rb({AMREX_D_DECL(0.0,0.0,0.0)}, {AMREX_D_DECL(1.0,1.0,1.0)});
        const int ip[AMREX_SPACEDIM] = {0};
        Geometry geom0(coarse_domain, &rb, 0, ip);
        Geometry geom1(fine_domain,   &rb, 0, ip);

        Vector<MultiFab> mfs(2);
        mfs[0].define(ba0, dm0, 1, 0);
        mfs[1].define(ba1, dm1, 1, 0);
        Vector<int> ref_ratios_2lev = {2};
        auto mask = analysis_util::get_covered_mf(mfs, ref_ratios_2lev);

        // Fine covers coarse cells [2..5]^3 (fine [4..11] / 2 = coarse [2..5])
        const long covered = AMREX_D_TERM(4, * 4, * 4);  // 4^3 coarse cells covered
        const long total_coarse = AMREX_D_TERM(8, * 8, * 8);
        ASSERT_EQ("T-COV-2 uncovered coarse cells",
                  mask[0].sum(0), static_cast<long>(total_coarse - covered));
        const long total_fine = AMREX_D_TERM(8, * 8, * 8);
        ASSERT_EQ("T-COV-2 all fine uncovered", mask[1].sum(0),
                  static_cast<long>(total_fine));
    }

    // T-COV-4: integrate() still calls get_covered_mf() (link-time + runtime check)
    {
        IntVect lo(AMREX_D_DECL(0,0,0)), hi(AMREX_D_DECL(7,7,7));
        Box domain(lo, hi);
        BoxArray ba(domain);
        DistributionMapping dm(ba);
        RealBox rb({AMREX_D_DECL(0.0,0.0,0.0)}, {AMREX_D_DECL(1.0,1.0,1.0)});
        const int ip[AMREX_SPACEDIM] = {0};
        Geometry geom(domain, &rb, 0, ip);

        Vector<MultiFab> mf(1);
        mf[0].define(ba, dm, 1, 0);
        mf[0].setVal(1.0);
        Vector<int> all_axes = {AMREX_D_DECL(0, 1, 2)};
        Vector<int> no_ratios;
        auto result = analysis_util::integrate(mf, {geom}, no_ratios, 0, 1, all_axes);
        ASSERT_TRUE("T-COV-4 integrate still works", result != nullptr);
    }
}

// ---------------------------------------------------------------------------
// Group 5: integrate regression
// ---------------------------------------------------------------------------
static void group_integrate()
{
    amrex::Print() << "\n=== Group 5: integrate Regression ===\n";

    IntVect lo(AMREX_D_DECL(0,0,0)), hi(AMREX_D_DECL(7,7,7));
    Box domain(lo, hi);
    BoxArray ba(domain);
    DistributionMapping dm(ba);
    RealBox rb({AMREX_D_DECL(0.0,0.0,0.0)}, {AMREX_D_DECL(1.0,1.0,1.0)});
    const int ip[AMREX_SPACEDIM] = {0};
    Geometry geom(domain, &rb, 0, ip);

    const Real val = 3.0;
    Vector<MultiFab> mf(1);
    mf[0].define(ba, dm, 1, 0);
    mf[0].setVal(val);
    Vector<int> no_ratios;

    // T-INT-1: all-axis integration of a constant field
    // Expected: val * domain_volume = 3.0 * 1.0 = 3.0
    {
        Vector<int> all_axes = {AMREX_D_DECL(0, 1, 2)};
        auto result = analysis_util::integrate(mf, {geom}, no_ratios, 0, 1, all_axes);
        AMREX_ALWAYS_ASSERT(result != nullptr);
        ASSERT_NEAR("T-INT-1 all-axis integral", (*result)[0], val * 1.0, 1e-10);
    }

    // T-INT-2: x-axis integration only — result has ny*nz entries
    {
        Vector<int> x_axis = {0};
        auto result = analysis_util::integrate(mf, {geom}, no_ratios, 0, 1, x_axis);
        AMREX_ALWAYS_ASSERT(result != nullptr);
        const int ny = AMREX_D_PICK(1, 8, 8);
        const int nz = AMREX_D_PICK(1, 1, 8);
        const int expected_size = ny * nz;
        ASSERT_EQ("T-INT-2 x-integral result size",
                  static_cast<int>(result->size()), expected_size);
        // integrate() accumulates field * cell_volume (3D), so x-integration
        // gives sum_i val*dV = nx * val * dx^3 = 8 * 3.0 * (1/8)^3 = 3.0/64
        ASSERT_NEAR("T-INT-2 x-integral value", (*result)[0], val / 64.0, 1e-10);
    }
}

// ---------------------------------------------------------------------------
// Abort-path entry point (for subprocess tests in the shell script)
// ---------------------------------------------------------------------------
static void run_abort_test(const std::string& test_name)
{
    amrex::Print() << "Running abort test: " << test_name << "\n";
    if (test_name == "find_var_index") {
        amrex::Vector<std::string> vars = {"a", "b"};
        analysis_util::find_var_index(vars, "c", /*abort=*/true);
    } else if (test_name == "read_plotfile_bad_var") {
        IntVect lo(AMREX_D_DECL(0,0,0)), hi(AMREX_D_DECL(3,3,3));
        Box domain(lo, hi);
        BoxArray ba(domain);
        DistributionMapping dm(ba);
        RealBox rb({AMREX_D_DECL(0.0,0.0,0.0)}, {AMREX_D_DECL(1.0,1.0,1.0)});
        const int ip[AMREX_SPACEDIM] = {0};
        Geometry geom(domain, &rb, 0, ip);
        Vector<MultiFab> mf(1);
        mf[0].define(ba, dm, 1, 0);
        mf[0].setVal(0.0);
        analysis_util::write_plotfile("/tmp/plt_abort_var", mf, {"f"}, {geom});
        analysis_util::read_plotfile("/tmp/plt_abort_var", {"nonexistent"});
    } else if (test_name == "read_plotfile_no_file") {
        analysis_util::read_plotfile("/does/not/exist/plt00000", {"f"});
    } else if (test_name == "read_mef_no_file") {
        analysis_util::read_mef("/does/not/exist/data.mef");
    } else {
        amrex::Abort("Unknown abort test: " + test_name);
    }
}

// ---------------------------------------------------------------------------
// Main
// ---------------------------------------------------------------------------
int main(int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    {
        ParmParse pp;

        // Abort-path mode: run one abort test and exit
        std::string abort_test;
        pp.query("abort_test", abort_test);
        if (!abort_test.empty()) {
            run_abort_test(abort_test);
            // Should not reach here if the abort fired
            amrex::Abort("abort_test did not abort as expected");
        }

        // Normal mode
        std::string run_dir = ".";
        pp.query("run_dir", run_dir);

        group_string_utils();
        group_plotfile_io(run_dir);
        group_mef_io(run_dir);
        group_get_covered_mf();
        group_integrate();

        amrex::Print() << "\n========================================\n";
        amrex::Print() << "  PASSED: " << s_pass << "\n";
        amrex::Print() << "  FAILED: " << s_fail << "\n";
        amrex::Print() << "========================================\n";
    }
    amrex::Finalize();
    return (s_fail > 0) ? 1 : 0;
}
