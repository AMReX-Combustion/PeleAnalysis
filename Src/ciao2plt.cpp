#include <map>
#include <string>
#include <iostream>

#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_VisMF.H>
#include <AMReX_Print.H>

#include "H5Cpp.h"

using namespace amrex;
using namespace H5;

static void
print_usage(int, char* argv[])
{
  std::cerr
    << "usage:\n"
    << "  " << argv[0] << " infile=<file.h5> [OPTIONS]\n\n"
    << "Options:\n"
    << "  outfile=<name>          output plt name [default: plt_<stem>]\n"
    << "  vars=\"v1 v2 ...\"        fields to read [default: all]\n"
    << "  max_grid_size=<N>       AMReX box size [default: 32]\n"
    << "  per=\"0 0 0\"             periodicity (fallback if not in file)\n";
  std::exit(1);
}

// Strip directory and extension: /path/to/data.h5 -> data
static std::string
stem(const std::string& path)
{
  auto sep = path.rfind('/');
  std::string base = (sep == std::string::npos) ? path : path.substr(sep + 1);
  auto dot = base.rfind('.');
  return (dot == std::string::npos) ? base : base.substr(0, dot);
}

// H5Lexists emits diagnostics in HDF5 >=1.14 when an intermediate path
// component is absent. Suppress the error stack for the duration of the check.
static bool
h5_exists(hid_t loc, const char* name)
{
  H5E_auto_t old_func;
  void* old_data;
  H5Eget_auto(H5E_DEFAULT, &old_func, &old_data);
  H5Eset_auto(H5E_DEFAULT, nullptr, nullptr);
  htri_t ret = H5Lexists(loc, name, H5P_DEFAULT);
  H5Eset_auto(H5E_DEFAULT, old_func, old_data);
  return ret > 0;
}

// Find the single mesh group (skip IO-information).
static std::string
discover_mesh_group(const H5File& file)
{
  Group root = file.openGroup("/");
  std::string mesh;
  hsize_t n = root.getNumObjs();
  for (hsize_t i = 0; i < n; ++i) {
    std::string name = root.getObjnameByIdx(i);
    if (root.getObjTypeByIdx(i) == H5G_GROUP && name != "IO-information") {
      if (!mesh.empty())
        Abort("Multiple mesh groups found in HDF5 file; only one is supported");
      mesh = name;
    }
  }
  if (mesh.empty())
    Abort("No mesh group found in HDF5 file");
  Print() << "Mesh group: " << mesh << "\n";
  return mesh;
}

// List all field names: cv_data_real datasets + scalars/SC attribute names.
static Vector<std::string>
list_fields(const H5File& file, const std::string& mesh)
{
  Vector<std::string> names;

  Group cvg = file.openGroup("/" + mesh + "/data/cv_data_real");
  hsize_t n = cvg.getNumObjs();
  for (hsize_t i = 0; i < n; ++i)
    if (cvg.getObjTypeByIdx(i) == H5G_DATASET)
      names.push_back(cvg.getObjnameByIdx(i));

  std::string scpath = "/" + mesh + "/data/scalars";
  if (h5_exists(file.getId(), scpath.c_str())) {
    Group scalg = file.openGroup(scpath);
    if (h5_exists(scalg.getId(), "SC")) {
      DataSet sc = scalg.openDataSet("SC");
      hsize_t dims[4] = {};
      sc.getSpace().getSimpleExtentDims(dims);
      for (int i = 1; i <= (int)dims[0]; ++i) {
        Attribute attr = sc.openAttribute("Index " + std::to_string(i));
        std::string val;
        attr.read(attr.getStrType(), val);
        names.push_back(val);
      }
    }
  }
  return names;
}

// Read simulation time from globals_r0/time (required).
static Real
read_time(const H5File& file, const std::string& mesh)
{
  std::string path = "/" + mesh + "/data/globals_r0/time";
  if (!h5_exists(file.getId(), path.c_str()))
    Abort("Time dataset not found at " + path);
  double t = 0.0;
  file.openDataSet(path).read(&t, PredType::NATIVE_DOUBLE);
  Print() << "Time = " << t << "\n";
  return (Real)t;
}

// Read periodicity from sd_info if present; fall back to per= input parameter.
static Vector<int>
read_periodicity(const H5File& file, const std::string& mesh)
{
  Vector<int> per(3, 0);
  std::string sdpath = "/" + mesh + "/geometry/sd_info";
  if (h5_exists(file.getId(), sdpath.c_str())) {
    Group sd = file.openGroup(sdpath);
    if (h5_exists(sd.getId(), "xper")) {
      auto ri = [&](const std::string& name) {
        int v = 0;
        sd.openDataSet(name).read(&v, PredType::NATIVE_INT);
        return v;
      };
      per[0] = ri("xper");
      per[1] = ri("yper");
      per[2] = ri("zper");
      Print() << "Periodicity (from file): " << per[0] << " " << per[1] << " "
              << per[2] << "\n";
      return per;
    }
  }
  ParmParse pp;
  pp.queryarr("per", per, 0, 3);
  return per;
}

// Read coordinate system type from sd_info/icyl if present; fall back to
// coord_sys= input.
static int
read_coord_sys(const H5File& file, const std::string& mesh)
{
  std::string sdpath = "/" + mesh + "/geometry/sd_info";
  if (h5_exists(file.getId(), sdpath.c_str())) {
    Group sd = file.openGroup(sdpath);
    if (h5_exists(sd.getId(), "icyl")) {
      int v = 0;
      sd.openDataSet("icyl").read(&v, PredType::NATIVE_INT);
      Print() << "Coord sys (from file): " << v << "\n";
      return v;
    }
  }
  int coord = 0;
  ParmParse pp;
  pp.query("coord_sys", coord);
  return coord;
}

// Build BoxArray and Geometry from grid node coordinates.
// Grid coords are node-based: N+1 nodes define N cells.
static void
setup_grid(
  const H5File& file,
  const std::string& mesh,
  const Vector<int>& per,
  int coord_sys,
  int max_grid_size,
  BoxArray& ba,
  Geometry& geom)
{
  std::string gp = "/" + mesh + "/geometry/grid/";

  auto read_coord = [&](const std::string& name) {
    DataSet ds = file.openDataSet(gp + name);
    hsize_t dim[1] = {};
    ds.getSpace().getSimpleExtentDims(dim);
    Vector<Real> v(dim[0]);
    ds.read(v.data(), PredType::NATIVE_DOUBLE);
    return v;
  };

  Vector<Real> x = read_coord("x");
  Vector<Real> y = read_coord("y");
  Vector<Real> z = read_coord("z");

  int nx = (int)x.size() - 1;
  int ny = (int)y.size() - 1;
  int nz = (int)z.size() - 1;

  Print() << "Grid: " << nx << " x " << ny << " x " << nz << " cells\n";

#if AMREX_SPACEDIM == 2
  if (nz != 1)
    Abort("DIM=2 but file has nz=" + std::to_string(nz) + " (expected 1)");
#endif

  Box domain(IntVect::Zero, IntVect(AMREX_D_DECL(nx - 1, ny - 1, nz - 1)));
  Real plo[AMREX_SPACEDIM] = {AMREX_D_DECL(x[0], y[0], z[0])};
  Real phi[AMREX_SPACEDIM] = {AMREX_D_DECL(x[nx], y[ny], z[nz])};
  RealBox rb(plo, phi);
#if AMREX_SPACEDIM != 2
  if (coord_sys != 0)
    Abort("coord_sys != 0 (cylindrical/RZ) requires a DIM=2 build");
#endif
  geom = Geometry(domain, &rb, coord_sys, per.data());

  ba = BoxArray(domain);
  ba.maxSize(max_grid_size);
}

// Read all requested fields into data.
// Fields come from cv_data_real (3D datasets) or scalars/SC (4D stacked array).
static void
read_fields(
  const H5File& file,
  const std::string& mesh,
  const Vector<std::string>& vars,
  MultiFab& data)
{
  std::string cvbase = "/" + mesh + "/data/cv_data_real/";
  std::string scpath = "/" + mesh + "/data/scalars/SC";
  bool has_sc = h5_exists(file.getId(), scpath.c_str());

  // Build scalar name -> 0-based index map.
  std::map<std::string, int> sc_map;
  if (has_sc) {
    DataSet sc = file.openDataSet(scpath);
    hsize_t dims[4] = {};
    sc.getSpace().getSimpleExtentDims(dims);
    for (int i = 1; i <= (int)dims[0]; ++i) {
      Attribute attr = sc.openAttribute("Index " + std::to_string(i));
      std::string name;
      attr.read(attr.getStrType(), name);
      sc_map[name] = i - 1;
    }
  }

  // Validate all vars exist before starting I/O.
  for (const auto& v : vars) {
    bool in_cv = h5_exists(file.getId(), (cvbase + v).c_str());
    bool in_sc = sc_map.count(v) > 0;
    if (!in_cv && !in_sc) {
      auto avail = list_fields(file, mesh);
      std::string msg = "Field '" + v + "' not found. Available:";
      for (const auto& a : avail)
        msg += " " + a;
      Abort(msg);
    }
  }

  for (int n = 0; n < (int)vars.size(); ++n) {
    const std::string& vname = vars[n];
    Print() << "Reading " << vname << "\n";

    bool is_cv = h5_exists(file.getId(), (cvbase + vname).c_str());
    DataSet ds = file.openDataSet(is_cv ? cvbase + vname : scpath);

    for (MFIter mfi(data); mfi.isValid(); ++mfi) {
      const Box& box = mfi.validbox();
      const int* lo = box.loVect();
      const int* hi = box.hiVect();

      int nx_b = hi[0] - lo[0] + 1;
      int ny_b = hi[1] - lo[1] + 1;
#if AMREX_SPACEDIM == 3
      int nz_b = hi[2] - lo[2] + 1;
      int lo2 = lo[2];
#else
      int nz_b = 1;
      int lo2 = 0;
#endif
      Vector<Real> buf(nx_b * ny_b * nz_b);
      DataSpace fsp = ds.getSpace();

      if (is_cv) {
        hsize_t cnt[3] = {(hsize_t)nz_b, (hsize_t)ny_b, (hsize_t)nx_b};
        hsize_t off[3] = {(hsize_t)lo2, (hsize_t)lo[1], (hsize_t)lo[0]};
        fsp.selectHyperslab(H5S_SELECT_SET, cnt, off);
        DataSpace msp(3, cnt);
        ds.read(buf.data(), PredType::NATIVE_DOUBLE, msp, fsp);
      } else {
        int si = sc_map.at(vname);
        hsize_t cnt[4] = {1, (hsize_t)nz_b, (hsize_t)ny_b, (hsize_t)nx_b};
        hsize_t off[4] = {
          (hsize_t)si, (hsize_t)lo2, (hsize_t)lo[1], (hsize_t)lo[0]};
        fsp.selectHyperslab(H5S_SELECT_SET, cnt, off);
        hsize_t mcnt[3] = {(hsize_t)nz_b, (hsize_t)ny_b, (hsize_t)nx_b};
        DataSpace msp(3, mcnt);
        ds.read(buf.data(), PredType::NATIVE_DOUBLE, msp, fsp);
      }

      Array4<Real> const& arr = data.array(mfi);
      int lo0 = lo[0], lo1 = lo[1], nxb = nx_b, nyb = ny_b;
      amrex::ParallelFor(
        box, [arr, &buf, n, lo0, lo1, lo2, nxb,
              nyb] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
          arr(i, j, k, n) =
            buf[(k - lo2) * nyb * nxb + (j - lo1) * nxb + (i - lo0)];
        });
    }
  }
}

int
main(int argc, char* argv[])
{
  amrex::Initialize(argc, argv);
  {
    const Real strt_time = amrex::second();

    if (argc < 2)
      print_usage(argc, argv);

    ParmParse pp;
    if (pp.contains("help"))
      print_usage(argc, argv);

    std::string infile;
    pp.get("infile", infile);

    std::string outfile = "plt_" + stem(infile);
    pp.query("outfile", outfile);

    int max_grid_size = 32;
    pp.query("max_grid_size", max_grid_size);

#ifdef AMREX_USE_MPI
    hid_t fapl_id = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(
      fapl_id, ParallelDescriptor::Communicator(), MPI_INFO_NULL);
    H5File file(
      infile.c_str(), H5F_ACC_RDONLY, H5::FileCreatPropList::DEFAULT,
      H5::FileAccPropList(fapl_id));
    H5Pclose(fapl_id);
#else
    H5File file(infile.c_str(), H5F_ACC_RDONLY);
#endif
    std::string mesh = discover_mesh_group(file);

    Vector<std::string> vars;
    if (pp.countval("vars") > 0) {
      pp.getarr("vars", vars);
    } else {
      vars = list_fields(file, mesh);
      Print() << vars.size() << " fields:";
      for (const auto& v : vars)
        Print() << " " << v;
      Print() << "\n";
    }

    Real time = read_time(file, mesh);
    Vector<int> per = read_periodicity(file, mesh);
    int coord_sys = read_coord_sys(file, mesh);

    BoxArray ba;
    Geometry geom;
    setup_grid(file, mesh, per, coord_sys, max_grid_size, ba, geom);

    MultiFab data(ba, DistributionMapping(ba), (int)vars.size(), 0);
    read_fields(file, mesh, vars, data);

    // Cap the number of plotfile data files via the n_files option (AMReX)
    int n_files = amrex::VisMF::GetNOutFiles();
    pp.query("n_files", n_files);
    amrex::VisMF::SetNOutFiles(n_files);

    Print() << "Writing " << outfile << "...\n";
    WriteSingleLevelPlotfile(outfile, data, vars, geom, time, 0);

    Real elapsed = amrex::second() - strt_time;
    ParallelDescriptor::ReduceRealMax(elapsed);
    Print() << "Done. Run time: " << elapsed << " s\n";
  }
  amrex::Finalize();
  return 0;
}
