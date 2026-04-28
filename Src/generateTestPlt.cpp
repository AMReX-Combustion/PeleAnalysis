#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_TagBox.H>
#include <AMReX_Cluster.H>
#include <algorithm>
#include <string>
#include <map>

using namespace amrex;

static void
print_usage(int, char* argv[])
{
  std::cerr
    << "Usage:\n"
    << "  " << argv[0] << " [OPTIONS]\n\n"

    << "Required arguments:\n"
    << "  amr.n_cell=N [N N]           Number of cells per dimension (level "
       "0)\n"
    << "  geometry.prob_lo=X [Y Z]     Physical lower bound of the domain\n"
    << "  geometry.prob_hi=X [Y Z]     Physical upper bound of the domain\n"
    << "  geometry.is_periodic=N [N N] Periodicity flags per dimension\n"
    << "  field.names=NAME [NAME ...]  Space-separated list of field names\n"
    << "  <name>.type=TYPE             Field type for each named field\n"
    << "  <name>.output_name=S         Variable name written to plotfile "
       "header\n"
    << "                               (default: same as field name; use to "
       "embed\n"
    << "                                '/' or other ParmParse-reserved "
       "characters)\n\n"

    << "Geometry options:\n"
    << "  geometry.coord_sys=N         0 = Cartesian (default), 1 = "
       "cylindrical/RZ\n\n"

    << "Field types (for <name>.type):\n"
    << "  constant                     Uniform value everywhere\n"
    << "  plane_step / plane_smooth    Step/smooth transition at a coordinate "
       "plane\n"
    << "                               axis=0|1|2  position=F  "
       "[smooth_width=F]\n"
    << "  double_plane_step/smooth     Slab between two parallel planes\n"
    << "                               axis=0|1|2  position=F  "
       "position_second=F\n"
    << "  sphere_step / sphere_smooth  Step/smooth sphere (3D) or circle (2D)\n"
    << "                               center=X [Y Z]  radius=F  "
       "[smooth_width=F]\n"
    << "  ring_step / ring_smooth      Spherical shell (annulus in 2D)\n"
    << "                               center=X [Y Z]  radius_inner=F  "
       "radius_outer=F\n"
    << "  cylinder_step / cylinder_smooth  Step/smooth infinite cylinder\n"
    << "                               axis=0|1|2 (default 2)  center=X [Y Z]\n"
    << "                               radius=F  [smooth_width=F]\n"
    << "                               Radial distance is measured from the "
       "named\n"
    << "                               axis line, not from a point.\n"
    << "  sine                         Separable sin*cos*sin wave\n"
    << "                               frequency=Fx [Fy Fz]  phase=Px [Py Pz]\n"
    << "                               amplitude=A  offset=B\n\n"

    << "AMR options (multilevel):\n"
    << "  amr.max_level=N              Maximum refinement level (default 0)\n"
    << "  amr.ref_ratio=N [N ...]      Refinement ratio per level (default 2)\n"
    << "  amr.grid_eff=F               Clustering efficiency threshold "
       "(default 0.7)\n"
    << "  amr.blocking_factor=N        Minimum box size (default 8)\n"
    << "  amr.n_error_buf=N [N ...]    Buffer cells around tagged region "
       "(default 1)\n"
    << "  amr.max_grid_size=N          Maximum box size (default 32)\n"
    << "  amr.refinement_indicators=NAME [NAME ...]\n"
    << "                               Space-separated list of refinement "
       "criteria\n"
    << "  amr.<name>.in_box_lo=X [Y Z] Physical lower bound of refinement "
       "region\n"
    << "  amr.<name>.in_box_hi=X [Y Z] Physical upper bound of refinement "
       "region\n"
    << "  amr.<name>.value_greater=F   Refine where field > F\n"
    << "  amr.<name>.value_less=F      Refine where field < F\n"
    << "  amr.<name>.adjacent_difference_greater=F  Refine where max neighbor "
       "diff > F\n"
    << "  amr.<name>.field_name=NAME   Field to evaluate (value-based "
       "criteria)\n"
    << "  amr.<name>.max_level=N       Max level for this indicator (default: "
       "amr.max_level)\n\n"

    << "  -h, --help                   Show this help message\n\n"

    << "Visit PeleAnalysis/Src/InputSamples for examples or refer to "
    << "the documentation.\n";

  std::exit(1);
}

// Enum for field types
enum class FieldType {
  Constant,
  PlaneStep,
  PlaneSmooth,
  DoublePlaneStep,
  DoublePlaneSmooth,
  CircleStep,
  CircleSmooth,
  RingStep,
  RingSmooth,
  CylinderStep,
  CylinderSmooth,
  Sine
};

// Structure to hold field configuration
struct FieldConfig
{
  char name[16];
  FieldType type;

  // Common parameters
  Real value_inside = 1.0;
  Real value_outside = 0.0;

  // Plane parameters
  int plane_axis = 0; // 0=x, 1=y, 2=z
  Real plane_position = 0.5;
  Real plane_position_second = 0.6;
  Real smooth_width = 0.1;

  // Circle/Sphere parameters
  GpuArray<Real, AMREX_SPACEDIM> center = {AMREX_D_DECL(0.5, 0.5, 0.5)};
  Real radius = 0.25;

  // Ring/Spherical shell parameters
  Real radius_inner = 0.25;
  Real radius_outer = 0.25;

  // Sine parameters
  GpuArray<Real, AMREX_SPACEDIM> frequency = {AMREX_D_DECL(1.0, 1.0, 1.0)};
  GpuArray<Real, AMREX_SPACEDIM> phase = {AMREX_D_DECL(0.0, 0.0, 0.0)};
  Real amplitude = 1.0;
  Real offset = 0.0;
};

// Enum for refinement criterion types
enum class RefinementType {
  InBox,
  ValueGreater,
  ValueLess,
  AdjacentDiffGreater
};

// Structure to hold refinement indicator configuration
struct RefinementRegion
{
  std::string name;
  int max_level = -1; // -1 resolved to global max_level after parsing
  RefinementType type = RefinementType::InBox;

  // InBox parameters
  GpuArray<Real, AMREX_SPACEDIM> lo = {AMREX_D_DECL(0.0, 0.0, 0.0)};
  GpuArray<Real, AMREX_SPACEDIM> hi = {AMREX_D_DECL(1.0, 1.0, 1.0)};

  // Value-based parameters
  std::string field_name;
  int field_comp = -1; // resolved after field name parsing
  Real threshold = 0.0;
};

// Helper function to parse field type
FieldType
parseFieldType(const std::string& type_str)
{
  static const std::map<std::string, FieldType> type_map = {
    {"constant", FieldType::Constant},
    {"plane_step", FieldType::PlaneStep},
    {"plane_smooth", FieldType::PlaneSmooth},
    {"double_plane_step", FieldType::DoublePlaneStep},
    {"double_plane_smooth", FieldType::DoublePlaneSmooth},
    {"circle_step", FieldType::CircleStep},
    {"circle_smooth", FieldType::CircleSmooth},
    {"sphere_step", FieldType::CircleStep},
    {"sphere_smooth", FieldType::CircleSmooth},
    {"ring_step", FieldType::RingStep},
    {"ring_smooth", FieldType::RingSmooth},
    {"spherical_shell_step", FieldType::RingStep},
    {"spherical_shell_smooth", FieldType::RingSmooth},
    {"cylinder_step", FieldType::CylinderStep},
    {"cylinder_smooth", FieldType::CylinderSmooth},
    {"sine", FieldType::Sine}};

  auto it = type_map.find(type_str);
  if (it != type_map.end()) {
    return it->second;
  }

  amrex::Abort("Unknown field type: " + type_str);
  return FieldType::Constant;
}

// Smooth step function
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Real
smoothstep(Real edge0, Real edge1, Real x)
{
  Real t = std::max(0.0, std::min(1.0, (x - edge0) / (edge1 - edge0)));
  return t * t * t * (t * (6.0 * t - 15.0) + 10.0);
}

// Function to evaluate field value at a point
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Real
evaluateField(
  const FieldConfig& config,
  Real x,
  Real y,
  Real z,
  const GpuArray<Real, AMREX_SPACEDIM>& dx,
  const GpuArray<Real, AMREX_SPACEDIM>& prob_lo)
{

  Real result = 0.0;

  switch (config.type) {
  case FieldType::Constant:
    result = config.value_inside;
    break;

  case FieldType::PlaneStep: {
    Real coord;
    if (config.plane_axis == 0)
      coord = x;
    else if (config.plane_axis == 1)
      coord = y;
    else
      coord = z;

    result = (coord < config.plane_position) ? config.value_inside
                                             : config.value_outside;
    break;
  }

  case FieldType::PlaneSmooth: {
    Real coord;
    if (config.plane_axis == 0)
      coord = x;
    else if (config.plane_axis == 1)
      coord = y;
    else
      coord = z;

    Real edge0 = config.plane_position - config.smooth_width / 2.0;
    Real edge1 = config.plane_position + config.smooth_width / 2.0;
    Real blend = smoothstep(edge0, edge1, coord);
    result = config.value_inside * (1.0 - blend) + config.value_outside * blend;
    break;
  }

  case FieldType::DoublePlaneStep: {
    Real coord;
    if (config.plane_axis == 0)
      coord = x;
    else if (config.plane_axis == 1)
      coord = y;
    else
      coord = z;

    if (
      (coord < config.plane_position) ||
      (coord > config.plane_position_second)) {
      result = config.value_outside;
    } else {
      result = config.value_inside;
    }
    break;
  }

  case FieldType::DoublePlaneSmooth: {
    Real coord;
    if (config.plane_axis == 0)
      coord = x;
    else if (config.plane_axis == 1)
      coord = y;
    else
      coord = z;

    // First plane transition
    Real edge0_first = config.plane_position - config.smooth_width / 2.0;
    Real edge1_first = config.plane_position + config.smooth_width / 2.0;
    Real blend_first = smoothstep(edge0_first, edge1_first, coord);

    // Second plane transition
    Real edge0_second =
      config.plane_position_second - config.smooth_width / 2.0;
    Real edge1_second =
      config.plane_position_second + config.smooth_width / 2.0;
    Real blend_second = smoothstep(edge0_second, edge1_second, coord);

    // Combine: outside if coord < first OR coord > second, inside if first <
    // coord < second Assumes plane_position < plane_position_second
    // blend_first: 0 when coord < first, 1 when coord > first
    // blend_second: 0 when coord < second, 1 when coord > second
    Real blend = blend_first * (1.0 - blend_second);

    result = config.value_inside * blend + config.value_outside * (1.0 - blend);
    break;
  }

  case FieldType::CircleStep: {
#if AMREX_SPACEDIM == 2
    Real dx_val = x - config.center[0];
    Real dy_val = y - config.center[1];
    Real r = std::sqrt(dx_val * dx_val + dy_val * dy_val);
#else
    Real dx_val = x - config.center[0];
    Real dy_val = y - config.center[1];
    Real dz_val = z - config.center[2];
    Real r = std::sqrt(dx_val * dx_val + dy_val * dy_val + dz_val * dz_val);
#endif
    result = (r < config.radius) ? config.value_inside : config.value_outside;
    break;
  }

  case FieldType::CircleSmooth: {
#if AMREX_SPACEDIM == 2
    Real dx_val = x - config.center[0];
    Real dy_val = y - config.center[1];
    Real r = std::sqrt(dx_val * dx_val + dy_val * dy_val);
#else
    Real dx_val = x - config.center[0];
    Real dy_val = y - config.center[1];
    Real dz_val = z - config.center[2];
    Real r = std::sqrt(dx_val * dx_val + dy_val * dy_val + dz_val * dz_val);
#endif
    Real edge0 = config.radius - config.smooth_width / 2.0;
    Real edge1 = config.radius + config.smooth_width / 2.0;
    Real blend = smoothstep(edge0, edge1, r);
    result = config.value_inside * (1.0 - blend) + config.value_outside * blend;
    break;
  }

  case FieldType::RingStep: {
#if AMREX_SPACEDIM == 2
    Real dx_val = x - config.center[0];
    Real dy_val = y - config.center[1];
    Real r = std::sqrt(dx_val * dx_val + dy_val * dy_val);
#else
    Real dx_val = x - config.center[0];
    Real dy_val = y - config.center[1];
    Real dz_val = z - config.center[2];
    Real r = std::sqrt(dx_val * dx_val + dy_val * dy_val + dz_val * dz_val);
#endif

    if ((r < config.radius_inner) || (r > config.radius_outer)) {
      result = config.value_outside;
    } else {
      result = config.value_inside;
    }
    break;
  }

  case FieldType::RingSmooth: {
#if AMREX_SPACEDIM == 2
    Real dx_val = x - config.center[0];
    Real dy_val = y - config.center[1];
    Real r = std::sqrt(dx_val * dx_val + dy_val * dy_val);
#else
    Real dx_val = x - config.center[0];
    Real dy_val = y - config.center[1];
    Real dz_val = z - config.center[2];
    Real r = std::sqrt(dx_val * dx_val + dy_val * dy_val + dz_val * dz_val);
#endif
    // Smooth transition at inner radius (outside -> inside)
    Real edge0_inner = config.radius_inner - config.smooth_width / 2.0;
    Real edge1_inner = config.radius_inner + config.smooth_width / 2.0;
    Real blend_inner = smoothstep(edge0_inner, edge1_inner, r);

    // Smooth transition at outer radius (inside -> outside)
    Real edge0_outer = config.radius_outer - config.smooth_width / 2.0;
    Real edge1_outer = config.radius_outer + config.smooth_width / 2.0;
    Real blend_outer = smoothstep(edge0_outer, edge1_outer, r);

    // Combine: outside if r < inner OR r > outer, inside if inner < r < outer
    // blend_inner: 0 when r < inner, 1 when r > inner
    // blend_outer: 0 when r < outer, 1 when r > outer
    Real blend = blend_inner * (1.0 - blend_outer);

    result = config.value_inside * blend + config.value_outside * (1.0 - blend);
    break;
  }

  case FieldType::CylinderStep: {
    Real r;
    if (config.plane_axis == 0) { // cylinder axis along x; transverse = y,z
#if AMREX_SPACEDIM == 3
      Real dy = y - config.center[1];
      Real dz = z - config.center[2];
      r = std::sqrt(dy * dy + dz * dz);
#else
      r = std::abs(y - config.center[1]); // 2D: only one transverse direction
#endif
    } else if (config.plane_axis == 1) { // cylinder axis along y; transverse =
                                         // x,z
#if AMREX_SPACEDIM == 3
      Real dx = x - config.center[0];
      Real dz = z - config.center[2];
      r = std::sqrt(dx * dx + dz * dz);
#else
      r = std::abs(x - config.center[0]); // 2D: only one transverse direction
#endif
    } else { // cylinder axis along z; transverse = x,y
      Real dx = x - config.center[0];
      Real dy = y - config.center[1];
      r = std::sqrt(dx * dx + dy * dy);
    }
    result = (r < config.radius) ? config.value_inside : config.value_outside;
    break;
  }

  case FieldType::CylinderSmooth: {
    Real r;
    if (config.plane_axis == 0) {
#if AMREX_SPACEDIM == 3
      Real dy = y - config.center[1];
      Real dz = z - config.center[2];
      r = std::sqrt(dy * dy + dz * dz);
#else
      r = std::abs(y - config.center[1]);
#endif
    } else if (config.plane_axis == 1) {
#if AMREX_SPACEDIM == 3
      Real dx = x - config.center[0];
      Real dz = z - config.center[2];
      r = std::sqrt(dx * dx + dz * dz);
#else
      r = std::abs(x - config.center[0]);
#endif
    } else {
      Real dx = x - config.center[0];
      Real dy = y - config.center[1];
      r = std::sqrt(dx * dx + dy * dy);
    }
    Real edge0 = config.radius - config.smooth_width / 2.0;
    Real edge1 = config.radius + config.smooth_width / 2.0;
    Real blend = smoothstep(edge0, edge1, r);
    result = config.value_inside * (1.0 - blend) + config.value_outside * blend;
    break;
  }

  case FieldType::Sine: {
    Real arg_x = 2.0 * M_PI * config.frequency[0] * x + config.phase[0];
    Real arg_y = 2.0 * M_PI * config.frequency[1] * y + config.phase[1];

#if AMREX_SPACEDIM == 3
    Real arg_z = 2.0 * M_PI * config.frequency[2] * z + config.phase[2];
#else
    Real arg_z = 0.0;
#endif
    result = config.offset + config.amplitude * std::sin(arg_x) *
                               std::cos(arg_y) * std::sin(arg_z);
    break;
  }
  }

  return result;
}

// Build the BoxArray for refinement level `lev` (lev >= 1).
// Tags coarse cells analytically, buffers, clusters, and refines to fine index
// space.
BoxArray
buildFineBoxArray(
  int lev,
  const Geometry& coarse_geom,
  const BoxArray& coarse_ba,
  const DistributionMapping& coarse_dm,
  const Vector<RefinementRegion>& regions,
  const FieldConfig* h_configs, // host pointer to field configs
  int ref_ratio_lev,
  int n_error_buf,
  Real grid_eff,
  int /*blocking_factor*/,
  int max_grid_size)
{
  const auto coarse_dx = coarse_geom.CellSizeArray();
  const auto prob_lo = coarse_geom.ProbLoArray();
  const Box coarse_domain = coarse_geom.Domain();

  TagBoxArray tba(coarse_ba, coarse_dm, n_error_buf);
  tba.setVal(TagBox::CLEAR);

  for (const auto& reg : regions) {
    if (reg.max_level < lev) {
      continue;
    }

    if (reg.type == RefinementType::InBox) {
      const GpuArray<Real, AMREX_SPACEDIM> rlo = reg.lo;
      const GpuArray<Real, AMREX_SPACEDIM> rhi = reg.hi;

      for (MFIter mfi(tba); mfi.isValid(); ++mfi) {
        const Box& box = mfi.validbox();
        auto tag = tba[mfi].array();
        amrex::Loop(box, [=](int i, int j, int k) noexcept {
          Real x = prob_lo[0] + (i + 0.5) * coarse_dx[0];
          Real y = prob_lo[1] + (j + 0.5) * coarse_dx[1];
#if AMREX_SPACEDIM == 3
          Real z = prob_lo[2] + (k + 0.5) * coarse_dx[2];
#else
          Real z = 0.0;
#endif
          if (AMREX_D_TERM(
                x >= rlo[0] && x <= rhi[0], &&y >= rlo[1] && y <= rhi[1],
                &&z >= rlo[2] && z <= rhi[2])) {
            tag(i, j, k) = TagBox::SET;
          }
        });
      }

    } else if (
      reg.type == RefinementType::ValueGreater ||
      reg.type == RefinementType::ValueLess) {
      const int comp = reg.field_comp;
      const Real thresh = reg.threshold;
      const bool do_greater = (reg.type == RefinementType::ValueGreater);
      const FieldConfig cfg = h_configs[comp];

      for (MFIter mfi(tba); mfi.isValid(); ++mfi) {
        const Box& box = mfi.validbox();
        auto tag = tba[mfi].array();
        amrex::Loop(box, [=](int i, int j, int k) noexcept {
          Real x = prob_lo[0] + (i + 0.5) * coarse_dx[0];
          Real y = prob_lo[1] + (j + 0.5) * coarse_dx[1];
#if AMREX_SPACEDIM == 3
          Real z = prob_lo[2] + (k + 0.5) * coarse_dx[2];
#else
          Real z = 0.0;
#endif
          Real val = evaluateField(cfg, x, y, z, coarse_dx, prob_lo);
          if ((do_greater && val > thresh) || (!do_greater && val < thresh)) {
            tag(i, j, k) = TagBox::SET;
          }
        });
      }

    } else { // AdjacentDiffGreater
      const int comp = reg.field_comp;
      const Real thresh = reg.threshold;
      const FieldConfig cfg = h_configs[comp];

      for (MFIter mfi(tba); mfi.isValid(); ++mfi) {
        const Box& box = mfi.validbox();
        auto tag = tba[mfi].array();
        amrex::Loop(box, [=](int i, int j, int k) noexcept {
          Real x = prob_lo[0] + (i + 0.5) * coarse_dx[0];
          Real y = prob_lo[1] + (j + 0.5) * coarse_dx[1];
#if AMREX_SPACEDIM == 3
          Real z = prob_lo[2] + (k + 0.5) * coarse_dx[2];
#else
          Real z = 0.0;
#endif
          Real val = evaluateField(cfg, x, y, z, coarse_dx, prob_lo);
          Real max_diff = 0.0;

          Real vxp =
            evaluateField(cfg, x + coarse_dx[0], y, z, coarse_dx, prob_lo);
          Real vxm =
            evaluateField(cfg, x - coarse_dx[0], y, z, coarse_dx, prob_lo);
          max_diff = std::max(max_diff, std::abs(vxp - val));
          max_diff = std::max(max_diff, std::abs(vxm - val));

#if AMREX_SPACEDIM >= 2
          Real vyp =
            evaluateField(cfg, x, y + coarse_dx[1], z, coarse_dx, prob_lo);
          Real vym =
            evaluateField(cfg, x, y - coarse_dx[1], z, coarse_dx, prob_lo);
          max_diff = std::max(max_diff, std::abs(vyp - val));
          max_diff = std::max(max_diff, std::abs(vym - val));
#endif

#if AMREX_SPACEDIM == 3
          Real vzp =
            evaluateField(cfg, x, y, z + coarse_dx[2], coarse_dx, prob_lo);
          Real vzm =
            evaluateField(cfg, x, y, z - coarse_dx[2], coarse_dx, prob_lo);
          max_diff = std::max(max_diff, std::abs(vzp - val));
          max_diff = std::max(max_diff, std::abs(vzm - val));
#endif

          if (max_diff > thresh) {
            tag(i, j, k) = TagBox::SET;
          }
        });
      }
    }
  }

  // Expand tagged cells by n_error_buf cells at coarse level
  tba.buffer(IntVect(n_error_buf));

  // Gather tagged cell positions to IO proc (collate is gather-to-root, not
  // allgather). Pattern mirrors AMReX AmrMesh::regrid (AMReX_AmrMesh.cpp ~line
  // 736-775).
  Gpu::PinnedVector<IntVect> pts;
  tba.collate(pts);
  // collate uses ReduceLongSum internally: if numtags==0 it calls clear() on
  // ALL ranks, so pts.empty() is a safe all-rank test for "nothing to refine".
  if (pts.empty()) {
    return BoxArray{};
  }

  // Cluster on IO proc only (non-IO procs have a single-element placeholder in
  // pts).
  BoxList fine_bl;
  if (ParallelDescriptor::IOProcessor()) {
    ClusterList clist(pts.data(), static_cast<Long>(pts.size()));
    clist.chop(grid_eff);
    BoxList coarse_bl = clist.boxList();
    coarse_bl.intersect(coarse_domain);
    coarse_bl.simplify();
    for (Box b : coarse_bl) {
      fine_bl.push_back(b.refine(ref_ratio_lev));
    }
    fine_bl.simplify();
  }

  // Broadcast the BoxList from IO proc to all ranks so every rank builds the
  // same BoxArray and DistributionMapping (required for MPI collectives).
  fine_bl.Bcast();

  if (fine_bl.isEmpty()) {
    return BoxArray{};
  }

  BoxArray ba(fine_bl);
  ba.maxSize(max_grid_size);
  return ba;
}

int
main(int argc, char* argv[])
{
  amrex::Initialize(argc, argv);
  {
    if (argc < 2) {
      print_usage(argc, argv);
    } else if (
      (std::strcmp(argv[1], "-h") == 0) ||
      (std::strcmp(argv[1], "--help") == 0)) {
      print_usage(argc, argv);
    }

    ParmParse pp;
    ParmParse ppamr("amr");
    ParmParse ppgeom("geometry");

    // Geometry: number of cells on level 0
    Vector<int> n_cell(AMREX_SPACEDIM);
    ppamr.getarr("n_cell", n_cell, 0, AMREX_SPACEDIM);

    // Domain bounds
    RealBox real_box;
    Vector<Real> pp_prob_x(AMREX_SPACEDIM, 0.0);
    ppgeom.getarr("prob_lo", pp_prob_x, 0, AMREX_SPACEDIM);
    GpuArray<Real, AMREX_SPACEDIM> prob_lo = {
      AMREX_D_DECL(pp_prob_x[0], pp_prob_x[1], pp_prob_x[2])};
    ppgeom.getarr("prob_hi", pp_prob_x, 0, AMREX_SPACEDIM);
    GpuArray<Real, AMREX_SPACEDIM> prob_hi = {
      AMREX_D_DECL(pp_prob_x[0], pp_prob_x[1], pp_prob_x[2])};

    for (int i = 0; i < AMREX_SPACEDIM; i++) {
      real_box.setLo(i, prob_lo[i]);
      real_box.setHi(i, prob_hi[i]);
    }

    // Periodicity
    IntVect pp_is_per;
    ppgeom.getarr("is_periodic", pp_is_per);
    Array<int, AMREX_SPACEDIM> is_per = {
      AMREX_D_DECL(pp_is_per[0], pp_is_per[1], pp_is_per[2])};

    // Coordinate system: 0 = Cartesian, 1 = cylindrical/RZ
    int coord = 0;
    ppgeom.query("coord_sys", coord);
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
      coord == 0 || coord == 1,
      "geometry.coord_sys must be 0 (Cartesian) or 1 (cylindrical/RZ)");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
      coord == 0 || AMREX_SPACEDIM == 2,
      "geometry.coord_sys=1 (cylindrical/RZ) requires a 2D build (DIM=2)");

    // Level-0 geometry
    IntVect domain_lo(AMREX_D_DECL(0, 0, 0));
    IntVect domain_hi(
      AMREX_D_DECL(n_cell[0] - 1, n_cell[1] - 1, n_cell[2] - 1));
    Box domain(domain_lo, domain_hi);
    Geometry geom0(domain, real_box, coord, is_per);

    // --- AMR parameters ---
    int max_level = 0;
    ppamr.query("max_level", max_level);
    int Nlev = max_level + 1;

    Vector<int> ref_ratio(std::max(max_level, 1), 2);
    if (ppamr.countval("ref_ratio") > 0) {
      Vector<int> rr;
      ppamr.getarr("ref_ratio", rr);
      for (int lev = 0; lev < max_level; ++lev) {
        ref_ratio[lev] =
          (lev < static_cast<int>(rr.size())) ? rr[lev] : rr.back();
      }
    }

    Real grid_eff = 0.7;
    ppamr.query("grid_eff", grid_eff);

    int blocking_factor = 8;
    ppamr.query("blocking_factor", blocking_factor);

    Vector<int> n_error_buf(std::max(max_level, 1), 1);
    if (ppamr.countval("n_error_buf") > 0) {
      Vector<int> neb;
      ppamr.getarr("n_error_buf", neb);
      for (int lev = 0; lev < max_level; ++lev) {
        n_error_buf[lev] =
          (lev < static_cast<int>(neb.size())) ? neb[lev] : neb.back();
      }
    }

    // Refinement indicators
    Vector<RefinementRegion> regions;
    {
      int n_ind = ppamr.countval("refinement_indicators");
      if (n_ind > 0) {
        Vector<std::string> ind_names(n_ind);
        ppamr.getarr("refinement_indicators", ind_names);
        regions.resize(n_ind);
        for (int i = 0; i < n_ind; ++i) {
          RefinementRegion& reg = regions[i];
          reg.name = ind_names[i];
          reg.max_level = max_level;

          ParmParse ppi("amr." + ind_names[i]);
          ppi.query("max_level", reg.max_level);

          if (ppi.countval("in_box_lo") > 0) {
            reg.type = RefinementType::InBox;
            Vector<Real> blo(AMREX_SPACEDIM), bhi(AMREX_SPACEDIM);
            ppi.getarr("in_box_lo", blo, 0, AMREX_SPACEDIM);
            ppi.getarr("in_box_hi", bhi, 0, AMREX_SPACEDIM);
            for (int d = 0; d < AMREX_SPACEDIM; ++d) {
              reg.lo[d] = blo[d];
              reg.hi[d] = bhi[d];
            }
          } else if (ppi.countval("value_greater") > 0) {
            reg.type = RefinementType::ValueGreater;
            ppi.get("value_greater", reg.threshold);
            ppi.get("field_name", reg.field_name);
          } else if (ppi.countval("value_less") > 0) {
            reg.type = RefinementType::ValueLess;
            ppi.get("value_less", reg.threshold);
            ppi.get("field_name", reg.field_name);
          } else if (ppi.countval("adjacent_difference_greater") > 0) {
            reg.type = RefinementType::AdjacentDiffGreater;
            ppi.get("adjacent_difference_greater", reg.threshold);
            ppi.get("field_name", reg.field_name);
          } else {
            amrex::Abort(
              "Refinement indicator '" + ind_names[i] +
              "' has no recognized criterion (in_box_lo, value_greater, "
              "value_less, or adjacent_difference_greater)");
          }
        }
      }
    }

    // Grid decomposition parameters
    int max_grid_size = 32;
    ppamr.query("max_grid_size", max_grid_size);
    int ngrow = 0;
    ppamr.query("ngrow", ngrow);

    // Parse field configurations
    int ncomp = pp.countval("field.names");
    Vector<std::string> field_names(ncomp);
    Vector<FieldConfig> field_configs(ncomp);

    if (ncomp > 0) {
      pp.getarr("field.names", field_names);
      Print() << "Fields:" << std::endl;

      for (int i = 0; i < ncomp; ++i) {
        ParmParse ppf(field_names[i]);
        FieldConfig& config = field_configs[i];
        std::strncpy(config.name, field_names[i].c_str(), 15);
        config.name[15] = '\0';

        std::string type_str;
        ppf.get("type", type_str);
        config.type = parseFieldType(type_str);

        Print() << "  Field " << i << ": " << field_names[i]
                << "  (type: " << type_str << ")" << std::endl;

        // Read parameters based on field type
        switch (config.type) {
        case FieldType::Constant: {
          ppf.query("value", config.value_inside);
          break;
        }

        case FieldType::PlaneStep:
        case FieldType::PlaneSmooth: {
          ppf.query("axis", config.plane_axis);
          ppf.query("position", config.plane_position);
          ppf.query("value_left", config.value_inside);
          ppf.query("value_right", config.value_outside);
          if (config.type == FieldType::PlaneSmooth) {
            ppf.query("smooth_width", config.smooth_width);
          }
          break;
        }

        case FieldType::DoublePlaneStep:
        case FieldType::DoublePlaneSmooth: {
          ppf.query("axis", config.plane_axis);
          ppf.query("position", config.plane_position);
          ppf.query("position_second", config.plane_position_second);
          ppf.query("value_inside", config.value_inside);
          ppf.query("value_outside", config.value_outside);
          if (config.type == FieldType::DoublePlaneSmooth) {
            ppf.query("smooth_width", config.smooth_width);
          }
          break;
        }

        case FieldType::CircleStep:
        case FieldType::CircleSmooth: {
          Vector<Real> pp_center(AMREX_SPACEDIM, 0.0);
          ppf.queryarr("center", pp_center);
          for (int dim = 0; dim < static_cast<int>(pp_center.size()); ++dim) {
            config.center[dim] = pp_center[dim];
          }
          ppf.query("radius", config.radius);
          ppf.query("value_inside", config.value_inside);
          ppf.query("value_outside", config.value_outside);
          if (config.type == FieldType::CircleSmooth) {
            ppf.query("smooth_width", config.smooth_width);
          }
          break;
        }

        case FieldType::RingStep:
        case FieldType::RingSmooth: {
          Vector<Real> pp_center(AMREX_SPACEDIM, 0.0);
          ppf.queryarr("center", pp_center);
          for (int dim = 0; dim < static_cast<int>(pp_center.size()); ++dim) {
            config.center[dim] = pp_center[dim];
          }
          ppf.query("radius_inner", config.radius_inner);
          ppf.query("radius_outer", config.radius_outer);
          ppf.query("value_inside", config.value_inside);
          ppf.query("value_outside", config.value_outside);
          if (config.type == FieldType::RingSmooth) {
            ppf.query("smooth_width", config.smooth_width);
          }
          break;
        }

        case FieldType::CylinderStep:
        case FieldType::CylinderSmooth: {
          config.plane_axis = 2; // default: cylinder axis along z
          ppf.query("axis", config.plane_axis);
          Vector<Real> pp_center(AMREX_SPACEDIM, 0.5);
          ppf.queryarr("center", pp_center);
          for (int dim = 0; dim < static_cast<int>(pp_center.size()); ++dim) {
            config.center[dim] = pp_center[dim];
          }
          ppf.query("radius", config.radius);
          ppf.query("value_inside", config.value_inside);
          ppf.query("value_outside", config.value_outside);
          if (config.type == FieldType::CylinderSmooth) {
            ppf.query("smooth_width", config.smooth_width);
          }
          break;
        }

        case FieldType::Sine: {
          Vector<Real> pp_frequency(AMREX_SPACEDIM, 0.0);
          ppf.queryarr("frequency", pp_frequency);
          for (int dim = 0; dim < static_cast<int>(pp_frequency.size());
               ++dim) {
            config.frequency[dim] = pp_frequency[dim];
          }
          Vector<Real> pp_phase(AMREX_SPACEDIM, 0.0);
          ppf.queryarr("phase", pp_phase);
          for (int dim = 0; dim < static_cast<int>(pp_phase.size()); ++dim) {
            config.phase[dim] = pp_phase[dim];
          }
          ppf.query("amplitude", config.amplitude);
          ppf.query("offset", config.offset);
          break;
        }
        }
      }
    } else {
      Abort("No fields specified...");
    }

    // Resolve field_comp for value-based refinement indicators
    for (auto& reg : regions) {
      if (reg.type != RefinementType::InBox) {
        auto it =
          std::find(field_names.begin(), field_names.end(), reg.field_name);
        if (it == field_names.end()) {
          amrex::Abort(
            "Refinement indicator '" + reg.name +
            "' references unknown field '" + reg.field_name + "'");
        }
        reg.field_comp = static_cast<int>(it - field_names.begin());
      }
    }

    // Build level-0 data structures
    Vector<Geometry> geoms(Nlev);
    Vector<BoxArray> grids(Nlev);
    Vector<DistributionMapping> dmaps(Nlev);
    Vector<MultiFab> mfs(Nlev);

    geoms[0] = geom0;
    grids[0] = BoxArray(domain);
    grids[0].maxSize(max_grid_size);
    dmaps[0] = DistributionMapping(grids[0]);
    mfs[0].define(grids[0], dmaps[0], ncomp, ngrow);

    // Lambda to fill one level analytically
    auto fillLevel = [&](int lev) {
      const auto lev_dx = geoms[lev].CellSizeArray();
      const auto lev_prob_lo = geoms[lev].ProbLoArray();

      Gpu::DeviceVector<FieldConfig> d_configs(ncomp);
      Gpu::copyAsync(
        Gpu::hostToDevice, field_configs.begin(), field_configs.end(),
        d_configs.begin());
      Gpu::streamSynchronize();
      FieldConfig* p_configs = d_configs.data();

      for (MFIter mfi(mfs[lev]); mfi.isValid(); ++mfi) {
        const Box& box = mfi.validbox();
        auto const& fab = mfs[lev].array(mfi);

        amrex::ParallelFor(
          box, ncomp, [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) {
            Real x = lev_prob_lo[0] + (i + 0.5) * lev_dx[0];
            Real y = lev_prob_lo[1] + (j + 0.5) * lev_dx[1];
#if AMREX_SPACEDIM == 3
            Real z = lev_prob_lo[2] + (k + 0.5) * lev_dx[2];
#else
            Real z = 0.0;
#endif
            fab(i, j, k, n) =
              evaluateField(p_configs[n], x, y, z, lev_dx, lev_prob_lo);
          });
      }
    };

    // Fill level 0 first (needed before tagging finer levels)
    fillLevel(0);

    // Build and fill fine levels
    for (int lev = 1; lev < Nlev; ++lev) {
      geoms[lev] = amrex::refine(geoms[lev - 1], ref_ratio[lev - 1]);

      BoxArray fine_ba = buildFineBoxArray(
        lev, geoms[lev - 1], grids[lev - 1], dmaps[lev - 1], regions,
        field_configs.data(), ref_ratio[lev - 1], n_error_buf[lev - 1],
        grid_eff, blocking_factor, max_grid_size);

      if (fine_ba.empty()) {
        amrex::Print() << "Warning: no tagged cells at level " << lev
                       << " — truncating to level " << lev - 1 << "\n";
        Nlev = lev;
        geoms.resize(Nlev);
        grids.resize(Nlev);
        dmaps.resize(Nlev);
        mfs.resize(Nlev);
        break;
      }

      grids[lev] = fine_ba;
      dmaps[lev] = DistributionMapping(grids[lev]);
      mfs[lev].define(grids[lev], dmaps[lev], ncomp, ngrow);

      amrex::Print() << "Level " << lev << ": " << grids[lev].size()
                     << " boxes, " << grids[lev].numPts() << " cells\n";
      fillLevel(lev);
    }

    // Variable names for plotfile (output_name overrides the ParmParse prefix)
    Vector<std::string> varnames(ncomp);
    for (int n = 0; n < ncomp; ++n) {
      varnames[n] = field_names[n];
      ParmParse(field_names[n]).query("output_name", varnames[n]);
    }

    // Write multilevel plot file
    std::string plotfile_name = "pltTestFile";
    pp.query("plotfile_name", plotfile_name);

    Vector<int> isteps(Nlev, 0);
    Vector<IntVect> refRatios(Nlev - 1);
    for (int lev = 0; lev < Nlev - 1; ++lev) {
      refRatios[lev] = IntVect(ref_ratio[lev]);
    }
    amrex::WriteMultiLevelPlotfile(
      plotfile_name, Nlev, GetVecOfConstPtrs(mfs), varnames, geoms, 0.0, isteps,
      refRatios);

    amrex::Print() << "Plotfile created: " << plotfile_name << "\n";
  }
  amrex::Finalize();
  return 0;
}
