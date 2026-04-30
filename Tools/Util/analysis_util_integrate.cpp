#include <analysis_util.H>

namespace analysis_util {

// Integrate the MF over selected axes, excluding fine-covered and EB-covered
// cells.  ref_ratios[k] is the refinement ratio between level k and level k+1.
std::unique_ptr<amrex::Gpu::ManagedVector<amrex::Real>>
integrate(
  const amrex::Vector<amrex::MultiFab>& a_mf,
  const amrex::Vector<amrex::Geometry>& geoms,
  const amrex::Vector<int>&             ref_ratios,
  const int                             scomp,
  const int                             ncomp,
  const amrex::Vector<int>&             axes_to_integrate
#ifdef AMREX_USE_EB
  ,
  const amrex::Vector<amrex::MultiFab>& vfrac_mf
#endif
)
{
  const int nlev = static_cast<int>(a_mf.size());
  amrex::Vector<amrex::iMultiFab> mask_mf =
    get_covered_mf(a_mf, ref_ratios);
  const int finest_level = nlev - 1;
  // Get the domain on the finest level
  amrex::Box probDomain = geoms[finest_level].Domain();
  // Domain size on each side
  amrex::IntVect sides = probDomain.length();

  const bool integrate_x =
    std::find(axes_to_integrate.begin(), axes_to_integrate.end(), 0) !=
    axes_to_integrate.end();
  const bool integrate_y =
    std::find(axes_to_integrate.begin(), axes_to_integrate.end(), 1) !=
    axes_to_integrate.end();
  const bool integrate_z =
    std::find(axes_to_integrate.begin(), axes_to_integrate.end(), 2) !=
    axes_to_integrate.end();

  int red_nx = integrate_x ? 1 : sides[0];
  int red_ny = integrate_y ? 1 : sides[1];
  int red_nz = integrate_z ? 1 : sides[2];

  size_t size = ncomp * red_nx * red_ny * red_nz;

  std::unique_ptr<amrex::Gpu::ManagedVector<amrex::Real>> integral =
    std::make_unique<amrex::Gpu::ManagedVector<amrex::Real>>(size, 0.0);

  int ratio = 1;

  for (int lev = finest_level; lev >= 0; --lev) {
    amrex::MultiFab volume(
      a_mf[lev].boxArray(), a_mf[lev].DistributionMap(), 1, 0);
    geoms[lev].GetVolume(volume);
    auto const& ma       = a_mf[lev].arrays();
    auto const& mask_ma  = mask_mf[lev].const_arrays();
    auto const& volume_ma = volume.const_arrays();
#ifdef AMREX_USE_EB
    auto const& vfrac_ma = vfrac_mf[lev]->const_arrays();
#endif

    int ratio_x = integrate_x ? 1 : ratio;
    int ratio_y = integrate_y ? 1 : ratio;
    int ratio_z = integrate_z ? 1 : ratio;

    int ratio_prod = ratio_x * ratio_y * ratio_z;
    amrex::Real* ptr = integral->data();
    amrex::ParallelFor(
      a_mf[lev], amrex::IntVect(0), ncomp,
      [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k, int n) noexcept {
        if (mask_ma[box_no](i, j, k) == 1) {
          amrex::Real val = ma[box_no](i, j, k, scomp + n) *
                            volume_ma[box_no](i, j, k) /
                            static_cast<amrex::Real>(ratio_prod);
#ifdef AMREX_USE_EB
          val *= vfrac_ma[box_no](i, j, k);
#endif
          for (int rx = 0; rx < ratio_x; rx++) {
            for (int ry = 0; ry < ratio_y; ry++) {
              for (int rz = 0; rz < ratio_z; rz++) {
                int ii = integrate_x ? 0 : i * ratio + rx;
                int jj = integrate_y ? 0 : j * ratio + ry;
                int kk = integrate_z ? 0 : k * ratio + rz;

                int idx = ((kk * red_ny + jj) * red_nx + ii) * ncomp + n;

                amrex::Gpu::Atomic::Add(&ptr[idx], val);
              }
            }
          }
        }
      });
    // Update ratio for the next (coarser) level.  ref_ratios[lev-1] is the
    // refinement factor between level lev-1 and level lev.
    if (lev > 0 && !ref_ratios.empty()) {
      ratio *= ref_ratios[lev - 1];
    }
  }
  amrex::ParallelDescriptor::ReduceRealSum(
    integral->data(), static_cast<int>(integral->size()));
  return integral;
}
} // namespace analysis_util
