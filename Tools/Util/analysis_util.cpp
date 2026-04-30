#include <analysis_util.H>

using namespace amrex;

namespace analysis_util {

std::string
get_file_root(const std::string& infile)
{
  std::vector<std::string> tokens = Tokenize(infile, std::string("/"));
  return tokens[tokens.size() - 1];
}

int
find_var_index(
  const Vector<std::string>& var_names,
  const std::string&         name,
  bool                       abort_if_not_found)
{
  for (int i = 0; i < static_cast<int>(var_names.size()); ++i) {
    if (var_names[i] == name)
      return i;
  }
  if (abort_if_not_found) {
    amrex::Abort(
      "analysis_util::find_var_index: variable '" + name +
      "' not found in plotfile");
  }
  return -1;
}

amrex::Vector<amrex::iMultiFab>
get_covered_mf(
  const amrex::Vector<amrex::MultiFab>& mf,
  const amrex::Vector<int>&             ref_ratios)
{
  const int nlev = static_cast<int>(mf.size());
  amrex::Vector<amrex::iMultiFab> mask_mf;
  mask_mf.reserve(nlev);
  for (int lev = 0; lev < nlev; ++lev) {
    mask_mf.emplace_back(mf[lev].boxArray(), mf[lev].DistributionMap(), 1, 0);
    mask_mf.back().setVal(1);
  }
  const int finest_level = nlev - 1;
  for (int lev = 0; lev < finest_level; ++lev) {
    amrex::BoxArray baf = mf[lev + 1].boxArray();
    baf.coarsen(ref_ratios[lev]);
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    {
      std::vector<std::pair<int, amrex::Box>> isects;
      for (amrex::MFIter mfi(mask_mf[lev], amrex::TilingIfNotGPU());
           mfi.isValid(); ++mfi) {
        auto const& mask = mask_mf[lev].array(mfi);
        baf.intersections(mf[lev].boxArray()[mfi.index()], isects);
        for (const auto& is : isects) {
          amrex::ParallelFor(
            is.second, [mask] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
              mask(i, j, k) = 0;
            });
        }
      }
    }
  }
  return mask_mf;
}

} // namespace analysis_util
