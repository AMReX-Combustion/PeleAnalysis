#include <analysis_util.H>

#include <AMReX_DataServices.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_VisMF.H>

#include <cmath>

using namespace amrex;

namespace analysis_util {

int
set_plot_nfiles()
{
  int n_files = VisMF::GetNOutFiles();
  ParmParse pp;
  pp.query("n_files", n_files);
  VisMF::SetNOutFiles(n_files);
  return VisMF::GetNOutFiles();
}

// Returns a BoxArray with at least nprocs boxes by applying maxSize so that
// WriteMultiLevelPlotfile and AmrData::FillVar do not deadlock in MPI mode
// when some ranks would otherwise own no FABs.
static BoxArray
ensure_min_boxes(const BoxArray& ba, int nprocs)
{
  if (static_cast<int>(ba.size()) >= nprocs)
    return ba;
  const int ncells_per_rank =
    std::max(1, static_cast<int>(ba.numPts()) / nprocs);
  const int max_size = std::max(
    1, static_cast<int>(
         std::pow(static_cast<double>(ncells_per_rank), 1.0 / AMREX_SPACEDIM)));
  BoxArray ba_split = ba;
  ba_split.maxSize(max_size);
  return ba_split;
}

PlotfileData
read_plotfile(
  const std::string& infile,
  const Vector<std::string>& var_names,
  int finest_level,
  int n_grow,
  const Vector<int>& is_per_in)
{
  DataServices::SetBatchMode();
  Amrvis::FileType fileType(Amrvis::NEWPLT);
  DataServices dataServices(infile, fileType);
  if (!dataServices.AmrDataOk()) {
    amrex::Abort(
      "analysis_util::read_plotfile: cannot open plotfile: " + infile);
  }
  AmrData& amrData = dataServices.AmrDataRef();

  finest_level = std::min(finest_level, amrData.FinestLevel());
  const int n_lev = finest_level + 1;

  const Vector<std::string>& all_vars = amrData.PlotVarNames();
  for (const auto& name : var_names) {
    find_var_index(all_vars, name, /*abort_if_not_found=*/true);
  }

  Vector<int> is_per(AMREX_SPACEDIM, 0);
  if (!is_per_in.empty()) {
    AMREX_ALWAYS_ASSERT(static_cast<int>(is_per_in.size()) == AMREX_SPACEDIM);
    is_per = is_per_in;
  }

  RealBox rb(&(amrData.ProbLo()[0]), &(amrData.ProbHi()[0]));
  const int coord = 0;
  const int ncomp = static_cast<int>(var_names.size());
  Vector<int> dest_comps(ncomp);
  for (int i = 0; i < ncomp; ++i)
    dest_comps[i] = i;

  const int nprocs = ParallelDescriptor::NProcs();

  PlotfileData result;
  result.n_lev = n_lev;
  result.time = amrData.Time();
  result.var_names = var_names;
  result.mf.resize(n_lev);
  result.geoms.resize(n_lev);

  for (int lev = 0; lev < n_lev; ++lev) {
    // Guard: if the file has fewer boxes than MPI ranks, FillVar can
    // deadlock because some ranks own no FABs in the collective I/O.
    // Redistribute by splitting the BoxArray before calling FillVar.
    const BoxArray ba = (nprocs > 1)
                          ? ensure_min_boxes(amrData.boxArray(lev), nprocs)
                          : amrData.boxArray(lev);
    const DistributionMapping dm(ba);
    result.geoms[lev] =
      Geometry(amrData.ProbDomain()[lev], &rb, coord, is_per.data());
    result.mf[lev].define(ba, dm, ncomp, n_grow);
    amrData.FillVar(result.mf[lev], lev, var_names, dest_comps);
  }

  result.ref_ratios.resize(finest_level);
  for (int lev = 0; lev < finest_level; ++lev) {
    result.ref_ratios[lev] = amrData.RefRatio()[lev];
  }

  return result;
}

void
write_plotfile(
  const std::string& outfile,
  const Vector<MultiFab>& mf,
  const Vector<std::string>& var_names,
  const Vector<Geometry>& geoms,
  Real time,
  const Vector<int>& ref_ratios_in)
{
  // Synchronize all ranks before the collective write to avoid deadlock
  // when some ranks arrive later (e.g. after a preceding read_plotfile).
  ParallelDescriptor::Barrier();

  // Honour the `n_files` ParmParse option to cap the number of output files.
  set_plot_nfiles();

  const int n_lev = static_cast<int>(mf.size());
  Vector<int> isteps(n_lev, 0);
  Vector<IntVect> ref_ratios(std::max(n_lev - 1, 0));
  for (int lev = 0; lev < n_lev - 1; ++lev) {
    const int r = ref_ratios_in.empty() ? 2 : ref_ratios_in[lev];
    ref_ratios[lev] = IntVect(AMREX_D_DECL(r, r, r));
  }

  // Guard: WriteMultiLevelPlotfile can deadlock in MPI mode when a rank
  // owns no FABs (e.g. single-box domain with many ranks). Redistribute so
  // every rank has at least one box before calling the collective write.
  const int nprocs = ParallelDescriptor::NProcs();
  bool needs_split = false;
  if (nprocs > 1) {
    for (int lev = 0; lev < n_lev; ++lev) {
      if (static_cast<int>(mf[lev].boxArray().size()) < nprocs) {
        needs_split = true;
        break;
      }
    }
  }

  if (needs_split) {
    Vector<MultiFab> mf_split(n_lev);
    for (int lev = 0; lev < n_lev; ++lev) {
      BoxArray ba = ensure_min_boxes(mf[lev].boxArray(), nprocs);
      DistributionMapping dm(ba);
      mf_split[lev].define(ba, dm, mf[lev].nComp(), 0);
      mf_split[lev].ParallelCopy(mf[lev]);
    }
    WriteMultiLevelPlotfile(
      outfile, n_lev, GetVecOfConstPtrs(mf_split), var_names, geoms, time,
      isteps, ref_ratios);
  } else {
    WriteMultiLevelPlotfile(
      outfile, n_lev, GetVecOfConstPtrs(mf), var_names, geoms, time, isteps,
      ref_ratios);
  }
  // Ensure the Header written by IOProcessor is visible to all ranks
  // before any rank returns and potentially tries to read the same file.
  ParallelDescriptor::Barrier();
}

} // namespace analysis_util
