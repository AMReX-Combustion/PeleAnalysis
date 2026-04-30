#include <AMReX_MLMG.H>
#include <AMReX_MLPoisson.H>
#include <AMReX_MLABecLaplacian.H>
#include <analysis_util.H>
#include <AMReX_ParmParse.H>

namespace analysis_util {

void
gradient(
  const amrex::Vector<amrex::MultiFab>& a_mf,
  const amrex::Vector<amrex::Geometry>& geoms,
  const int                             scomp,
  const int                             ncomp,
  amrex::Vector<amrex::MultiFab>&       grad_mf)
{
  // Need AMREX_SPACEDIM components of gradients
  AMREX_ALWAYS_ASSERT(grad_mf[0].nComp() == AMREX_SPACEDIM * ncomp);
  // Need the data in to have at least one grow cell
  AMREX_ALWAYS_ASSERT(a_mf[0].nGrowVect() > amrex::IntVect(0));

  const int nlev = static_cast<int>(a_mf.size());
  amrex::Vector<amrex::BoxArray>          grids(nlev);
  amrex::Vector<amrex::DistributionMapping> dmap(nlev);
  for (int lev = 0; lev < nlev; ++lev) {
    grids[lev] = a_mf[lev].boxArray();
    dmap[lev]  = a_mf[lev].DistributionMap();
  }

  amrex::LPInfo info;
  info.setAgglomeration(1);
  info.setConsolidation(1);
  info.setMetricTerm(false);
  info.setMaxCoarseningLevel(0);
  amrex::MLPoisson poisson({geoms}, {grids}, {dmap}, info);
  poisson.setMaxOrder(4);
  amrex::Array<amrex::LinOpBCType, AMREX_SPACEDIM> lo_bc;
  amrex::Array<amrex::LinOpBCType, AMREX_SPACEDIM> hi_bc;

  amrex::ParmParse pp;
  amrex::Vector<int> sym_dir(AMREX_SPACEDIM, 0);
  pp.queryarr("sym_dir", sym_dir, 0, AMREX_SPACEDIM);

  const amrex::IntVect periodicity = geoms[0].periodicity().intVect();
  for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
    if (periodicity[idim] == 1) {
      lo_bc[idim] = hi_bc[idim] = amrex::LinOpBCType::Periodic;
    } else {
      if (sym_dir[idim] == 1) {
        lo_bc[idim] = hi_bc[idim] = amrex::LinOpBCType::reflect_odd;
      } else {
        lo_bc[idim] = hi_bc[idim] = amrex::LinOpBCType::Neumann;
      }
    }
  }
  poisson.setDomainBC(lo_bc, hi_bc);
  constexpr int nGrowGrad = 1;
  for (int n = 0; n < ncomp; ++n) {
    amrex::Vector<amrex::Array<amrex::MultiFab, AMREX_SPACEDIM>> grad(nlev);
    amrex::Vector<std::unique_ptr<amrex::MultiFab>>              phi(nlev);
    amrex::Vector<amrex::MultiFab>                               lap;
    lap.reserve(nlev);
    for (int lev = 0; lev < nlev; ++lev) {
      for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        grad[lev][idim].define(
          amrex::convert(grids[lev], amrex::IntVect::TheDimensionVector(idim)),
          dmap[lev], 1, nGrowGrad);
      }
      phi[lev] = std::make_unique<amrex::MultiFab>(
        a_mf[lev], amrex::make_alias, scomp + n, 1);
      poisson.setLevelBC(lev, phi[lev].get());
      lap.emplace_back(grids[lev], dmap[lev], 1, 1);
    }

    amrex::MLMG mlmg(poisson);
    mlmg.apply(GetVecOfPtrs(lap), GetVecOfPtrs(phi));
    mlmg.getFluxes(
      GetVecOfArrOfPtrs(grad), GetVecOfPtrs(phi),
      amrex::MLMG::Location::FaceCenter);
    for (int lev = 0; lev < nlev; ++lev) {
      // Convert to cell avg gradient
      amrex::MultiFab gradAlias(
        grad_mf[lev], amrex::make_alias, n * AMREX_SPACEDIM, AMREX_SPACEDIM);
      average_face_to_cellcenter(gradAlias, 0, GetArrOfConstPtrs(grad[lev]));
      gradAlias.mult(-1.0);
    }
  }
}

} // namespace analysis_util
