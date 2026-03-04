#include <AMReX_Geometry.H>
#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Utility.H>
#include <AMReX_MultiFab.H>

#include "StreamPC.H"

using namespace amrex;

StreamParticleContainer::StreamParticleContainer(
  int a_nPtsOnStrm,
  const Vector<Geometry>& a_geoms,
  const Vector<DistributionMapping>& a_dmaps,
  const Vector<BoxArray>& a_bas,
  const Vector<int>& a_rrs,
  const int a_nVar,
  const IntVect a_is_per,
  const Vector<std::string>& a_outVarNames)
  : ParticleContainer<0, 3, 0, 0>(a_geoms, a_dmaps, a_bas, a_rrs)
{
  Nlev = a_geoms.size();
  nPtsOnStrm = a_nPtsOnStrm;
  fcomp = a_nVar;
  pcomp = AMREX_SPACEDIM + fcomp;
  sizeOfRealStreamData = nPtsOnStrm * pcomp;
  is_per = a_is_per;
  outVarNames = a_outVarNames;
  ResizeRuntimeRealComp(sizeOfRealStreamData, true);
}

void
StreamParticleContainer::InitParticles(const Vector<Vector<Real>>& locs)
{
  BL_PROFILE("StreamParticleContainer::InitParticles");

  const int nProc = ParallelDescriptor::NProcs();
  const int myProc = ParallelDescriptor::MyProc();
  int nLocs = locs.size();
  auto& particle_tile = DefineAndReturnParticleTile(Nlev - 1, myProc, 0);
  for (int n = 0; n < nLocs; n++) {
    // lets just initialise on a random processor and redistribute afterwards
    // trying to find the right processor here doesn't work well with particles
    // near boundaries
    if (n % nProc == myProc) {
      // Create IDs for the two particles, NextID only works process-locally, so
      // lets manually assign instead
      const int p1_id = 2 * n + 1;
      const int p2_id = 2 * n + 2;
      // create a pair of particles
      std::pair<ParticleType, ParticleType> ppair;
      // Give them the two IDs
      ppair.first.id() = p1_id;
      ppair.second.id() = p2_id;
      // Initialise on same process
      ppair.first.cpu() = myProc;
      ppair.second.cpu() = myProc;
      // Initialise at same position
      AMREX_D_EXPR(
        ppair.first.pos(0) = locs[n][0], ppair.first.pos(1) = locs[n][1],
        ppair.first.pos(2) = locs[n][2]);
      AMREX_D_EXPR(
        ppair.second.pos(0) = locs[n][0], ppair.second.pos(1) = locs[n][1],
        ppair.second.pos(2) = locs[n][2]);
      // Initialise position at 0
      ppair.first.idata(0) = 0;
      ppair.second.idata(0) = 0;
      // First particle going forwards, second going backwards
      ppair.first.idata(1) = 1;
      ppair.second.idata(1) = -1;
      // Each particle needs to know the other one's ID
      ppair.first.idata(2) = p2_id;
      ppair.second.idata(2) = p1_id;
      // Add each particle to the tile and allocate real space in order
      particle_tile.push_back(ppair.first);
      for (int i = 0; i < NumRuntimeRealComps(); i++) {
        particle_tile.push_back_real(i, i < AMREX_SPACEDIM ? locs[n][i] : 0);
      }
      particle_tile.push_back(ppair.second);
      for (int i = 0; i < NumRuntimeRealComps(); i++) {
        particle_tile.push_back_real(i, i < AMREX_SPACEDIM ? locs[n][i] : 0);
      }
    }
  }
  Redistribute();
}

void
StreamParticleContainer::SetParticleLocation(const int a_streamLoc)
{
  BL_PROFILE("StreamParticleContainer::SetParticleLocation");
  // offset to find coordinates
  int offset = pcomp * a_streamLoc;
  dim3 newpos;
  for (int lev = 0; lev < Nlev; ++lev) {
    for (MyParIter pti(*this, lev); pti.isValid(); ++pti) {
      auto& aos = pti.GetArrayOfStructs();
      auto& soa = pti.GetStructOfArrays();

      for (size_t pindex = 0; pindex < aos.size(); ++pindex) {
        ParticleType& p = aos[pindex];
        if (p.id() > 0) {
          newpos = {AMREX_D_DECL(
            soa.GetRealData(offset + RealData::xloc)[pindex],
            soa.GetRealData(offset + RealData::yloc)[pindex],
            soa.GetRealData(offset + RealData::zloc)[pindex])};

          AMREX_D_EXPR(
            p.pos(0) = newpos[0], p.pos(1) = newpos[1], p.pos(2) = newpos[2]);
        }
      }
    }
  }
  Redistribute();
}

static void
VectorNormalize(Vector<Real>& vec, int dir)
{
  static Real eps = 1.e24;
  Real mag = std::sqrt(
    AMREX_D_TERM(vec[0] * vec[0], +vec[1] * vec[1], +vec[2] * vec[2]));
  if (mag < eps) {
    for (int i = 0; i < AMREX_SPACEDIM; ++i)
      vec[i] *= dir / mag;
  } else {
    vec = {AMREX_D_DECL(0, 0, 0)};
  }
}

static bool
InterpolateVector(
  const dim3 x,
  const FArrayBox& gfab,
  const Real* dx,
  const Real* plo,
  const Real* phi,
  Vector<Real>& u,
  int nComp)
{
  int3 b;
  dim3 n;

  for (int d = 0; d < AMREX_SPACEDIM; ++d) {
    b[d] = (int)std::floor((x[d] - plo[d]) / dx[d] - 0.5);
    n[d] = (x[d] - ((b[d] + 0.5) * dx[d] + plo[d])) / dx[d];
    n[d] = std::max(0., std::min(1., n[d]));
  }

  const auto& gbx = gfab.box();
  const auto& glo = gbx.smallEnd();
  const auto& ghi = gbx.bigEnd();
  for (int d = 0; d < AMREX_SPACEDIM; ++d) {
    if (b[d] < glo[d] || b[d] > ghi[d] - 1) {
      std::cout << "dir: " << d << std::endl;
      std::cout << "d,b,glo,ghi: " << d << " " << b[d] << " " << glo[d] << " "
                << ghi[d] << std::endl;
      std::cout << "x,plo,phi " << x[d] << " " << plo[d] << " " << phi[d]
                << std::endl;
      std::cout << "boxlo,boxhi " << plo[d] + glo[d] * dx[d] << " "
                << plo[d] + (ghi[d] + 1) * dx[d] << std::endl;
      return false;
    }
  }
  const auto& g = gfab.array();
  for (int i = 0; i < nComp; ++i) {
#if AMREX_SPACEDIM == 2
    u[i] = +n[0] * n[1] * g(b[0] + 1, b[1] + 1, 0, i) +
           n[0] * (1 - n[1]) * g(b[0] + 1, b[1], 0, i) +
           (1 - n[0]) * n[1] * g(b[0], b[1] + 1, 0, i) +
           (1 - n[0]) * (1 - n[1]) * g(b[0], b[1], 0, i);
#else
    u[i] = +n[0] * n[1] * n[2] * g(b[0] + 1, b[1] + 1, b[2] + 1, i) +
           n[0] * (1 - n[1]) * n[2] * g(b[0] + 1, b[1], b[2] + 1, i) +
           n[0] * n[1] * (1 - n[2]) * g(b[0] + 1, b[1] + 1, b[2], i) +
           n[0] * (1 - n[1]) * (1 - n[2]) * g(b[0] + 1, b[1], b[2], i) +
           (1 - n[0]) * n[1] * n[2] * g(b[0], b[1] + 1, b[2] + 1, i) +
           (1 - n[0]) * (1 - n[1]) * n[2] * g(b[0], b[1], b[2] + 1, i) +
           (1 - n[0]) * n[1] * (1 - n[2]) * g(b[0], b[1] + 1, b[2], i) +
           (1 - n[0]) * (1 - n[1]) * (1 - n[2]) * g(b[0], b[1], b[2], i);
#endif
  }
  return true;
}

bool
StreamParticleContainer::RungeKutta4(
  dim3& x,
  Real hrk,
  const FArrayBox& v,
  const Real* dx,
  const Real* plo,
  const Real* phi,
  int dir,
  int cSpace,
  Real dxFine)
{
  Vector<Real> vec(AMREX_SPACEDIM);
  dim3 k1, k2, k3, k4;
  dim3 xx = x;
  if (!InterpolateVector(xx, v, dx, plo, phi, vec, AMREX_SPACEDIM))
    return false;
  // before normalising vec, lets grab |grad vec|

  // in physical space, dt = hrk * dxFine (hrk <= 0.5)
  // in cspace, dt = min(0.95*dxLev,hrk/|gradC|) - min stops going too far in
  // one step and breaking things

  Real dt;
  if (cSpace == 0) {
    dt = hrk * dxFine;
  } else {
    Real abs_grad = std::sqrt(
      AMREX_D_TERM(vec[0] * vec[0], +vec[1] * vec[1], +vec[2] * vec[2]));
    dt = min(0.95 * dx[0], hrk / abs_grad);
  }

  // convert gradC to gradC / |gradC| (will also clip to zero if barely any
  // gradient)
  VectorNormalize(vec, dir);

  for (int d = 0; d < AMREX_SPACEDIM; ++d) {
    k1[d] = vec[d] * dt;
    xx[d] = x[d] + k1[d] * 0.5;
  }
  if (!InterpolateVector(xx, v, dx, plo, phi, vec, AMREX_SPACEDIM))
    return false;
  // update dt if using cspace
  if (cSpace != 0) {
    dt = min(
      0.95 * dx[0],
      hrk / std::sqrt(AMREX_D_TERM(
              vec[0] * vec[0], +vec[1] * vec[1], +vec[2] * vec[2])));
  }
  VectorNormalize(vec, dir);

  for (int d = 0; d < AMREX_SPACEDIM; ++d) {
    k2[d] = vec[d] * dt;
    xx[d] = x[d] + k2[d] * 0.5;
  }
  if (!InterpolateVector(xx, v, dx, plo, phi, vec, AMREX_SPACEDIM))
    return false;
  if (cSpace != 0) {
    dt = min(
      0.95 * dx[0],
      hrk / std::sqrt(AMREX_D_TERM(
              vec[0] * vec[0], +vec[1] * vec[1], +vec[2] * vec[2])));
  }
  VectorNormalize(vec, dir);

  for (int d = 0; d < AMREX_SPACEDIM; ++d) {
    k3[d] = vec[d] * dt;
    xx[d] = x[d] + k3[d];
  }
  if (!InterpolateVector(xx, v, dx, plo, phi, vec, AMREX_SPACEDIM))
    return false;

  if (cSpace != 0) {
    dt = min(
      0.95 * dx[0],
      hrk / std::sqrt(AMREX_D_TERM(
              vec[0] * vec[0], +vec[1] * vec[1], +vec[2] * vec[2])));
  }
  VectorNormalize(vec, dir);

  const Real third = 1. / 3.;
  const Real sixth = 1. / 6.;
  dim3 delta;
  for (int d = 0; d < AMREX_SPACEDIM; ++d) {
    k4[d] = vec[d] * dt;
    delta[d] = (k1[d] + k4[d]) * sixth + (k2[d] + k3[d]) * third;
  }

  // cut step length to keep in domain (FIXME: Deal with periodic - hopefully
  // fixed?)
  Real scale = 1;
  for (int d = 0; d < AMREX_SPACEDIM; ++d) {
    if ((x[d] + delta[d] <= plo[d] + dx[d]) && (is_per[d] == 0)) {
      // std::cout << "Too low, clipping path: x=" << x[d] << ", delta = " <<
      // delta[d] << std::endl;
      scale = 0.0; // stop moving
      // scale = std::min(scale, std::abs((x[d] - (plo[d]+dx[d]))/delta[d]));
    }
    if ((x[d] + delta[d] >= phi[d] - dx[d]) && (is_per[d] == 0)) {
      // std::cout << "Too high, clipping path: x=" << x[d] << ", delta = " <<
      // delta[d] << std::endl;
      scale = 0.0; // just stop
      // scale = std::min(scale, std::abs(((phi[d]-dx[d]) - x[d])/delta[d]));
    }
  }
  for (int d = 0; d < AMREX_SPACEDIM; ++d) {
    x[d] += scale * delta[d];
    x[d] = std::min(
      phi[d] - 1.e-10,
      std::max(plo[d] + 1.e-10, x[d])); // Deal with precision issues
  }
  return true;
}

void
StreamParticleContainer::ComputeNextLocation(
  int a_fromLoc,
  const Real a_hRK,
  const Vector<MultiFab>& a_vectorField,
  const int a_cSpace)
{
  BL_PROFILE("StreamParticleContainer::ComputeNextLocation");

  SetParticleLocation(a_fromLoc);
  // Real dt = a_hRK;
  Real dxFine = Geom(Nlev - 1).CellSize()[0];

  const int new_loc_id = a_fromLoc + 1;
  int offset = pcomp * new_loc_id;

  for (int lev = 0; lev < Nlev; ++lev) {
    const auto& geom = Geom(lev);
    const auto& dx = geom.CellSize();
    const auto& plo = geom.ProbLo();
    const auto& phi = geom.ProbHi();

    for (MyParIter pti(*this, lev); pti.isValid(); ++pti) {
      auto& aos = pti.GetArrayOfStructs();
      auto& soa = pti.GetStructOfArrays();
      const FArrayBox& v = a_vectorField[lev][pti];

      for (size_t pindex = 0; pindex < aos.size(); ++pindex) {
        ParticleType& p = aos[pindex];
        const int dir = p.idata(1);
        dim3 x = {AMREX_D_DECL(p.pos(0), p.pos(1), p.pos(2))};
        if (p.id() > 0) {
          if (!RungeKutta4(x, a_hRK, v, dx, plo, phi, dir, a_cSpace, dxFine)) {
            Abort("bad RK");
          }
        }
        AMREX_D_EXPR(
          soa.GetRealData(offset + RealData::xloc)[pindex] = x[0],
          soa.GetRealData(offset + RealData::yloc)[pindex] = x[1],
          soa.GetRealData(offset + RealData::zloc)[pindex] = x[2]);
      }
    }
  }
}

void
StreamParticleContainer::InterpDataAtLocation(
  int a_fromLoc, const Vector<MultiFab>& a_vectorField)
{
  BL_PROFILE("StreamParticleContainer::InterpDataAtLocation");

  SetParticleLocation(a_fromLoc);

  int offset = pcomp * a_fromLoc; // components on particle

  for (int lev = 0; lev < Nlev; ++lev) {
    const auto& geom = Geom(lev);
    const auto& dx = geom.CellSize();
    const auto& plo = geom.ProbLo();
    const auto& phi = geom.ProbHi();

    dim3 Lx;
    for (int iComp = 0; iComp < AMREX_SPACEDIM; iComp++)
      Lx[iComp] = phi[iComp] - plo[iComp];

    for (MyParIter pti(*this, lev); pti.isValid(); ++pti) {
      auto& aos = pti.GetArrayOfStructs();
      auto& soa = pti.GetStructOfArrays();
      const FArrayBox& v = a_vectorField[lev][pti];

      for (size_t pindex = 0; pindex < aos.size(); ++pindex) {
        ParticleType& p = aos[pindex];
        if (p.id() > 0) {
          // where's the particle?
          dim3 x = {AMREX_D_DECL(p.pos(0), p.pos(1), p.pos(2))};
          Vector<Real> ntrpvOut(fcomp); // components in infile

          // interpolate all data to particle location
          InterpolateVector(
            x, v, dx, plo, phi, ntrpvOut, fcomp); // components in infile

          // copy the interpolated data to the particle
          // first DIM components are particle location
          // next DIM components are stream location w/o adjusting for
          // periodicity then we have the interpolated data we want
          for (int iComp = AMREX_SPACEDIM; iComp < fcomp; ++iComp) {
            int idxOnPart = offset + iComp + AMREX_SPACEDIM;
            soa.GetRealData(idxOnPart)[pindex] = ntrpvOut[iComp];
          }
          // let's figure out the location on the stream w/o adjusting for
          // periodicity
          if (a_fromLoc == 0) { // nothing to do on the surface; just take a
                                // copy of the location
            for (int d = 0; d < AMREX_SPACEDIM; ++d) {
              int idx = offset + d;
              soa.GetRealData(idx + AMREX_SPACEDIM)[pindex] =
                soa.GetRealData(idx)[pindex];
            }
          } else { // set new location adjusting for periodicity
            // calculate the change in position delta
            // add to the old position
            // adjust delta if it's affected by periodicity (i.e. delta too big)
            for (int d = 0; d < AMREX_SPACEDIM; ++d) {
              int idx = offset + d;
              int idxOld = idx - pcomp;
              Real xnew = soa.GetRealData(idx)[pindex];
              Real xold = soa.GetRealData(idxOld)[pindex];
              Real delta = xnew - xold;
              if (fabs(delta) > dx[d]) { // has been adjusted for periodicity
                // printf("%i %i %e %e %e",pindex,d,xnew,xold,delta);
                if (delta < 0.)
                  delta += Lx[d];
                else
                  delta -= Lx[d];
                // printf(" --> %e\n",delta);
              }
              // store periodicity-adjusted copy of location at idx+SPACEDIM
              Real sold = soa.GetRealData(idxOld + AMREX_SPACEDIM)[pindex];
              Real snew = sold + delta;
              soa.GetRealData(idx + AMREX_SPACEDIM)[pindex] = snew;
            }
          }
          // hack last component to AMR level for debugging
          // soa.GetRealData(offset + DEF_FCOMP-1)[pindex] = (Real)lev;
        }
      }
    }
  }
}

void
StreamParticleContainer::WriteStreamAsTecplot(const std::string& outFile)
{
  // Set location to first point on stream to guarantee partner line is local
  SetParticleLocation(0);

  // Create a folder and have each processor write their own data, one file per
  // streamline
  auto myProc = ParallelDescriptor::MyProc();

  if (!amrex::UtilCreateDirectory(outFile, 0755))
    amrex::CreateDirectoryFailed(outFile);
  ParallelDescriptor::Barrier();

  bool will_write = false;
  for (int lev = 0; lev < Nlev && !will_write; ++lev) {
    for (MyParIter pti(*this, lev); pti.isValid() && !will_write; ++pti) {
      auto& aos = pti.GetArrayOfStructs();

      for (size_t pindex = 0; pindex < aos.size() && !will_write; ++pindex) {
        ParticleType& p = aos[pindex];
        will_write |= (p.id() > 0);
      }
    }
  }

  if (will_write) {
    std::string fileName = outFile + "/str_";
    fileName = Concatenate(fileName, myProc) + ".dat";
    std::ofstream ofs(fileName.c_str());
    ofs << "VARIABLES = \"";
    for (int iComp = 0; iComp < fcomp - 1; ++iComp)
      ofs << outVarNames[iComp] << "\" \"";
    ofs << outVarNames[fcomp - 1] << "\"\n";

    for (int lev = 0; lev < Nlev; ++lev) {
      for (MyParIter pti(*this, lev); pti.isValid(); ++pti) {
        auto& aos = pti.GetArrayOfStructs();
        auto& soa = pti.GetStructOfArrays();

        for (size_t pindex = 0; pindex < aos.size(); ++pindex) {
          ofs << "ZONE I=1 J=" << nPtsOnStrm << " K=1 FORMAT=POINT\n";
          for (int j = 0; j < nPtsOnStrm; ++j) {
            // by including the spacedim offset, we use locations w/o
            // periodicity adjustments
            int offset = j * pcomp + AMREX_SPACEDIM;
            // std::cout << offset << std::endl;
            for (int iComp = 0; iComp < fcomp; ++iComp) {
              ofs << soa.GetRealData(offset + iComp)[pindex] << " ";
            }
            ofs << '\n';
          }
        }
      }
    }
    ofs.close();
  }
}

#if AMREX_DEBUG // AJA inspect function
void
StreamParticleContainer::InspectParticles(const int nStreamPairs)
{
  BL_PROFILE("StreamParticleContainer::InspectParticles");

  SetParticleLocation(0);

  int myProc = ParallelDescriptor::MyProc();
  int nProcs = ParallelDescriptor::NProcs();

  Vector<Vector<int>> partIndexing(nProcs);
  Vector<Vector<int>> ptiIndexing(nProcs);

  for (int iProc = 0; iProc < nProcs; iProc++) {
    partIndexing[iProc].resize(2 * nStreamPairs + 1);
    ptiIndexing[iProc].resize(2 * nStreamPairs + 1);
    for (int iStream = 0; iStream <= 2 * nStreamPairs; iStream++) {
      partIndexing[iProc][iStream] = -1;
      ptiIndexing[iProc][iStream] = -1;
    }
  }

  int ptiCounter;

  int minId = 100000000;
  int maxId = -minId;

  for (int lev = 0; lev < Nlev; ++lev) {
    ptiCounter = 0;

    for (MyParIter pti(*this, lev); pti.isValid(); ++pti, ++ptiCounter) {

      auto& aos = pti.GetArrayOfStructs();

      for (size_t pindex = 0; pindex < aos.size(); ++pindex) {
        ParticleType& p = aos[pindex];

        if (p.id() > 0) {
          int pId = p.id();
          minId = min(minId, pId);
          maxId = max(maxId, pId);
          if (pId > 2 * nStreamPairs) {
            std::cout << "pId > 2*nStreamPairs " << pId << " > "
                      << (2 * nStreamPairs) << std::endl;
          }
          partIndexing[myProc][pId] = pindex;
          ptiIndexing[myProc][pId] = ptiCounter;
        } else {
          Print() << "p.id() = " << p.id() << " !> 0" << std::endl;
        } // p.id

      } // pindex
    } // pti

    // Now let's see if the pair is on the same processor and pti

    ptiCounter = 0; // reset

    for (MyParIter pti(*this, lev); pti.isValid(); ++pti, ++ptiCounter) {
      auto& aos = pti.GetArrayOfStructs();
      int goodCount = 0;
      for (int iStream = 1; iStream <= 2 * nStreamPairs; iStream++) {
        int pindex1 = partIndexing[myProc][iStream];
        int ptiindex1 = ptiIndexing[myProc][iStream];
        // if ( (pindex1>-1) && (pindex1<aos.size()) ) {
        if ((pindex1 > -1) && (ptiindex1 == ptiCounter)) {
          ParticleType& p1 = aos[pindex1];
          int myId1 = p1.id();
          int pairId1 = p1.idata(2);
          int pindex2 = partIndexing[myProc][pairId1];
          ParticleType& p2 = aos[pindex2];
          int myId2 = p2.id();
          int pairId2 = p2.idata(2);
          if ((myId1 != pairId2) || (pairId1 != myId2)) {
            std::cout << "myProc / pindex1 / myId1 / pairId1 = " << myProc
                      << " / " << pindex1 << " / " << myId1 << " / " << pairId1
                      << " / " << std::endl;
            std::cout << "myProc / pindex2 / myId2 / pairId2 = " << myProc
                      << " / " << pindex2 << " / " << myId2 << " / " << pairId2
                      << " / " << std::endl;
          } else {
            // looks good
            goodCount += 2;
          }
        }
      }
      if (goodCount != (int)aos.size()) {
        std::cout << "aos.size() != goodCount : " << aos.size()
                  << " != " << goodCount << std::endl;
      }
    }

  } // lev

  ParallelDescriptor::ReduceIntMin(minId);
  ParallelDescriptor::ReduceIntMax(maxId);
  if (ParallelDescriptor::IOProcessor())
    printf("InspectParticles: minId / maxId = %i / %i\n", minId, maxId);

  int duplicate(0);
  int missing(0);
  for (int iStream = 1; iStream <= 2 * nStreamPairs; iStream++) {
    int pIdx = partIndexing[myProc][iStream];
    int ptiIdx = ptiIndexing[myProc][iStream];
    int gotcha = 0;
    if ((pIdx != -1) && (ptiIdx != -1)) {
      gotcha++;
    }
    ParallelDescriptor::ReduceIntSum(gotcha);
    if (gotcha == 0) {
      Print() << "iStream(" << iStream << ") == 0" << std::endl;
      missing++;
    }
    if (gotcha > 1) {
      Print() << "iStream(" << iStream << ") > 1" << std::endl;
      duplicate++;
    }
  }
  if (missing > 0)
    Print() << "Missing " << missing << " particles !!!" << std::endl;
  else {
    Print() << "No particles missing :-)" << std::endl;
  }
  if (duplicate > 0)
    Print() << "There are " << duplicate << " duplicates !!!" << std::endl;
  else {
    Print() << "No duplicate particles :-)" << std::endl;
  }
}

#endif

void
StreamParticleContainer::WriteStreamAsBinary(
  const std::string& outFile, Vector<int>& faceData)
{
  // Set location to first point on stream to guarantee partner line is local
  SetParticleLocation(0);

  // Create a folder and have each processor write their own data, one file per
  // streamline
  auto myProc = ParallelDescriptor::MyProc();
  auto nProcs = ParallelDescriptor::NProcs();

  if (!amrex::UtilCreateDirectory(outFile, 0755))
    amrex::CreateDirectoryFailed(outFile);
  ParallelDescriptor::Barrier();

  bool will_write = false;
  for (int lev = 0; lev < Nlev && !will_write; ++lev) {
    for (MyParIter pti(*this, lev); pti.isValid() && !will_write; ++pti) {
      auto& aos = pti.GetArrayOfStructs();

      for (size_t pindex = 0; pindex < aos.size() && !will_write; ++pindex) {
        ParticleType& p = aos[pindex];
        will_write |= (p.id() > 0);
      }
    }
  }

  // Need to count the total number of streams to be written
  // by all ptiters on all levels on this processor
  int nStreams = 0;
  int nStreamsCheck = 0;

  for (int lev = 0; lev < Nlev; ++lev) {
    for (MyParIter pti(*this, lev); pti.isValid(); ++pti) {
      auto& aos = pti.GetArrayOfStructs();
      for (size_t pindex = 0; pindex < aos.size(); ++pindex) {
        ParticleType& p = aos[pindex];
        if ((p.id() > 0) && (p.idata(1) == 1)) {
          nStreams++;
        }
      }
    }
  }

  // write to a binary file
  std::string rootName = outFile + "/str_";
  std::string fileName = Concatenate(rootName, myProc) + ".bin";
  // std::string headName = Concatenate(rootName,myProc) + ".head";
  FILE* file = fopen(fileName.c_str(), "w");
  // FILE *head=fopen(headName.c_str(),"w");
  //  total number of streams in file
  fwrite(&(nStreams), sizeof(int), 1, file);

  int minId = 100000000;
  int maxId = -minId;
  for (int lev = 0; lev < Nlev; ++lev) {

    for (MyParIter pti(*this, lev); pti.isValid(); ++pti) {
      auto& aos = pti.GetArrayOfStructs();
      auto& soa = pti.GetStructOfArrays();

      int aosSize = aos.size();

      // construct the pair mapping for this pti
      // Vector<int> pIdMap(2*nStreamPairs+1);
      std::map<int, int> pid_to_pindex;
      for (int pindex = 0; pindex < aosSize; ++pindex) {
        ParticleType& p = aos[pindex];
        int pId = p.id();
        if (pId > 0) {
          pid_to_pindex[pId] = pindex;
        }
      }

      for (int pindex = 0; pindex < aosSize; ++pindex) {
        ParticleType& p = aos[pindex];
        int pId = p.id();
        int dir = p.idata(1);
        int pairId = p.idata(2);

        if ((pId > 0) && (dir == 1)) {
          // write info about this stream and its pair
          int pindexPair = pid_to_pindex[pairId];
          ParticleType& pPair = aos[pindexPair];
          int pIdPair = pPair.id();
          int pairIdPair = pPair.idata(2);

          // sanity check
          if ((pIdPair != pairId) || (pId != pairIdPair)) {
            std::cout << pIdPair << std::endl;
            std::cout << pairId << std::endl;
            std::cout << pId << std::endl;
            std::cout << pairIdPair << std::endl;
            Abort("Bad pair mapping");
          }
          // back out the original id from the surface (both count from 1)
          int pIdInv = (pId + 1) / 2;
          fwrite(&(pIdInv), sizeof(int), 1, file); // pair id

          // fprintf(head,"%i %i %i\n",pId,dir,pairId);

          minId = min(minId, pId);
          maxId = max(maxId, pId);

          // if we write in the order paricle->component->position on surface,
          // then end up with a single-component stream together in memory
          // by including spacedim in offset, we use locations w/o periodicity
          // adjustments

          for (int iComp = 0; iComp < fcomp; iComp++) {
            // first write pair (dir=-1) backwards (without double writing
            // surface)
            for (int j = nPtsOnStrm - 1; j > 0; j--) {
              int offset = j * pcomp + iComp + AMREX_SPACEDIM;
              fwrite(
                &(soa.GetRealData(offset)[pindexPair]), sizeof(Real), 1, file);
            }
            // now write this particle's (dir=1) data (with surface)
            for (int j = 0; j < nPtsOnStrm; j++) {
              int offset = j * pcomp + iComp + AMREX_SPACEDIM;
              fwrite(&(soa.GetRealData(offset)[pindex]), sizeof(Real), 1, file);
            }
          }

          nStreamsCheck++; // sanity check to make sure we wrote number of
                           // streams anticipated
        }
      }
    }
  }
  ParallelDescriptor::ReduceIntMin(minId);
  ParallelDescriptor::ReduceIntMax(maxId);

  if (nStreams != nStreamsCheck)
    std::cout << "(nStreams!=nStreamsCheck) : " << nStreams
              << " != " << nStreamsCheck << std::endl;

  fclose(file);

  // Write header file with everything consistent across all processors
  ParallelDescriptor::ReduceIntSum(nStreams);
  if (ParallelDescriptor::IOProcessor()) {
    fileName = outFile + "/Header";
    std::ofstream ofs(fileName.c_str());
    ofs << "Even odder-ball replacement for sampled streams" << std::endl;
    ofs << nProcs << std::endl;   // translates to number of files to read
    ofs << nStreams << std::endl; // total number of streams
    ofs << 2 * nPtsOnStrm - 1 << std::endl; // number of points
    ofs << fcomp << std::endl;              // number of variables
    for (int iComp = 0; iComp < fcomp; ++iComp)
      ofs << outVarNames[iComp] << " ";
    ofs << std::endl;
    ofs << faceData.size() << std::endl;
    ofs.write((char*)faceData.dataPtr(), sizeof(int) * faceData.size());
    ofs << '\n';
    ofs.close();
  }
}
