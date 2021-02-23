#include <AMReX_Geometry.H>
#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Utility.H>
#include <AMReX_MultiFab.H>

#include "StreamPC.H"

using namespace amrex;

typedef Array<Real,AMREX_SPACEDIM> dim3;
typedef Array<int,AMREX_SPACEDIM> int3;

StreamParticleContainer::
StreamParticleContainer(int                                 a_nPtsOnStrm,
                        const Vector<Geometry>            & a_geoms,
                        const Vector<DistributionMapping> & a_dmaps,
                        const Vector<BoxArray>            & a_bas,
                        const Vector<int>                 & a_rrs)
  : ParticleContainer<0, 3, 0, 0> (a_geoms, a_dmaps, a_bas, a_rrs)
{
  Nlev = a_geoms.size();
  nPtsOnStrm = a_nPtsOnStrm;
  sizeOfRealStreamData = nPtsOnStrm * RealData::ncomp;
  for (int i=0; i<sizeOfRealStreamData; ++i) {
    AddRealComp(true);
  }
  for (int lev=0; lev<numLevels(); ++lev)
  {
    for (MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi)
    {
      auto& particle_tile = DefineAndReturnParticleTile(lev, mfi.index(), mfi.LocalTileIndex());
    }
  }
}

void
StreamParticleContainer::
InitParticles (const Vector<Vector<Real>>& locs)
{
  BL_PROFILE("StreamParticleContainer::InitParticles");

  int streamLoc = 0;
  int offset = RealData::ncomp*streamLoc;

  int lev = 0;
  int grid_id = 0;
  int tile_id = 0;
  int owner = ParticleDistributionMap(lev)[0];
  auto& particle_tile = GetParticles(0)[std::make_pair(grid_id,tile_id)];

  if (ParallelDescriptor::MyProc() == owner)
  {
    for (const auto& loc : locs)
    {
      // Keep track of pairs of lines
      Array<Long,2> ppair = {ParticleType::NextID(), ParticleType::NextID()};

      for (int i_part=0; i_part<2; i_part++)
      {
        ParticleType p;
        p.id()  = ppair[i_part];
        p.cpu() = ParallelDescriptor::MyProc();

        AMREX_D_EXPR(p.pos(0) = loc[0],
                     p.pos(1) = loc[1],
                     p.pos(2) = loc[2]);

        p.idata(0) = streamLoc;                  // Current position
        p.idata(1) = i_part==0 ? +1 : -1;        // Direction of integration
        p.idata(2) = ppair[ i_part==0 ? 1 : 0];  // Other line from this seed

        particle_tile.push_back(p);

        auto& soa = particle_tile.GetStructOfArrays();
        for (int i=0; i<NumRuntimeRealComps(); ++i)
        {
          soa.GetRealData(i).push_back(i<AMREX_SPACEDIM ? loc[i] : 0);
        }
      }
    }
  }
  Redistribute();
}

void
StreamParticleContainer::
SetParticleLocation(int a_streamLoc, int a_nGrow)
{
  BL_PROFILE("StreamParticleContainer::SetParticleLocation");

  AMREX_ALWAYS_ASSERT(a_nGrow > 0);
  bool redist = false;
  int offset = RealData::ncomp * a_streamLoc;
  dim3 newpos, blo, bhi;
  for (int lev = 0; lev < Nlev; ++lev)
  {
    const auto& geom = Geom(lev);
    const auto& dx = geom.CellSize();
    const auto& plo = geom.ProbLo();

    for (MyParIter pti(*this, lev); pti.isValid(); ++pti)
    {
      auto& aos = pti.GetArrayOfStructs();
      auto& soa = pti.GetStructOfArrays();
      const auto vbx = grow(pti.validbox(),a_nGrow-1);
      const auto& vse = vbx.smallEnd();
      const auto& vbe = vbx.bigEnd();
      blo = {AMREX_D_DECL(plo[0] + vse[0]*dx[0],
                          plo[1] + vse[1]*dx[1],
                          plo[2] + vse[2]*dx[2])};

      bhi = {AMREX_D_DECL(plo[0] + (vbe[0]+1)*dx[0],
                          plo[1] + (vbe[1]+1)*dx[1],
                          plo[2] + (vbe[2]+1)*dx[2])};

      for (int pindex=0; pindex<aos.size(); ++pindex)
      {
        ParticleType& p = aos[pindex];
        if (p.id() > 0)
        {
          newpos = {AMREX_D_DECL(soa.GetRealData(offset + RealData::xloc)[pindex],
                                 soa.GetRealData(offset + RealData::yloc)[pindex],
                                 soa.GetRealData(offset + RealData::zloc)[pindex])};

          AMREX_D_EXPR(p.pos(0) = newpos[0], p.pos(1) = newpos[1], p.pos(2) = newpos[2]);

          for (int d=0; d<AMREX_SPACEDIM; ++d)
          {
            redist |= (newpos[d]<blo[d] || newpos[d]>bhi[d]);
          }
        }
      }
    }
  }
  ParallelDescriptor::ReduceBoolOr(redist);
  if (redist) {
    //Print() << "  redistributing" << std::endl;
    Redistribute();
  }
}

static void vnrml(dim3& vec, int dir)
{
  static Real eps = 1.e12;
  Real sum = AMREX_D_TERM(vec[0] * vec[0],
                          + vec[1] * vec[1],
                          + vec[2] * vec[2]);
  dim3 u = vec;
  if (sum < eps) {
    sum = 1. / std::sqrt(sum);
    for (int i=0; i<AMREX_SPACEDIM; ++i) vec[i] *= dir * sum;
  }
  else {
    vec = {AMREX_D_DECL(0, 0, 0)};
  }
}

static bool ntrpv(const dim3& x,const FArrayBox& gfab,
                  const Real* dx,const Real* plo,const Real* phi,dim3& u)
{
  int3 b;
  dim3 n;

  for (int d=0; d<AMREX_SPACEDIM; ++d) {
    b[d] = std::floor( (x[d] - plo[d]) / dx[d] - 0.5 );
    n[d] = ( x[d] - ( (b[d] + 0.5 ) * dx[d] + plo[d] ) )/dx[d];
    n[d] = std::max(0., std::min(1.,n[d]));
  }

  const auto& gbx = gfab.box();
  const auto& glo = gbx.smallEnd();
  const auto& ghi = gbx.bigEnd();
  for (int d=0; d<AMREX_SPACEDIM; ++d) {
    if (b[d] < glo[d] ||  b[d] > ghi[d]-1) {
      Print() << "dir: " << d << std::endl;
      Print() << "d,b,glo,ghi: " << d << " " << b[d] << " " << glo[d] << " " << ghi[d] << std::endl;
      Print() << "x,plo,phi " << x[d] << " " << plo[d] << " " << phi[d] << std::endl;
      Print() << "boxlo,boxhi " << plo[d]+glo[d]*dx[d] << " " << plo[d]+(ghi[d]+1)*dx[d] << std::endl;
      return false;
    }
  }

  const auto& g = gfab.array();
  for (int i=0; i<AMREX_SPACEDIM; ++i) {
#if AMREX_SPACEDIM == 2
    u[i] =
      +   n[0]   *   n[1]   * g(b[0]+1,b[1]+1,0,i)
      +   n[0]   * (1-n[1]) * g(b[0]+1,b[1]  ,0,i)
      + (1-n[0]) *   n[1]   * g(b[0]  ,b[1]+1,0,i)
      + (1-n[0]) * (1-n[1]) * g(b[0]  ,b[1]  ,0,i);
#else
    u[i] =
      +    n[0]   *    n[1]  *    n[2]  * g(b[0]+1,b[1]+1,b[2]+1,i)
      +    n[0]   * (1-n[1]) *    n[2]  * g(b[0]+1,b[1]  ,b[2]+1,i)
      +    n[0]   *    n[1]  * (1-n[2]) * g(b[0]+1,b[1]+1,b[2]  ,i)
      +    n[0]   * (1-n[1]) * (1-n[2]) * g(b[0]+1,b[1]  ,b[2]  ,i)
      +  (1-n[0]) *    n[1]  *    n[2]  * g(b[0]  ,b[1]+1,b[2]+1,i)
      +  (1-n[0]) * (1-n[1]) *    n[2]  * g(b[0]  ,b[1]  ,b[2]+1,i)
      +  (1-n[0]) *    n[1]  * (1-n[2]) * g(b[0]  ,b[1]+1,b[2]  ,i)
      +  (1-n[0]) * (1-n[1]) * (1-n[2]) * g(b[0]  ,b[1]  ,b[2]  ,i);
#endif
  }
  return true;
}

static bool
RK4(dim3 & x,Real dt,const FArrayBox& v,const Real* dx,const Real* plo,const Real* phi,int dir)
{
  dim3 vec, k1, k2, k3, k4;
  dim3 xx = x;
  if (!ntrpv(xx,v,dx,plo,phi,vec)) return false;
  vnrml(vec,dir);

  for (int d=0; d<AMREX_SPACEDIM; ++d) {
    k1[d] = vec[d] * dt;
    xx[d] = x[d] + k1[d] * 0.5;
  }
  if (!ntrpv(xx,v,dx,plo,phi,vec)) return false;
  vnrml(vec,dir);

  for (int d=0; d<AMREX_SPACEDIM; ++d) {
    k2[d] = vec[d] * dt;
    xx[d] = x[d] + k2[d] * 0.5;
  }
  if (!ntrpv(xx,v,dx,plo,phi,vec)) return false;
  vnrml(vec,dir);

  for (int d=0; d<AMREX_SPACEDIM; ++d) {
    k3[d] = vec[d] * dt;
    xx[d] = x[d] + k3[d];
  }
  if (!ntrpv(xx,v,dx,plo,phi,vec)) return false;
  vnrml(vec,dir);

  const Real third = 1./3.;
  const Real sixth = 1./6.;
  dim3 delta;
  for (int d=0; d<AMREX_SPACEDIM; ++d) {
    k4[d] = vec[d] * dt;
    delta[d] = (k1[d] + k4[d])*sixth + (k2[d] + k3[d])*third;
  }

  // cut step length to keep in domain (FIXME: Deal with periodic)
  Real scale = 1;
  for (int d=0; d<AMREX_SPACEDIM; ++d) {
    if (x[d]+delta[d] < plo[d]) {
      scale = std::min(scale, std::abs((x[d] - plo[d])/delta[d]));
    }
    if (x[d]+delta[d] > plo[d]) {
      scale = std::min(scale, std::abs((phi[d] - x[d])/delta[d]));
    }
  }
  for (int d=0; d<AMREX_SPACEDIM; ++d) {
    x[d] += scale * delta[d];
    x[d] = std::min(phi[d]-1.e-10, std::max(plo[d]+1.e-10, x[d]) ); // Deal with precision issues
  }
  return true;
}

void
StreamParticleContainer::
ComputeNextLocation(int                      a_fromLoc,
                    Real                     a_delta_t,
                    const Vector<MultiFab> & a_vectorField)
{
  BL_PROFILE("StreamParticleContainer::ComputeNextLocation");

  const int nGrow = a_vectorField[0].nGrow();
  SetParticleLocation(a_fromLoc,nGrow);

  const int new_loc_id = a_fromLoc + 1;
  int offset = RealData::ncomp * new_loc_id;

  for (int lev = 0; lev < Nlev; ++lev)
  {
    const auto& geom = Geom(lev);
    const auto& dx = geom.CellSize();
    const auto& plo = geom.ProbLo();
    const auto& phi = geom.ProbHi();

    for (MyParIter pti(*this, lev); pti.isValid(); ++pti)
    {
      auto& aos = pti.GetArrayOfStructs();
      auto& soa = pti.GetStructOfArrays();
      const FArrayBox& v = a_vectorField[lev][pti];

      for (int pindex=0; pindex<aos.size(); ++pindex)
      {
        ParticleType& p = aos[pindex];
        const int dir = p.idata(1);
        dim3 x = {AMREX_D_DECL(p.pos(0), p.pos(1), p.pos(2))};
        if (p.id() > 0)
        {
          if (!RK4(x,a_delta_t,v,dx,plo,phi,dir))
          {
            Abort("bad RK");
          }
        }
        AMREX_D_EXPR(soa.GetRealData(offset + RealData::xloc)[pindex] = x[0],
                     soa.GetRealData(offset + RealData::yloc)[pindex] = x[1],
                     soa.GetRealData(offset + RealData::zloc)[pindex] = x[2]);
      }
    }
  }
}

void
StreamParticleContainer::
WriteStreamAsTecplot(const std::string& outFile)
{
  // Set location to first point on stream to guarantee partner line is local
  SetParticleLocation(0,1);

  // Create a folder and have each processor write their own data, one file per streamline
  auto myProc = ParallelDescriptor::MyProc();
  auto nProcs = ParallelDescriptor::NProcs();

  if (!amrex::UtilCreateDirectory(outFile, 0755))
    amrex::CreateDirectoryFailed(outFile);
  ParallelDescriptor::Barrier();

  bool will_write = false;
  for (int lev = 0; lev < Nlev && !will_write; ++lev)
  {
    for (MyParIter pti(*this, lev); pti.isValid() && !will_write; ++pti)
    {
      auto& aos = pti.GetArrayOfStructs();
      auto& soa = pti.GetStructOfArrays();

      for (int pindex=0; pindex<aos.size() && !will_write; ++pindex)
      {
        ParticleType& p = aos[pindex];
        will_write |= (p.id() > 0);
      }
    }
  }

  if (will_write)
  {
    std::string fileName = outFile + "/str_";
    fileName = Concatenate(fileName,myProc) + ".dat";
    std::ofstream ofs(fileName.c_str());
    ofs << "VARIABLES = " << AMREX_D_TERM("X ", "Y ", "Z") << '\n';

    for (int lev = 0; lev < Nlev; ++lev)
    {
      for (MyParIter pti(*this, lev); pti.isValid(); ++pti)
      {
        auto& aos = pti.GetArrayOfStructs();
        auto& soa = pti.GetStructOfArrays();

        for (int pindex=0; pindex<aos.size(); ++pindex)
        {
          ofs << "ZONE I=1 J=" << nPtsOnStrm << " k=1 FORMAT=POINT\n";
          for (int j=0; j<nPtsOnStrm; ++j)
          {
            int offset = j*RealData::ncomp;
            dim3 vals = {AMREX_D_DECL(soa.GetRealData(offset + RealData::xloc)[pindex],
                                      soa.GetRealData(offset + RealData::yloc)[pindex],
                                      soa.GetRealData(offset + RealData::zloc)[pindex])};
            for (int d=0; d<AMREX_SPACEDIM; ++d)
            {
              ofs << vals[d] << " ";
            }
            ofs << '\n';
          }
        }
      }
    }
    ofs.close();
  }
}
