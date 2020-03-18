#include <AMReX_Geometry.H>
#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Utility.H>
#include <AMReX_MultiFab.H>

#include "StreamPC.H"

using namespace amrex;

StreamParticleContainer::
StreamParticleContainer(const Vector<Geometry>            & a_geoms,
                        const Vector<DistributionMapping> & a_dmaps,
                        const Vector<BoxArray>            & a_bas,
                        const Vector<int>                 & a_rrs)
  : ParticleContainer<0, 2, RealData::sizeOfRealStreamData, 0> (a_geoms, a_dmaps, a_bas, a_rrs),
  Nlev(a_geoms.size())
{
}

void
StreamParticleContainer::
InitParticles()
{
  BL_PROFILE("StreamParticleContainer::InitParticles");

  int num_ppc = 2;
  int streamLoc = 0;
  int offset = RealData::ncomp*streamLoc;

  int finestLevel = Nlev - 1;
  std::vector< std::pair<int,Box> > isects;
  FArrayBox mask;
  for (int lev=0; lev<Nlev; ++lev)
  {
    const auto& geom = Geom(lev);
    const auto& dx = geom.CellSize();
    const auto& plo = geom.ProbLo();

    BoxArray baf;
    if (lev < finestLevel) {
      baf = BoxArray(ParticleBoxArray(lev+1)).coarsen(this->GetParGDB()->refRatio(lev));
    }

    for (MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi)
    {
      const Box& tile_box  = mfi.tilebox();
      if (BL_SPACEDIM<3 || tile_box.contains(IntVect(D_DECL(0,50,107)))) {
      const RealBox tile_realbox{tile_box, geom.CellSize(), geom.ProbLo()};
      const int grid_id = mfi.index();
      const int tile_id = mfi.LocalTileIndex();
      auto& particle_tile = GetParticles(lev)[std::make_pair(grid_id,tile_id)];

      mask.resize(tile_box,1);
      mask.setVal(1);
      if (lev < finestLevel) {
        isects = baf.intersections(tile_box);
        for (const auto& p : isects) {
          mask.setVal(0,p.second,0,1);
        }
      }
      mask.setVal(0);

      for (IntVect iv = tile_box.smallEnd(); iv <= tile_box.bigEnd(); tile_box.next(iv))
      {
        Array<Real,BL_SPACEDIM> tloc = {AMREX_D_DECL(plo[0] + (iv[0] + 0.5)*dx[0],
                                                     plo[1] + (iv[1] + 0.5)*dx[1],
                                                     plo[2] + (iv[2] + 0.5)*dx[2])};
        if (lev==finestLevel && tloc[1] > .003 && tloc[1] < .005) mask(iv,0) = 1;
        
        if (mask(iv,0) > 0)
        {
          for (int i_part=0; i_part<num_ppc; i_part++)
          {
            Array<Real,BL_SPACEDIM> loc = {AMREX_D_DECL(plo[0] + (iv[0] + 0.5)*dx[0],
                                                        plo[1] + (iv[1] + 0.5)*dx[1],
                                                        plo[2] + (iv[2] + 0.5)*dx[2])};

            ParticleType p;
            p.id()  = ParticleType::NextID();
            p.cpu() = ParallelDescriptor::MyProc();

            AMREX_D_EXPR(p.pos(0) = loc[0],
                         p.pos(1) = loc[1],
                         p.pos(2) = loc[2]);

            p.idata(0) = streamLoc;
            p.idata(1) = i_part==0 ? +1 : -1;

            std::array<double, RealData::sizeOfRealStreamData> real_attribs;
            for (int i=0; i<real_attribs.size(); ++i) real_attribs[i] = 0;
          
            AMREX_D_EXPR(real_attribs[offset + RealData::xloc] = loc[0],
                         real_attribs[offset + RealData::yloc] = loc[1],
                         real_attribs[offset + RealData::zloc] = loc[2]);

            AMREX_ASSERT(this->Index(p, lev) == iv);

            particle_tile.push_back(p);
            particle_tile.push_back_real(real_attribs);
          }
        }
      }}
    }
  }
}

void
StreamParticleContainer::
SetParticleLocation(int a_streamLoc)
{    
  BL_PROFILE("StreamParticleContainer::SetParticleLocation");
  
  int offset = RealData::ncomp * a_streamLoc;
  for (int lev = 0; lev < Nlev; ++lev)
  {
    for (MyParIter pti(*this, lev); pti.isValid(); ++pti)
    {
      auto& aos = pti.GetArrayOfStructs();
      auto& soa = pti.GetStructOfArrays();

      for (int pindex=0; pindex<aos.size(); ++pindex)
      {
        ParticleType& p = aos[pindex];
        if (p.id() > 0)
        {
          AMREX_D_EXPR(p.pos(0) = soa.GetRealData(offset + RealData::xloc)[pindex],
                       p.pos(1) = soa.GetRealData(offset + RealData::yloc)[pindex],
                       p.pos(2) = soa.GetRealData(offset + RealData::zloc)[pindex]);
        }
      }
    }
  }
}

typedef Array<Real,AMREX_SPACEDIM> dim3;
typedef Array<int,AMREX_SPACEDIM> int3;

static bool IsOK(const dim3& x, const Real* plo, const Real* phi)
{
  bool ok = true;
  for (int i=0; i<AMREX_SPACEDIM; ++i) {
    ok |= x[i] < plo[i] || x[i] > phi[i];
  }
  return ok;
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

  if (!IsOK(x,plo,phi)) return false;

  for (int d=0; d<AMREX_SPACEDIM; ++d) {
    b[d] = int( (x[d] - plo[d]) / dx[d] - 0.5 );
    n[d] = ( x[d] - ( (b[d] + 0.5 ) * dx[d] + plo[d] ) )/dx[d];
    n[d] = std::max(0., std::min(1.,n[d]));
  }

  const auto& gbx = gfab.box();
  const auto& glo = gbx.smallEnd();
  const auto& ghi = gbx.bigEnd();
  for (int d=0; d<AMREX_SPACEDIM; ++d) {
    if (b[d] < glo[d] ||  b[d] > ghi[d]-1) {
      Print() << "d,b,glo,ghi: " << d << " " << b[d] << " " << glo[d] << " " << ghi[d] << std::endl;
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

static void
RK4(dim3 & x,Real dt,const FArrayBox& v,const Real* dx,const Real* plo,const Real* phi,int dir, bool& ok)
{
  dim3 vec, k1, k2, k3, k4;
  dim3 xx = x;
  ok = ntrpv(xx,v,dx,plo,phi,vec);
  if ( !ok ) return;
  vnrml(vec,dir);

  for (int d=0; d<AMREX_SPACEDIM; ++d) {
    k1[d] = vec[d] * dt;
    xx[d] = x[d] + k1[d] * 0.5;
  }
  ok = ntrpv(xx,v,dx,plo,phi,vec);
  if ( !ok ) return;
  vnrml(vec,dir);

  for (int d=0; d<AMREX_SPACEDIM; ++d) {
    k2[d] = vec[d] * dt;
    xx[d] = x[d] + k2[d] * 0.5;
  }
  ok = ntrpv(xx,v,dx,plo,phi,vec);
  if ( !ok ) return;
  vnrml(vec,dir);

  for (int d=0; d<AMREX_SPACEDIM; ++d) {
    k3[d] = vec[d] * dt;
    xx[d] = x[d] + k3[d];
  }
  ok = ntrpv(xx,v,dx,plo,phi,vec);
  if ( !ok ) return;
  vnrml(vec,dir);

  const Real third = 1./3.;
  const Real sixth = 1./6.;
  for (int d=0; d<AMREX_SPACEDIM; ++d) {
    k4[d] = vec[d] * dt;
    x[d] += (k1[d] + k4[d])*sixth + (k2[d] + k3[d])*third;
  }
}

void
StreamParticleContainer::
ComputeNextLocation(int                      a_fromLoc,
                    Real                     a_delta_t,
                    const Vector<MultiFab> & a_vectorField)
{    
  BL_PROFILE("StreamParticleContainer::ComputeNextLocation");

  SetParticleLocation(a_fromLoc);
  
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
        if (p.id() > 0)
        {
          const int dir = p.idata(1);
          dim3 x = {AMREX_D_DECL(p.pos(0), p.pos(1), p.pos(2))};
          bool ok = true;
          RK4(x,a_delta_t,v,dx,plo,phi,dir,ok);
          if (!ok) {
            Print() << "pindex, pid " << pindex << " " << p.id() << std::endl; 
            Abort("bad RK");
          }
          AMREX_D_EXPR(soa.GetRealData(offset + RealData::xloc)[pindex] = x[0],
                       soa.GetRealData(offset + RealData::yloc)[pindex] = x[1],
                       soa.GetRealData(offset + RealData::zloc)[pindex] = x[2]);
        }
      }
    }
  }
}

void
StreamParticleContainer::
WriteStreamAsTecplot(const std::string& outFile)
{
  // Create a folder and have each processor write their own data, one file per streamline
  auto myProc = ParallelDescriptor::MyProc();
  auto nProcs = ParallelDescriptor::NProcs();
  int cnt = 0;
  
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
          ParticleType& p = aos[pindex];
          if (p.id() > 0)
          {
            ofs << "ZONE I=1 J=" << RealData::nPointOnStream << " k=1 FORMAT=POINT\n";
            for (int j=0; j<RealData::nPointOnStream; ++j)
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
    }
    ofs.close();
  }
}
