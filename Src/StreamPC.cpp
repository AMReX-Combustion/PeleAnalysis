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
  
  for (int lev=0; lev<Nlev; ++lev)
  {
    const Geometry& geom = Geom(lev);
    const Real* dx = geom.CellSize();
    const Real* plo = geom.ProbLo();

    for (MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi)
    {
      const Box& tile_box  = mfi.tilebox();
      if (BL_SPACEDIM<3 || tile_box.contains(IntVect(D_DECL(0,50,107)))) {
      const RealBox tile_realbox{tile_box, geom.CellSize(), geom.ProbLo()};
      const int grid_id = mfi.index();
      const int tile_id = mfi.LocalTileIndex();
      auto& particle_tile = GetParticles(lev)[std::make_pair(grid_id,tile_id)];

      for (IntVect iv = tile_box.smallEnd(); iv <= tile_box.bigEnd(); tile_box.next(iv))
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
      }}
    }
  }
}

void
StreamParticleContainer::
SetParticleLocation(int a_streamLoc)
{    
  BL_PROFILE("StreamParticleContainer::SetParticleLocation");
  
  for (int lev = 0; lev < Nlev; ++lev)
  {
    for (MyParIter pti(*this, lev); pti.isValid(); ++pti)
    {
      auto& aos = pti.GetArrayOfStructs();
      auto& soa = pti.GetStructOfArrays();

      for (auto& p : aos)
      {
        auto id = p.id();
        int offset = RealData::ncomp * a_streamLoc;
        AMREX_D_EXPR(p.pos(0) = soa.GetRealData(offset + RealData::xloc)[id],
                     p.pos(1) = soa.GetRealData(offset + RealData::yloc)[id],
                     p.pos(2) = soa.GetRealData(offset + RealData::zloc)[id]);
      }
    }
  }
}
