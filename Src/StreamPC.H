#ifndef STREAM_PARTICLE_CONTAINER_H_
#define STREAM_PARTICLE_CONTAINER_H_

#include <AMReX_Particles.H>

struct RealData
{
  enum {
    AMREX_D_DECL(xloc = 0,
                 yloc,
                 zloc),
    ncomp
  };
  static const int nPointOnStream = 500;
  static const int sizeOfRealStreamData = nPointOnStream * ncomp;
};

class StreamParticleContainer
  : public amrex::ParticleContainer<0, 3, RealData::sizeOfRealStreamData, 0>
{
public:

  using MyParIter = amrex::ParIter<0, 3, RealData::sizeOfRealStreamData, 0>;
    
  StreamParticleContainer (const amrex::Vector<amrex::Geometry>            & a_geoms,
                           const amrex::Vector<amrex::DistributionMapping> & a_dmaps,
                           const amrex::Vector<amrex::BoxArray>            & a_bas,
                           const amrex::Vector<int>                        & a_rrs);

  void InitParticles(const amrex::Vector<amrex::Vector<amrex::Real>>& a_locs);
  
  void SetParticleLocation(int a_streamPos);
  
  void ComputeNextLocation(int                                   a_fromLoc,
                           amrex::Real                           a_delta_t,
                           const amrex::Vector<amrex::MultiFab>& a_vectorField);

  void WriteStreamAsTecplot(const std::string& outfile);

  int NumLevels () const {return Nlev;}

protected:
  int Nlev;

private:
  
};

#endif
