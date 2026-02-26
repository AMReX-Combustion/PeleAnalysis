#ifndef STREAM_PARTICLE_CONTAINER_H_
#define STREAM_PARTICLE_CONTAINER_H_

#include <AMReX_Particles.H>

typedef amrex::Array<amrex::Real, AMREX_SPACEDIM> dim3;
typedef amrex::Array<int, AMREX_SPACEDIM> int3;

struct RealData
{
  enum {
    AMREX_D_DECL(xloc = 0, yloc, zloc),
  };
};

class StreamParticleContainer : public amrex::ParticleContainer<0, 3, 0, 0>
{
public:
  using MyParIter = amrex::ParIter<0, 3, 0, 0>;

  StreamParticleContainer(
    int a_nPtsOnStrm,
    const amrex::Vector<amrex::Geometry>& a_geoms,
    const amrex::Vector<amrex::DistributionMapping>& a_dmaps,
    const amrex::Vector<amrex::BoxArray>& a_bas,
    const amrex::Vector<int>& a_rrs,
    const int a_nVar,
    const amrex::IntVect a_is_per,
    const amrex::Vector<std::string>& a_outVarNames);

  void InitParticles(const amrex::Vector<amrex::Vector<amrex::Real>>& a_locs);

  void SetParticleLocation(int a_streamPos);

  void ComputeNextLocation(
    int a_fromLoc,
    const amrex::Real a_hRK,
    const amrex::Vector<amrex::MultiFab>& a_vectorField,
    const int a_cSpace);

  void InterpDataAtLocation(
    int a_fromLoc, const amrex::Vector<amrex::MultiFab>& a_vectorField);

  void
  WriteStreamAsBinary(const std::string& outfile, amrex::Vector<int>& faceData);

  void WriteStreamAsTecplot(const std::string& outfile);

  bool RK4(
    dim3& x,
    amrex::Real hrk,
    const amrex::FArrayBox& v,
    const amrex::Real* dx,
    const amrex::Real* plo,
    const amrex::Real* phi,
    int dir,
    int cSpace,
    amrex::Real dxFine);

#if AMREX_DEBUG // AJA inspect function
  void InspectParticles(const int nStreamPairs);
#endif

protected:
  int Nlev;
  int nPtsOnStrm;
  int sizeOfRealStreamData;
  int pcomp;
  int fcomp;
  amrex::IntVect is_per;
  amrex::Vector<std::string> outVarNames;

private:
};

#endif
