

struct distFuncInterp
{
  distFuncInterp(const MultiFab& a_distance) : distanct(a_distance) {}
  Real operator() (Real x[AMREX_SPACEDIM]) {
    IntVect iv;
    for (int i=0; i<AMREX_SPACEDIM; ++i) {
      iv[i] = index(x[i],i);
    }
    const auto& ba = distanceMF.boxArray();
    int boxID = -1;
    for (int i=0; i<ba.size() && boxID<0; ++i) {
      if (ba[i].contains(iv)) {
        boxID = i;
      }
    }

    Real result;
    if (boxID < 0) {
      // Get distance from coarser representation
    }
    else
    {
      // Using the FArrayBox in distanceMF[boxID], interpolate between iv and iv+IntVect(1,1,1)
    }
    return result;
  }

  int index(Real loc, int dir) {
    return floor((loc - dlo[dir])/dx[dir]);
  }
protected:
  Real dx[AMREX_SPACEDIM];
  Real dlo[AMREX_SPACEDIM];
  const MultiFab& distanceMF;
};

  
