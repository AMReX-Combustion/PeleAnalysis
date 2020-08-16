#include <AMReX_Array.H>
#include <AMReX_MultiFab.H>
#include <AMReX_EB2.H>
#include <AMReX_EB2_IF_Base.H>
#include <vector>
#include <string.h>
#include <AMReX_MultiFab.H>
#include <AMReX_EB2_IF_Triangulated.H>

using namespace amrex;

int
main (int   argc,
      char* argv[])
{
  Initialize(argc,argv);
  {
    std::string isoFile("flatplain");
    EB2::TriangulatedIF Tri(isoFile);    
    Real A[3] = {0,0,0};
    Real distance=(Tri.distanceInterpolation_->distance(A));
    Print() << "Distance is " << distance << std::endl;
  }
  Finalize();
}
