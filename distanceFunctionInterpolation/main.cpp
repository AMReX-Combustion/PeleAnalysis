#include <AMReX_Array.H>
#include <AMReX_MultiFab.H>
#include <AMReX_EB2.H>
#include <AMReX_EB2_IF_Base.H>
#include <vector>
#include <string.h>
#include <AMReX_MultiFab.H>
#include <AMReX_EB2_IF_Triangulated.H>
void main(int argc, char* argv[])
{
   amrex::Initialize(argc,argv);
   {

       char* isoFile="flatplain";

       amrex::EB2::TriangulatedIF Tri(isoFile);

       double A[3];
       A[0]=0;
       A[1]=0;
       A[2]=0;

       double distance=(Tri.distanceInterpolation_->distance(A));

   }
   amrex::Finalize();
}
