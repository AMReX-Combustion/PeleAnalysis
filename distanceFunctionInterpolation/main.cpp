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
    ParmParse pp;    
    int nCell = 64; 
    pp.query("nCell",nCell);
    int max_grid_size = 32; pp.query("max_grid_size",max_grid_size);
    Box domain(IntVect(D_DECL(0,0,0)),
               //IntVect(D_DECL(nCell-1,nCell-1,nCell-1)));
               IntVect(D_DECL(64-1,64-1,192-1)));
    BoxArray grids(domain);
    grids.maxSize(max_grid_size);
    RealBox probDomain({D_DECL(0.000,0.000,0.0000)},{D_DECL(0.03,0.03,.09)});  
  
    Array<int,AMREX_SPACEDIM> is_periodic = {D_DECL(0,0,0)};
    Geometry geom(domain,probDomain,0,is_periodic);
 


    std::string isoFile("flatplain");
    std::string typeName("mef");    

    DistributionMapping dm = DistributionMapping(grids);
    EB2::TriangulatedIF Tri(isoFile,typeName);    
    Tri.finalize(geom,grids,dm);   


    Real A[3] = {0,0,0};
    Real distance=(Tri.distanceInterpolation_->distance(A));
    Print() << "Distance is " << distance << std::endl;
  }
  Finalize();
}
