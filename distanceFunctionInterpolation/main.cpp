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

    std::string isoFileName, isoFileType;
    pp.get("isoFileName",isoFileName);
    pp.get("isoFileType",isoFileType);
//    isoFileName="Stldirtdum.stl";
//    isoFileType="stl_binary";

    EB2::TriangulatedIF Tri(isoFileName, isoFileType);     
   
    Print()<<"finish loading"<<std::endl; 
  
    Real timer_init=0.,timer_end=0.;
 
    timer_init = amrex::second();

    Vector<int> nCell(AMREX_SPACEDIM);
    pp.getarr("nCell",nCell,0,AMREX_SPACEDIM);
    Box domain(IntVect(D_DECL(0,0,0)),
               IntVect(D_DECL(nCell[0]-1, nCell[1]-1, nCell[2]-1)));

    Vector<Real> probLo(AMREX_SPACEDIM), probHi(AMREX_SPACEDIM);
    pp.getarr("probLo",probLo,0,AMREX_SPACEDIM);
    pp.getarr("probHi",probHi,0,AMREX_SPACEDIM);
    

    //RealBox probDomain(&(probLo[0]),&(probHi[0]));

    RealBox probDomain({AMREX_D_DECL(probLo[0],probLo[1],probLo[2])}, {AMREX_D_DECL(probHi[0],probHi[1],probHi[2])});


    //std::cout<<"probLo="<<"("<<probLo[0]<<","<<probLo[1]<<","<<probLo[2]<<")"<<std::endl;

    //std::cout<<"probLo="<<"("<<probHi[0]<<","<<probHi[1]<<","<<probHi[2]<<")"<<std::endl;

    Array<int,AMREX_SPACEDIM> is_periodic = {D_DECL(0,0,0)};
    Geometry geom(domain,probDomain,0,is_periodic);
    

    BoxArray grids(domain);
    int max_grid_size = 32; pp.query("max_grid_size",max_grid_size);
    grids.maxSize(max_grid_size);
    DistributionMapping dm = DistributionMapping(grids);

    Print()<<"start finalizing"<<std::endl;

    Tri.finalize(geom,grids,dm);   

    Print()<<"finish finalizing"<<std::endl;

    Vector<Real> pt(AMREX_SPACEDIM);
    pp.getarr("point",pt,0,AMREX_SPACEDIM);

    Print() << "Distance is " << Tri(AMREX_D_DECL(pt[0],pt[1],pt[2])) << std::endl;
    timer_end = amrex::second();

    Print()<<"time cost"<<(timer_end - timer_init)<< std::endl;
  }
  Finalize();
}
