#include <string>
#include <iostream>
#include <set>
#include <cmath>
#include <AMReX.H>
#include <AMReX_EB2.H>
#include <AMReX_EB2_IF.H>
//#include <AMReX_EB2_IF.H>

#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>

#include <AMReX_PlotFileUtil.H>
#include <AMReX_FillPatchUtil.H>
#include "makelevelset3.h"
#include "AMReX_EB2_IF_Triangulated.H"
using namespace amrex;
using std::list;
using std::vector;
using std::string;
using std::endl;
using std::cerr;



void EB2::TriangulatedIF::buildDistance(/*int   argc,char* argv[]const std::string& isoFile*/)
{
//  amrex::Initialize(/*argc,argv*/);
//  {
   //   ParmParse pp;

    // Read in isosurface
      //std::string isoFile; 
    //  pp.get("isoFile",isoFile);
    //std::string isoFile = IsoFile;
    
    
    for (MFIter mfi(distanceMF); mfi.isValid(); ++mfi)
    {
   /*   std::vector<Vec3r> pointList;
	std::vector<Vec3ui> triList;
	std::vector<Vec3r> nList;
      Vec3r ni;
   */
  /*    for (int node=0; node<nNodes; ++node) {
        const IntVect iv(D_DECL(node,0,0));
        vertList.push_back(Vec3r(D_DECL(nodes(iv,0),nodes(iv,1),nodes(iv,2))));
      }
      for int elt=0; elt<nElts; ++elt) {
      int offset = elt * nodesPerElt;
        faceList.push_back(Vec3ui(D_DECL(faceData[offset]-1,faceData[offset+1]-1,faceData[offset+2]-1)));
        normalList.push_back(normal(vertList[faceData[offset]-1],vertList[faceData[offset+1]-1],vertList[faceData[offset+2]-1]));
      }
    std::cout<<"loaded data"<<std::endl;*/
      // Copy into member data
/*      for (int node=0; node<nNodes; ++node) 
      {
        vertList.push_back(std::vector<Real>() );
        for(int j=0;j<3;j++) {
          vertList[node].push_back(pointList[node][j]);
        }
          
      }
      for (int elt=0; elt<nElts; ++elt) 
      {
        faceList.push_back(std::vector<int>() );
        normalList.push_back(std::vector<Real>() );
        for(int j=0;j<3;j++)
        {
            faceList[elt].push_back(triList[elt][j]);
            normalList[elt].push_back(nList[elt][j]);
        }
      }
*/
      const Real* dx = this->geom().CellSize();
      const Real* plo = this->geom().ProbLo();


      const Box& vbox = this->grids()[mfi.index()];
      Vec3r local_origin(plo[0] + vbox.smallEnd()[0]*dx[0],
                         plo[1] + vbox.smallEnd()[1]*dx[1],
                         plo[2] + vbox.smallEnd()[2]*dx[2]);
      Array3r phi_grid;
      double dx1 = double(dx[0]);

      make_level_set3(faceList, vertList, normalList,local_origin, dx1,
                      vbox.length(0),vbox.length(1),vbox.length(2), phi_grid);

      const auto& d = distanceMF.array(mfi);
      const int* lo = vbox.loVect();
      const int* hi = vbox.hiVect();
      for (int k=lo[2]; k<=hi[2]; ++k)
      {
        int kL=k-lo[2];
        for (int j=lo[1]; j<=hi[1]; ++j)
        {
          int jL=j-lo[1];
          for (int i=lo[0]; i<=hi[0]; ++i)
          {
            int iL=i-lo[0];
            d(i,j,k) = phi_grid(iL,jL,kL);
          }
        }
      }
    }
    //amrex::WriteSingleLevelPlotfile("plt", distanceMF, {"den"}, geom, 0.0, 0);
    //VisMF::Write(distance,"distance");    
}
