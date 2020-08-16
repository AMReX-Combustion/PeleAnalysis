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

Vec3d normal(const Vec3d &x0, const Vec3d &x1, const Vec3d &x2)
{
    Vec3d x12(x1-x0),x23(x2-x1);
    Vec3d normal;

    normal[0] = x12[1]*x23[2]-x12[2]*x23[1];
    normal[1] = x12[2]*x23[0]-x12[0]*x23[2];
    normal[2] = x12[0]*x23[1]-x12[1]*x23[0];

    float r = std::sqrt(normal[0]*normal[0]+normal[1]*normal[1]+normal[2]*normal[2]);
    normal[0] /= r;
    normal[1] /= r;
    normal[2] /= r;
    return normal;
}

static
std::vector<std::string> parseVarNames(std::istream& is)
{
  std::string line;
  std::getline(is,line);
  return Tokenize(line,std::string(", "));
}

static std::string parseTitle(std::istream& is)
{
  std::string line;
  std::getline(is,line);
  return line;
}

std::string
getFileRoot(const std::string& infile)
{
  vector<std::string> tokens = Tokenize(infile,std::string("/"));
  return tokens[tokens.size()-1];
}

void EB2::TriangulatedIF::buildDistance(/*int   argc,char* argv[]*/const std::string& IsoFile)
{
//  amrex::Initialize(/*argc,argv*/);
//  {
      ParmParse pp;

    // Read in isosurface
      //std::string isoFile; 
    //  pp.get("isoFile",isoFile);
    std::string isoFile = IsoFile;
    
    if (ParallelDescriptor::IOProcessor()) {
      std::cerr << "Reading isoFile... " << isoFile << std::endl;
    }
    
    FArrayBox nodes;
    Vector<int> faceData;
    int nElts;
    int nodesPerElt;
    int nNodes, nCompNodes;
    std::vector<std::string> surfNames;
  
    std::ifstream ifs;
    ifs.open(isoFile.c_str(),std::ios::in|std::ios::binary);

    const std::string title = parseTitle(ifs);
    surfNames = parseVarNames(ifs);
    nCompNodes = surfNames.size();

    ifs >> nElts;
    ifs >> nodesPerElt;

    Print() << "nelts: " << nElts << std::endl;
    Print() << "nodesperelt: " << nodesPerElt << std::endl;

    FArrayBox tnodes;
    tnodes.readFrom(ifs);
    nNodes = tnodes.box().numPts();

    // transpose the data so that the components are 'in the right spot for fab data'
    nodes.resize(tnodes.box(),nCompNodes);
    Real** np = new Real*[nCompNodes];
    for (int j=0; j<nCompNodes; ++j)
      np[j] = nodes.dataPtr(j);

    Real* ndat = tnodes.dataPtr();
    for (int i=0; i<nNodes; ++i)
    {
      for (int j=0; j<nCompNodes; ++j)
      {
        np[j][i] = ndat[j];
      }
      ndat += nCompNodes;
    }
    delete [] np;
    tnodes.clear();

    faceData.resize(nElts*nodesPerElt,0);
    ifs.read((char*)faceData.dataPtr(),sizeof(int)*faceData.size());
    ifs.close();

    Print() << "Read " << nElts << " elements and " << nNodes << " nodes" << std::endl;

    int nCell = 64; 
    pp.query("nCell",nCell);
    int max_grid_size = 1; pp.query("max_grid_size",max_grid_size);
    Box domain(IntVect(D_DECL(0,0,0)),
               //IntVect(D_DECL(nCell-1,nCell-1,nCell-1)));
               IntVect(D_DECL(3-1,3-1,9-1)));
    BoxArray grids(domain);
    grids.maxSize(max_grid_size);
    RealBox probDomain({D_DECL(0.000,0.000,0.0000)},{D_DECL(0.03,0.03,.09)});  
  
    Array<int,AMREX_SPACEDIM> is_periodic = {D_DECL(0,0,0)};
    geom.define(domain,probDomain,0,is_periodic);
    const Real* dx = geom.CellSize();
    const Real* plo = geom.ProbLo();

    distanceMF.define(grids,DistributionMapping(grids),1,0);
    for (MFIter mfi(distanceMF); mfi.isValid(); ++mfi)
    {
      std::vector<Vec3d> pointList;
      std::vector<Vec3ui> triList;
      std::vector<Vec3d> nList;
      Vec3d ni;

      for (int node=0; node<nNodes; ++node) {
        const IntVect iv(D_DECL(node,0,0));
        pointList.push_back(Vec3d(D_DECL(nodes(iv,0),nodes(iv,1),nodes(iv,2))));
      }
      for (int elt=0; elt<nElts; ++elt) {
        int offset = elt * nodesPerElt;
        triList.push_back(Vec3ui(D_DECL(faceData[offset]-1,faceData[offset+1]-1,faceData[offset+2]-1)));
        nList.push_back(normal(pointList[faceData[offset]-1],pointList[faceData[offset+1]-1],pointList[faceData[offset+2]-1]));
      }

      // Copy into member data
      for (int node=0; node<nNodes; ++node) 
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

      const Box& vbox = grids[mfi.index()];
      Vec3d local_origin(plo[0] + vbox.smallEnd()[0]*dx[0],
                         plo[1] + vbox.smallEnd()[1]*dx[1],
                         plo[2] + vbox.smallEnd()[2]*dx[2]);
      Array3d phi_grid;
      double dx1 = double(dx[0]);

      make_level_set3(triList, pointList, nList,local_origin, dx1,
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
