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

Vec3f normal(const Vec3f &x0, const Vec3f &x1, const Vec3f &x2)
{
    Vec3f x12(x1-x0),x23(x2-x1);
    Vec3f normal;

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

void EB2::TriangulatedIF::buildDistance(/*int   argc,char* argv[]*/char* IsoFile)
{
//  amrex::Initialize(/*argc,argv*/);
//  {
      ParmParse pp;

    // Read in isosurface
      //std::string isoFile; 
    //  pp.get("isoFile",isoFile);
    std::string isoFile = IsoFile;
    
    if (ParallelDescriptor::IOProcessor())
      std::cerr << "Reading isoFile... " << isoFile << std::endl;
    
    Real strt_io = ParallelDescriptor::second();

    FArrayBox nodes;
    Vector<long> faceData;
    long nElts;
    long nodesPerElt;
    long nNodes, nCompNodes;
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
    for (long j=0; j<nCompNodes; ++j)
      np[j] = nodes.dataPtr(j);

    Real* ndat = tnodes.dataPtr();
    for (long i=0; i<nNodes; ++i)
    {
      for (long j=0; j<nCompNodes; ++j)
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
   // grids(domain);
     grids.maxSize(max_grid_size);
    //RealBox probDomain({D_DECL(0,0,0)},{D_DECL(1,1,1)});
//    RealBox probDomain({D_DECL(-0.0033,-0.0033,-0.0099)},{D_DECL(0.0033,0.0033,.0099)});
    RealBox probDomain({D_DECL(0.000,0.000,0.0000)},{D_DECL(0.03,0.03,.09)});  
  
    Array<int,AMREX_SPACEDIM> is_periodic = {D_DECL(0,0,0)};
    Geometry geom(domain,probDomain,0,is_periodic);
  //  this->geom(domain,probDomain,0,is_periodic);
    const Real* dx = geom.CellSize();
    const Real* plo = geom.ProbLo();

    MultiFab distance(grids,DistributionMapping(grids),1,0);
  //  this->Dfab(grids,DistributionMapping(grids),1,0);

    for (MFIter mfi(Dfab); mfi.isValid(); ++mfi) {

      // TODO: Prune surface

 /*     std::vector<Vec3f> vertList;
      std::vector<Vec3ui> faceList;
      std::vector<Vec3f> normalList;
      Vec3f ni;
*/
      for (long node=0; node<nNodes; ++node) {
        const IntVect iv(D_DECL(node,0,0));
        vertList.push_back(Vec3f(D_DECL(nodes(iv,0),nodes(iv,1),nodes(iv,2))));
      }
      for (long elt=0; elt<nElts; ++elt) {
        int offset = elt * nodesPerElt;
        faceList.push_back(Vec3ui(D_DECL(faceData[offset]-1,faceData[offset+1]-1,faceData[offset+2]-1)));
        ni = normal(vertList[faceData[offset]-1],vertList[faceData[offset+1]-1],vertList[faceData[offset+2]-1]);
         
         normalList.push_back(ni);
        
   /*       Print()<<"===========================Element No  "<<elt<<std::endl; 
         Print()<<"normal"<<std::endl;*/
      }

      const Box& vbox = grids[mfi.index()];
      Vec3f local_origin(plo[0] + vbox.smallEnd()[0]*dx[0],
                         plo[1] + vbox.smallEnd()[1]*dx[1],
                         plo[2] + vbox.smallEnd()[2]*dx[2]);
      Array3f phi_grid;
      float dx1 = float(dx[0]);

      make_level_set3(faceList, vertList, normalList,local_origin, dx1,
                      vbox.length(0),vbox.length(1),vbox.length(2), phi_grid);

   //   vertList.clear();
   //   faceList.clear();

      const auto& d = Dfab.array(mfi);
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
    //amrex::EB_WriteSingleLevelPlotfile("plt", distance, {"den"}, geom, 0.0, 0);
    //VisMF::Write(distance,"distance");
    
    
  //  WriteSingleLevelPlotfile("Distance.out",distance,{"distance"},geom,0.0,0);

  //  TriangulatedIF::distanceInterpolation(distance);
//  }
//  amrex::Finalize();
//  return 0;
}
