#include <string>
#include <iostream>
#include <set>

#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>

#include <AMReX_PlotFileUtil.H>
#include <AMReX_FillPatchUtil.H>

using namespace amrex;
using std::set;
using std::vector;
using std::string;
using std::endl;
using std::cerr;


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

struct Edge
{
  Edge(int aL, int aR) noexcept : mL(aL), mR(aR) {}
  int L() const noexcept {return mL;}
  int R() const noexcept {return mR;}
  Edge reverse() const noexcept {return Edge(mR,mL);}
protected:
  int mL, mR;
};

struct Compare
{
  bool operator() (const Edge& L, const Edge& R) const noexcept {
    int Lmin = std::min(L.L(),L.R());
    int Lmax = std::max(L.L(),L.R());
    int Rmin = std::min(R.L(),R.R());
    int Rmax = std::max(R.L(),R.R());
    if (Lmin == Rmin) {
      return Lmax < Rmax;
    }
    return Lmin < Rmin;
  }
};

int
main (int   argc,
      char* argv[])
{
  amrex::Initialize(argc,argv);
  {
    ParmParse pp;

    // Read in isosurface
    std::string isoFile; pp.get("isoFile",isoFile);
    if (ParallelDescriptor::IOProcessor())
      std::cerr << "Reading isoFile... " << isoFile << std::endl;
    
    Real strt_io = ParallelDescriptor::second();

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

    set<Edge,Compare> edgeSet;
    for (int elt=0; elt<nElts; ++elt) {
      int offset = elt * nodesPerElt;

      Edge e1(faceData[offset+0],faceData[offset+1]);
      std::pair<set<Edge,Compare>::const_iterator,bool> it1 = edgeSet.insert(e1);
      if (!(it1.second)) {
        AMREX_ALWAYS_ASSERT(edgeSet.find(e1.reverse())!=edgeSet.end());
      }
      Edge e2(faceData[offset+1],faceData[offset+2]);
      std::pair<set<Edge,Compare>::const_iterator,bool> it2 = edgeSet.insert(e2);
      if (!(it2.second)) {
        AMREX_ALWAYS_ASSERT(edgeSet.find(e2.reverse())!=edgeSet.end());
      }
      Edge e3(faceData[offset+2],faceData[offset+0]);
      std::pair<set<Edge,Compare>::const_iterator,bool> it3 = edgeSet.insert(e3);
      if (!(it3.second)) {
        AMREX_ALWAYS_ASSERT(edgeSet.find(e3.reverse())!=edgeSet.end());
      }
    }
    Print() << "Found " << edgeSet.size() << " edges (nElts * 3 = " << nElts*3 << ")" << std::endl;
    Print() << "All shared edges are consistently numbered." << std::endl;
  }
  amrex::Finalize();
  return 0;
}
