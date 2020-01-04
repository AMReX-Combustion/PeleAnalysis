
#include <iostream>
#include <list>
#include <vector>
#include <map>
#include <set>
#include <fstream>
#include <iterator>
#include <iomanip>

using std::cout;
using std::cerr;
using std::endl;
using std::list;
using std::vector;
using std::set;
using std::map;
using std::pair;
using std::string;
using std::istream;
using std::ostream;
using std::ofstream;
using std::ifstream;

#include "AMReX_ParmParse.H"
#include "AMReX_ParallelDescriptor.H"
#include "AMReX_Utility.H"
#include "AMReX_FArrayBox.H"

using namespace amrex;

typedef Vector<Real> Point;
typedef list<Point> PointList;
typedef Vector<Point> PointVec;
typedef pair<int,int> Edge;
typedef map<Edge, pair<Point,int> > PMap;
typedef PMap::iterator PMapIt;
static PMap vertCache;

static Real epsilon_DEF = 1.e-8;

static string parseTitle(istream& is);
static vector<string> parseVarNames(istream& is);
static string rootName(const string in);

struct Segment
{
  Segment() : p(2), mLength(-1) {}
  Real Length ();
  const PMapIt& operator[] (int n) const { return p[n]; }
  PMapIt& operator[] (int n) { return p[n]; }
  int ID_l () const {return (*p[0]).second.second;}
  int ID_r () const {return (*p[1]).second.second;}
  void flip ();
  static int xComp,yComp,zComp;
  Vector<PMapIt> p;
    
private:
  void my_length();
  Real mLength;
};

typedef list<Segment> SegList;

bool operator<(const Edge& lhs, const Edge& rhs)
{
    if (std::min(lhs.first,lhs.second) != std::min(rhs.first,rhs.second))
    {
        return (std::min(lhs.first,lhs.second) < std::min(rhs.first,rhs.second));
    }
    return std::max(lhs.first,lhs.second) < std::max(rhs.first,rhs.second);
}

struct EdgeLess
{
    bool operator()(const Edge& lhs, const Edge& rhs) const
        { return operator<(lhs,rhs); }
};

bool operator==(const Edge& lhs, const Edge& rhs)
{
    return (lhs.first == rhs.first && lhs.second == rhs.second) ||
        (lhs.first == rhs.second && lhs.second == rhs.first);
}

bool operator!=(const Edge& lhs, const Edge& rhs)
{
    return !operator==(lhs,rhs);
}


bool operator==(const Segment& lhs, const Segment& rhs)
{
    return lhs.ID_l() == rhs.ID_l() && lhs.ID_r() == rhs.ID_r();
}

bool operator!=(const Segment& lhs, const Segment& rhs)
{
    return !operator==(lhs,rhs);
}

PMapIt VertexInterp(Real isoVal,int isoComp,const Vector<Point>& pts,int p1,int p2);

Vector<Segment> Segmentise(const Vector<Point>&  pts,
                          const Vector<int>&    elt,
                          Real                  isoVal,
                          int                   isoComp);


SegList::iterator FindMySeg(SegList& segs, int idx, Vector<PMapIt>& vertVec)
{
    const PMapIt& vertToFind = vertVec[idx];
    for (SegList::iterator it=segs.begin(); it!=segs.end(); ++it)
    {
        if ( ((*it)[0] == vertToFind) || ((*it)[1] == vertToFind) )
            return it;
    }
    return segs.end();
}


int main (int argc,
	  char* argv[])
{
  Initialize(argc,argv);
  {
    ParmParse pp;

    string infile; pp.get("infile",infile);

    int isoComp; pp.get("isoComp",isoComp);
    Real isoVal; pp.get("isoVal",isoVal);

    ifstream is(infile.c_str());

    auto title = parseTitle(is);
    auto names = parseVarNames(is);
    int nComp = names.size();

    int nElts;
    int MYLEN;
    is >> nElts;
    is >> MYLEN;

    FArrayBox nodeFab;
    nodeFab.readFrom(is);
    AMREX_ALWAYS_ASSERT(nComp==nodeFab.nComp());
    Real* nodeData = nodeFab.dataPtr();
    const int nNodes = nodeFab.box().length(0);
    Vector<Point> GridPts(nNodes);
    for (int i=0; i<nNodes; ++i) {
      GridPts[i].resize(nComp);
      for (int j=0; j<nComp; ++j) {
        GridPts[i][j] = nodeData[i*nComp+j];
      }
    }
    nodeFab.clear(); // Reclaim uneeded space
    cerr << nNodes << " nodes read in with "
         << nComp << " states per node" << endl;

    Vector<int> connData(nElts*MYLEN,0);
    is.read((char*)connData.dataPtr(),sizeof(int)*connData.size());
    Vector<Vector<int> > GridElts(nElts);

    for (int i=0; i<nElts; ++i)
    {
      GridElts[i].resize(MYLEN);
      for (int j=0; j<MYLEN; ++j)
        GridElts[i][j] = connData[i*MYLEN+j] - 1;
    }
    connData.clear(); // Reclaim uneeded space
    cerr << nElts << " elements read in with "
         << MYLEN << " nodes per element" << endl;

    Real Length = 0;
    SegList segments;

    for (int i=0; i<nElts; ++i) {

      Vector<Segment> eltSegs = Segmentise(GridPts,GridElts[i],isoVal,isoComp);
      if (eltSegs.size() > 0) {
    

        for (int j=0; j<eltSegs.size(); ++j)
        {
          Length += eltSegs[j].Length();
          segments.push_back(eltSegs[j]);
        }
      }
    }

    int nSegments = segments.size();
    cout << "Found " << nSegments << " segments "  << endl;

    // Number the isosurface vertices
    int cnt=0;
    for (PMapIt it=vertCache.begin(); it!=vertCache.end(); ++it)
      (*it).second.second = cnt++;

    // Build a reverse Vector into the vertCache
    Vector<PMapIt> vertVec(vertCache.size());
    cnt=0;
    for (PMapIt it=vertCache.begin(); it!=vertCache.end(); ++it)
      vertVec[cnt++] = it;

    // Find a segment with the specified vertex as one of its endpoints,
    // then assemble the list of segments to form the contour line.  If
    // we finish, and segments remain, start a new line.
    SegList segList = segments;
    list<SegList> cLines;
    if (segList.size()>0)
    {
      int idx = segList.front().ID_l();
      int newIdx;
      cLines.push_back(SegList());
        
      while (segList.begin() != segList.end())
      {
        SegList::iterator segIt = FindMySeg(segList,idx,vertVec);
        if (segIt != segList.end())
        {
          int idx_l = (*segIt).ID_l();
          int idx_r = (*segIt).ID_r();
          if ( idx_l == idx )
          {
            newIdx = idx_r;
            cLines.back().push_back(*segIt);
          }
          else
          {
            newIdx = idx_l;
            Segment newSeg = Segment(*segIt); newSeg.flip();
            cLines.back().push_back(newSeg);
          }
                
          segList.erase(segIt);
        }
        else
        {
          cLines.push_back(SegList());
          newIdx = segList.front().ID_l();
        }
            
        idx = newIdx;
      }
        
      // Connect up the line fragments as much as possible
      bool changed;
      do
      {
        changed = false;
        for (std::list<SegList>::iterator it = cLines.begin(); it!=cLines.end(); ++it)
        {
          if (!(*it).empty())
          {
            const int idx_l = (*it).front().ID_l();
            const int idx_r = (*it).back().ID_r();
            for (std::list<SegList>::iterator it1 = cLines.begin(); it1!=cLines.end(); ++it1)
            {
              if (!(*it1).empty() && (*it).front()!=(*it1).front())
              {
                if (idx_r == (*it1).front().ID_l())
                {
                  (*it).splice((*it).end(),*it1);
                  changed = true;
                }
                else if (idx_r == (*it1).back().ID_r())
                {
                  (*it1).reverse();
                  for (SegList::iterator it2=(*it1).begin(); it2!=(*it1).end(); ++it2)
                    (*it2).flip();
                  (*it).splice((*it).end(),*it1);
                  changed = true;
                }
                else if (idx_l == (*it1).front().ID_l())
                {
                  (*it1).reverse();
                  for (SegList::iterator it2=(*it1).begin(); it2!=(*it1).end(); ++it2)
                    (*it2).flip();
                  (*it).splice((*it).begin(),*it1);
                  changed = true;
                }
              }
            }
          }
        }
      } while(changed);
    }
    
    for (std::list<SegList>::iterator it = cLines.begin(); it!=cLines.end();)
    {
      if ((*it).empty())
        cLines.erase(it++);
      else
        it++;
    }
    cerr << "  number of contours " << cLines.size() << endl;

    ofstream os("out.dat");
    os << "VARIABLES =";
    for (int i=0; i<names.size(); ++i) {
      os << " " << names[i];
    }
    os << '\n';

    // One zone per contour line
    for (std::list<SegList>::iterator it = cLines.begin(); it!=cLines.end(); ++it) {
      int nSegNodes = it->size() + 1;
      os << "ZONE ZONETYPE=FELINESEG DATAPACKING=POINT N=" << nSegNodes << " E=" << it->size() << "\n";
      for (SegList::iterator it2=it->begin(); it2!=it->end(); ++it2) {
	const Point& p0 = (*it2)[0]->second.first;
        for (int n=0 ; n <nComp ; ++n )
          os << p0[n] << " ";
        os << "\n" ;
      }
      const auto& p1 = (*it).back()[1]->second.first;
      for (int n=0 ; n <nComp ; ++n )
        os << p1[n] << " ";
      os << "\n" ;

      int cnt=1;
      for (SegList::iterator it2=it->begin(); it2!=it->end(); ++it2) {
        os << cnt++ << " " << cnt << '\n'; // local numbering
      }
    }
    os.close();
  }
  Finalize();
  return 0;
}

static
vector<string> parseVarNames(istream& is)
{
    string line;
    std::getline(is,line);
    return Tokenize(line,string(", "));
}

static string parseTitle(istream& is)
{
    string line;
    std::getline(is,line);
    return line;
}

static string rootName(const string inStr)
{
    const string dirSep("/");
    const vector<string>& res  = Tokenize(inStr,dirSep);
    const vector<string>& nres = Tokenize(res[res.size()-1],string("."));
    string result = nres[0];
    for (int i=1; i<nres.size()-1; ++i)
        result = result + string(".") + nres[i];
    return result;
}



/*
   Given a grid cell and an isoVal, calculate the line segments
   required to represent the contour through the cell.
   Return an array of (at most 2) line segments
*/
Vector<Segment> Segmentise(const Vector<Point>&  pts,
                          const Vector<int>&    elt,
                          Real                  isoVal,
                          int                   isoComp)
{
   BL_ASSERT(elt.size()==3);
   Vector<PMapIt> vertlist(3);

   const int p0 = elt[0];
   const int p1 = elt[1];
   const int p2 = elt[2];

   bool lo_0 = pts[p0][isoComp] < isoVal;
   bool lo_1 = pts[p1][isoComp] < isoVal;
   bool lo_2 = pts[p2][isoComp] < isoVal;

   int count = 0;
   if (lo_0 ^ lo_1) {
     vertlist[count++] = VertexInterp(isoVal,isoComp,pts,p0,p1);
   }
   if (lo_1 ^ lo_2) {
     vertlist[count++] = VertexInterp(isoVal,isoComp,pts,p1,p2);
   }
   if (lo_2 ^ lo_0) {
     vertlist[count++] = VertexInterp(isoVal,isoComp,pts,p2,p0);
   }

   Vector<Segment> segments(0);
   if (count > 0) {
     BL_ASSERT(count == 2);
     segments.resize(1);
     segments[0][0] = vertlist[0];
     segments[0][1] = vertlist[1];
   }

   return segments;
}

/*
   Linearly interpolate the position where an isosurface cuts
   an edge between two vertices, each with their own scalar value
*/
Point VI_doIt(Real isoVal,int isoComp,const Vector<Point>& pts,int p1,int p2)
{
    BL_ASSERT(p1<pts.size() && p2<pts.size());
    const Point& pt1 = pts[p1];
    const Point& pt2 = pts[p2];

    BL_ASSERT(isoComp!=pts.size() && isoComp<pt1.size() && pt1.size()==pt2.size());

    const Real valp1 = pt1[isoComp];
    const Real valp2 = pt2[isoComp];
    
    if (std::abs(isoVal-valp1) < epsilon_DEF)
        return pt1;
    if (std::abs(isoVal-valp2) < epsilon_DEF)
        return pt2;
    if (std::abs(valp1-valp2) < epsilon_DEF)
        return pt1;
    
    Point res(pt1.size()+2);
//      Point res(pt1.size());
    const Real mu = (isoVal - valp1) / (valp2 - valp1);
    for (int j=0; j<res.size()-2; ++j)
        res[j] = pt1[j] + mu * (pt2[j] - pt1[j]);
    res[pt1.size()] = std::pow(res[0]*res[0] +res[1]*res[1],0.5);
    Real tmp = atan(std::abs(res[1]/(res[0]+1e-5)))*180/3.14;
    res[pt1.size()+1] = tmp;
    if (res[0]<0.0 && res[1] >0.0) res[pt1.size()+1] = 90+tmp;
    if (res[0]<0.0 && res[1] <0.0) res[pt1.size()+1] = 180+tmp;
    if (res[0]>0.0 && res[1] <0.0) res[pt1.size()+1] = 270+tmp;
    return res;
}

PMapIt VertexInterp(Real isoVal,int isoComp,const Vector<Point>& pts,int p1,int p2)
{
    Edge ppair(p1,p2);
    PMapIt fwd,rev;
    fwd = vertCache.find(ppair);
    if (fwd == vertCache.end())
    {
        rev = vertCache.find(Edge(p2,p1));
        
        if (rev == vertCache.end())
        {
            return vertCache.insert(
                std::make_pair(ppair,
                std::make_pair(VI_doIt(isoVal,isoComp,pts,p1,p2),0))).first;
        }
        else
        {
            return rev;
        }
    }
    return fwd;
}

int Segment::xComp = 0;
int Segment::yComp = 1;
int Segment::zComp = 2;

Real
Segment::Length()
{
    if (mLength<0)
        my_length();
    return mLength;
}

void
Segment::my_length()
{
  const Point& p0 = (*p[0]).second.first;
  const Point& p1 = (*p[1]).second.first;
  mLength = std::sqrt(((p1[xComp] - p0[xComp])*(p1[xComp] - p0[xComp]))
		      +((p1[yComp] - p0[yComp])*(p1[yComp] - p0[yComp]))
		      +((p1[zComp] - p0[zComp])*(p1[zComp] - p0[zComp])));
}    

void
Segment::flip()
{
    PMapIt ptmp = p[0];
    p[0] = p[1];
    p[1] = ptmp;
}    

