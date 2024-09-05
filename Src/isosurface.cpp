#include <string>
#include <iostream>
#include <set>
#include <list>

#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>

#include <AMReX_BLFort.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_FillPatchUtil.H>
#include <AMReX_PlotFileUtil.H>
#include "makelevelset3.h"

using namespace amrex;
using std::list;
using std::vector;
using std::string;
using std::endl;
using std::cerr;

static Real isoVal_DEF = 1090.;
static string isoCompName_DEF = "temp";
static Real epsilon_DEF = 1.e-15;

static
void
print_usage (int,
             char* argv[])
{
  std::cerr << "usage:\n";
  std::cerr << argv[0] << " inputs infile=<s> isoCompName=<s> isoVal=<v> [options] \n\tOptions:\n";
  std::cerr << "\t     infile=<s> where <s> is a pltfile\n";
  std::cerr << "\t     isoCompName=<s> where <s> is the quantity being contoured\n";
  std::cerr << "\t     isoVal=<v> where <v> is an isopleth value\n";
  std::cerr << "\t     Choosing quantities to interp to surface: \n";
  std::cerr << "\t       comps=int comp list [overrides sComp/nComp]\n";
  std::cerr << "\t       sComp=start comp[DEF->0]\n";
  std::cerr << "\t       nComp=number of comps[DEF->all]\n";
  std::cerr << "\t     finestLevel=<n> finest level to use in pltfile[DEF->all]\n";
  std::cerr << "\t     writeSurf=<1,0> output surface [DEF->1]\n";
  std::cerr << "\t     surfFormat=<MEF,XDMF> output surface format [DEF->MEF]\n";
  std::cerr << "\t     outfile_base=<s> base name of output file [DEF->gen'd]\n";
  std::cerr << "\t     build_distance_function=<t,f> create cc signed distance function [DEF->f]\n";
  std::cerr << "\t     rm_external_elements=<t,f> remove elts beyond what is needed for watertight surface [DEF->t]\n";
exit(1);
}

// A struct defining an edge as two IntVects (left & right)
struct Edge
{
  // Constructor automatically catching which is left & right
  Edge(const IntVect& lhs, const IntVect& rhs) {
    if (lhs < rhs) {
      IV_l = lhs;
      IV_r = rhs;
    } else {
      IV_l = rhs;
      IV_r = lhs;
    }
  }

  // Relational operator ==
  bool operator== (const Edge& rhs) const
    {
      return IV_l==rhs.IV_l && IV_r==rhs.IV_r;
    }

  // Relational operator <
  // Comparing first the left IntVect
  bool operator< (const Edge& rhs) const
    {
      if (IV_l==rhs.IV_l) {
        return IV_r < rhs.IV_r;
      } else {
        return IV_l < rhs.IV_l;
      }
    }

  IntVect IV_l;
  IntVect IV_r;
};

// Point in space
typedef Vector<Real> Point;

// Map of points associated to an edge
typedef std::map<Edge, Point> PMap;

typedef PMap::iterator PMapIt;

struct PMapItCompare
{
  bool operator() (const PMapIt& lhs, const PMapIt& rhs) const {
    const Point& l = lhs->second;
    const Point& r = rhs->second;
    const void* vl=&l;
    const void* vr=&r;
    return vl<vr;
  }
};

// Define the basic 2D element
struct Segment
{
  Segment() : p(2), mLength(-1) {}
  Real Length ();
  const Vector<Real>& normalVector ();
  const PMapIt& operator[] (int n) const { return p[n]; }
  PMapIt& operator[] (int n) { return p[n]; }
  void flip ();
  static int xComp,yComp;
  Vector<PMapIt> p;

private:
  void my_length();
  Real mLength;
  void my_normVector();
  Vector<Real> mnormVector;
};

int Segment::xComp = 0;
int Segment::yComp = 1;

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
#ifndef NDEBUG
  AMREX_ASSERT(xComp>=0 && yComp>=0);
  for (int i=0; i<p.size(); ++i)
  {
    AMREX_ASSERT(xComp<(*p[i]).second.size());
    AMREX_ASSERT(yComp<(*p[i]).second.size());
  }
#endif
  const Point& p0 = (*p[0]).second;
  const Point& p1 = (*p[1]).second;
  mLength = std::sqrt(((p1[xComp] - p0[xComp])*(p1[xComp] - p0[xComp]))
                      +((p1[yComp] - p0[yComp])*(p1[yComp] - p0[yComp])));
}

void
Segment::my_normVector()
{
    const Point& p0 = (*p[0]).second;
    const Point& p1 = (*p[1]).second;
    mnormVector[0] = (p0[yComp]-p1[yComp])/Length();
    mnormVector[1] = (p1[xComp]-p0[xComp])/Length();
}
const Vector<Real>&
Segment::normalVector()
{
    my_normVector();
    return mnormVector;
}

void
Segment::flip()
{
  PMapIt ptmp = p[0];
  p[0] = p[1];
  p[1] = ptmp;
}

std::ostream&
operator<< (std::ostream& os, const Segment& seg)
{
  const Point& p0 = seg.p[0]->second;
  const Point& p1 = seg.p[1]->second;

  os << '[' << p0[0];
  for (int i = 1; i < AMREX_SPACEDIM; i++)
  {
    os << ", " << p0[i];
  }
  os << "] ... ";

  os << '[' << p1[0];
  for (int i = 1; i < AMREX_SPACEDIM; i++)
  {
    os << ", " << p1[i];
  }
  os << "]";

  return os;
}

typedef list<Segment> SegList;

// Define the basic 3D element
struct Triangle
{
  Triangle() : p(3), mArea(-1) {}
  Real Area();
  const PMapIt& operator[] (int n) const { return p[n]; }
  PMapIt& operator[] (int n) { return p[n]; }
  static int xComp,yComp,zComp;
  Vector<PMapIt> p;

private:
  void my_area();
  Real mArea;
};

int Triangle::xComp = 0;
int Triangle::yComp = 1;
int Triangle::zComp = 2;

Real
Triangle::Area()
{
  if (mArea<0)
    my_area();
  return mArea;
}

void
Triangle::my_area()
{
#ifndef NDEBUG
  for (int i=0; i<p.size(); ++i)
  {
    AMREX_ASSERT(xComp>=0 && xComp<(*p[i]).second.size());
    AMREX_ASSERT(yComp>=0 && yComp<(*p[i]).second.size());
    AMREX_ASSERT(zComp>=0 && zComp<(*p[i]).second.size());
  }
#endif

  const Point& p0 = (*p[0]).second;
  const Point& p1 = (*p[1]).second;
  const Point& p2 = (*p[2]).second;
  mArea = 0.5*sqrt(
    pow(  ( p1[yComp] - p0[yComp])*(p2[zComp]-p0[zComp])
          -(p1[zComp] - p0[zComp])*(p2[yComp]-p0[yComp]), 2)

    + pow(( p1[zComp] - p0[zComp])*(p2[xComp]-p0[xComp])
          -(p1[xComp] - p0[xComp])*(p2[zComp]-p0[zComp]), 2)

    + pow(( p1[xComp] - p0[xComp])*(p2[yComp]-p0[yComp])
          -(p1[yComp] - p0[yComp])*(p2[xComp]-p0[xComp]), 2));
}

typedef list<Triangle> TriList;

/*
  Linearly interpolate the position where an isosurface cuts
  an edge between two vertices, each with their own scalar value
*/
Point
VI_doIt(Real isoVal, int isoComp, const Point& p1, const Point& p2)
{
  AMREX_ASSERT(p1.size()>isoComp && p2.size()>isoComp);

  const Real valp1 = p1[isoComp];
  const Real valp2 = p2[isoComp];

  if (std::abs(isoVal-valp1) < epsilon_DEF)
    return p1;
  if (std::abs(isoVal-valp2) < epsilon_DEF)
    return p2;
  if (std::abs(valp1-valp2) < epsilon_DEF)
    return p1;

  Point res(p1.size());
  const Real mu = (isoVal - valp1) / (valp2 - valp1);
  for (int j=0; j<res.size(); ++j)
    res[j] = p1[j] + mu * (p2[j] - p1[j]);

  return res;
}

PMapIt VertexInterp(Real isoVal, int isoComp,
                    const IntVect& p1,const Point& p1d,
                    const IntVect& p2,const Point& p2d,
                    PMap& vertCache)
{
  PMapIt fwd,rev;
  Edge edge(p1,p2);
  fwd = vertCache.find(edge);
  if (fwd == vertCache.end()) {
    rev = vertCache.find(Edge(p2,p1));

    if (rev == vertCache.end()) {
      std::pair<Edge,Point> ent(edge,VI_doIt(isoVal,isoComp,p1d,p2d));
      std::pair<PMapIt,bool> it = vertCache.insert(ent);
      AMREX_ASSERT(it.second);
      return it.first;
    } else {
      return rev;
    }
  }
  return fwd;
}

#if AMREX_SPACEDIM==2
/*
  Given a grid cell and an isoVal, calculate the line segments
  required to represent the contour through the cell.
  Return an array of (at most 2) line segments
*/

Vector<Segment> Segmentise(const FArrayBox& pts,
                           const FArrayBox& mask,
                           PMap&            vertCache,
                           const IntVect&   baseIV,
                           Real             isoVal,
                           int              isoComp)
{
  Vector<PMapIt> vertlist(4);
  Vector<Segment> segments;

  const IntVect& p0 = baseIV;
  const IntVect  p1 = p0 + amrex::BASISV(0);
  const IntVect  p2 = p1 + amrex::BASISV(1);
  const IntVect  p3 = p0 + amrex::BASISV(1);
  //
  // Bail if any of the point are masked out.
  //
  if (mask(p0)<0 || mask(p1)<0 || mask(p2)<0 || mask(p3)<0)
    return segments;

  const int nComp = pts.nComp();

  Point p0d(nComp); pts.getVal(p0d.dataPtr(),p0);
  Point p1d(nComp); pts.getVal(p1d.dataPtr(),p1);
  Point p2d(nComp); pts.getVal(p2d.dataPtr(),p2);
  Point p3d(nComp); pts.getVal(p3d.dataPtr(),p3);

  int segCase = 0;

  if (p0d[isoComp] < isoVal) segCase |= 1;
  if (p1d[isoComp] < isoVal) segCase |= 2;
  if (p2d[isoComp] < isoVal) segCase |= 4;
  if (p3d[isoComp] < isoVal) segCase |= 8;

  if (segCase==0 || segCase==15)
    return segments;

  if (segCase==5 || segCase==10)
  {
    segments.resize(2);
  } else
  {
    segments.resize(1);
  }

  switch (segCase)
  {
  case 1:
  case 14:
    vertlist[0] = VertexInterp(isoVal,isoComp,p0,p0d,p1,p1d,vertCache);
    vertlist[1] = VertexInterp(isoVal,isoComp,p3,p3d,p0,p0d,vertCache);
    break;
  case 2:
  case 13:
    vertlist[0] = VertexInterp(isoVal,isoComp,p0,p0d,p1,p1d,vertCache);
    vertlist[1] = VertexInterp(isoVal,isoComp,p1,p1d,p2,p2d,vertCache);
    break;
  case 3:
  case 12:
    vertlist[0] = VertexInterp(isoVal,isoComp,p1,p1d,p2,p2d,vertCache);
    vertlist[1] = VertexInterp(isoVal,isoComp,p3,p3d,p0,p0d,vertCache);
    break;
  case 4:
  case 11:
    vertlist[0] = VertexInterp(isoVal,isoComp,p1,p1d,p2,p2d,vertCache);
    vertlist[1] = VertexInterp(isoVal,isoComp,p2,p2d,p3,p3d,vertCache);
    break;
  case 6:
  case 9:
    vertlist[0] = VertexInterp(isoVal,isoComp,p0,p0d,p1,p1d,vertCache);
    vertlist[1] = VertexInterp(isoVal,isoComp,p2,p2d,p3,p3d,vertCache);
    break;
  case 7:
  case 8:
    vertlist[0] = VertexInterp(isoVal,isoComp,p2,p2d,p3,p3d,vertCache);
    vertlist[1] = VertexInterp(isoVal,isoComp,p3,p3d,p0,p0d,vertCache);
    break;
  case 5:
  case 10:
    vertlist[0] = VertexInterp(isoVal,isoComp,p0,p0d,p1,p1d,vertCache);
    vertlist[1] = VertexInterp(isoVal,isoComp,p1,p1d,p2,p2d,vertCache);
    vertlist[2] = VertexInterp(isoVal,isoComp,p2,p2d,p3,p3d,vertCache);
    vertlist[3] = VertexInterp(isoVal,isoComp,p3,p3d,p0,p0d,vertCache);
    break;
  }

  segments[0][0] = vertlist[0];
  segments[0][1] = vertlist[1];

  if (segCase==5 || segCase==10)
  {
    segments[1][0] = vertlist[2];
    segments[1][1] = vertlist[3];
  }

  return segments;
}

#else

/*
  Given a grid cell and an isoVal, calculate the triangular
  facets required to represent the isosurface through the cell.
  Return an array of (at most 5) triangular facets
*/
Vector<Triangle> Polygonise(const FArrayBox& pts,
                            const FArrayBox& mask,
                            PMap&            vertCache,
                            const IntVect&   baseIV,
                            Real             isoVal,
                            int              isoComp)
{
  int cubeindex;
  Vector<PMapIt> vertlist(12);
  Vector<Triangle> triangles; // result

  const IntVect& p0 = baseIV;
  const IntVect p1 = p0 + amrex::BASISV(0);
  const IntVect p2 = p1 + amrex::BASISV(1);
  const IntVect p3 = p0 + amrex::BASISV(1);
  const IntVect p4 = p0 + amrex::BASISV(2);
  const IntVect p5 = p4 + amrex::BASISV(0);
  const IntVect p6 = p5 + amrex::BASISV(1);
  const IntVect p7 = p4 + amrex::BASISV(1);

  // Bail if any of the point are masked out
  if (mask(p0)<0 || mask(p1)<0 || mask(p2)<0 || mask(p3)<0 ||
      mask(p4)<0 || mask(p5)<0 || mask(p6)<0 || mask(p7)<0 )
    return triangles;

  // Get coordinates of the cell vertices + any field variables
  const int nComp = pts.nComp();
  Point p0d(nComp); pts.getVal(p0d.dataPtr(),p0);
  Point p1d(nComp); pts.getVal(p1d.dataPtr(),p1);
  Point p2d(nComp); pts.getVal(p2d.dataPtr(),p2);
  Point p3d(nComp); pts.getVal(p3d.dataPtr(),p3);
  Point p4d(nComp); pts.getVal(p4d.dataPtr(),p4);
  Point p5d(nComp); pts.getVal(p5d.dataPtr(),p5);
  Point p6d(nComp); pts.getVal(p6d.dataPtr(),p6);
  Point p7d(nComp); pts.getVal(p7d.dataPtr(),p7);

  int edgeTable[256]={
    0x0  , 0x109, 0x203, 0x30a, 0x406, 0x50f, 0x605, 0x70c,
    0x80c, 0x905, 0xa0f, 0xb06, 0xc0a, 0xd03, 0xe09, 0xf00,
    0x190, 0x99 , 0x393, 0x29a, 0x596, 0x49f, 0x795, 0x69c,
    0x99c, 0x895, 0xb9f, 0xa96, 0xd9a, 0xc93, 0xf99, 0xe90,
    0x230, 0x339, 0x33 , 0x13a, 0x636, 0x73f, 0x435, 0x53c,
    0xa3c, 0xb35, 0x83f, 0x936, 0xe3a, 0xf33, 0xc39, 0xd30,
    0x3a0, 0x2a9, 0x1a3, 0xaa , 0x7a6, 0x6af, 0x5a5, 0x4ac,
    0xbac, 0xaa5, 0x9af, 0x8a6, 0xfaa, 0xea3, 0xda9, 0xca0,
    0x460, 0x569, 0x663, 0x76a, 0x66 , 0x16f, 0x265, 0x36c,
    0xc6c, 0xd65, 0xe6f, 0xf66, 0x86a, 0x963, 0xa69, 0xb60,
    0x5f0, 0x4f9, 0x7f3, 0x6fa, 0x1f6, 0xff , 0x3f5, 0x2fc,
    0xdfc, 0xcf5, 0xfff, 0xef6, 0x9fa, 0x8f3, 0xbf9, 0xaf0,
    0x650, 0x759, 0x453, 0x55a, 0x256, 0x35f, 0x55 , 0x15c,
    0xe5c, 0xf55, 0xc5f, 0xd56, 0xa5a, 0xb53, 0x859, 0x950,
    0x7c0, 0x6c9, 0x5c3, 0x4ca, 0x3c6, 0x2cf, 0x1c5, 0xcc ,
    0xfcc, 0xec5, 0xdcf, 0xcc6, 0xbca, 0xac3, 0x9c9, 0x8c0,
    0x8c0, 0x9c9, 0xac3, 0xbca, 0xcc6, 0xdcf, 0xec5, 0xfcc,
    0xcc , 0x1c5, 0x2cf, 0x3c6, 0x4ca, 0x5c3, 0x6c9, 0x7c0,
    0x950, 0x859, 0xb53, 0xa5a, 0xd56, 0xc5f, 0xf55, 0xe5c,
    0x15c, 0x55 , 0x35f, 0x256, 0x55a, 0x453, 0x759, 0x650,
    0xaf0, 0xbf9, 0x8f3, 0x9fa, 0xef6, 0xfff, 0xcf5, 0xdfc,
    0x2fc, 0x3f5, 0xff , 0x1f6, 0x6fa, 0x7f3, 0x4f9, 0x5f0,
    0xb60, 0xa69, 0x963, 0x86a, 0xf66, 0xe6f, 0xd65, 0xc6c,
    0x36c, 0x265, 0x16f, 0x66 , 0x76a, 0x663, 0x569, 0x460,
    0xca0, 0xda9, 0xea3, 0xfaa, 0x8a6, 0x9af, 0xaa5, 0xbac,
    0x4ac, 0x5a5, 0x6af, 0x7a6, 0xaa , 0x1a3, 0x2a9, 0x3a0,
    0xd30, 0xc39, 0xf33, 0xe3a, 0x936, 0x83f, 0xb35, 0xa3c,
    0x53c, 0x435, 0x73f, 0x636, 0x13a, 0x33 , 0x339, 0x230,
    0xe90, 0xf99, 0xc93, 0xd9a, 0xa96, 0xb9f, 0x895, 0x99c,
    0x69c, 0x795, 0x49f, 0x596, 0x29a, 0x393, 0x99 , 0x190,
    0xf00, 0xe09, 0xd03, 0xc0a, 0xb06, 0xa0f, 0x905, 0x80c,
    0x70c, 0x605, 0x50f, 0x406, 0x30a, 0x203, 0x109, 0x0   };

  int triTable[256][16] =
    {{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {0, 1, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {1, 8, 3, 9, 8, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {1, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {0, 8, 3, 1, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {9, 2, 10, 0, 2, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {2, 8, 3, 2, 10, 8, 10, 9, 8, -1, -1, -1, -1, -1, -1, -1},
     {3, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {0, 11, 2, 8, 11, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {1, 9, 0, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {1, 11, 2, 1, 9, 11, 9, 8, 11, -1, -1, -1, -1, -1, -1, -1},
     {3, 10, 1, 11, 10, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {0, 10, 1, 0, 8, 10, 8, 11, 10, -1, -1, -1, -1, -1, -1, -1},
     {3, 9, 0, 3, 11, 9, 11, 10, 9, -1, -1, -1, -1, -1, -1, -1},
     {9, 8, 10, 10, 8, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {4, 3, 0, 7, 3, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {0, 1, 9, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {4, 1, 9, 4, 7, 1, 7, 3, 1, -1, -1, -1, -1, -1, -1, -1},
     {1, 2, 10, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {3, 4, 7, 3, 0, 4, 1, 2, 10, -1, -1, -1, -1, -1, -1, -1},
     {9, 2, 10, 9, 0, 2, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1},
     {2, 10, 9, 2, 9, 7, 2, 7, 3, 7, 9, 4, -1, -1, -1, -1},
     {8, 4, 7, 3, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {11, 4, 7, 11, 2, 4, 2, 0, 4, -1, -1, -1, -1, -1, -1, -1},
     {9, 0, 1, 8, 4, 7, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1},
     {4, 7, 11, 9, 4, 11, 9, 11, 2, 9, 2, 1, -1, -1, -1, -1},
     {3, 10, 1, 3, 11, 10, 7, 8, 4, -1, -1, -1, -1, -1, -1, -1},
     {1, 11, 10, 1, 4, 11, 1, 0, 4, 7, 11, 4, -1, -1, -1, -1},
     {4, 7, 8, 9, 0, 11, 9, 11, 10, 11, 0, 3, -1, -1, -1, -1},
     {4, 7, 11, 4, 11, 9, 9, 11, 10, -1, -1, -1, -1, -1, -1, -1},
     {9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {9, 5, 4, 0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {0, 5, 4, 1, 5, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {8, 5, 4, 8, 3, 5, 3, 1, 5, -1, -1, -1, -1, -1, -1, -1},
     {1, 2, 10, 9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {3, 0, 8, 1, 2, 10, 4, 9, 5, -1, -1, -1, -1, -1, -1, -1},
     {5, 2, 10, 5, 4, 2, 4, 0, 2, -1, -1, -1, -1, -1, -1, -1},
     {2, 10, 5, 3, 2, 5, 3, 5, 4, 3, 4, 8, -1, -1, -1, -1},
     {9, 5, 4, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {0, 11, 2, 0, 8, 11, 4, 9, 5, -1, -1, -1, -1, -1, -1, -1},
     {0, 5, 4, 0, 1, 5, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1},
     {2, 1, 5, 2, 5, 8, 2, 8, 11, 4, 8, 5, -1, -1, -1, -1},
     {10, 3, 11, 10, 1, 3, 9, 5, 4, -1, -1, -1, -1, -1, -1, -1},
     {4, 9, 5, 0, 8, 1, 8, 10, 1, 8, 11, 10, -1, -1, -1, -1},
     {5, 4, 0, 5, 0, 11, 5, 11, 10, 11, 0, 3, -1, -1, -1, -1},
     {5, 4, 8, 5, 8, 10, 10, 8, 11, -1, -1, -1, -1, -1, -1, -1},
     {9, 7, 8, 5, 7, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {9, 3, 0, 9, 5, 3, 5, 7, 3, -1, -1, -1, -1, -1, -1, -1},
     {0, 7, 8, 0, 1, 7, 1, 5, 7, -1, -1, -1, -1, -1, -1, -1},
     {1, 5, 3, 3, 5, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {9, 7, 8, 9, 5, 7, 10, 1, 2, -1, -1, -1, -1, -1, -1, -1},
     {10, 1, 2, 9, 5, 0, 5, 3, 0, 5, 7, 3, -1, -1, -1, -1},
     {8, 0, 2, 8, 2, 5, 8, 5, 7, 10, 5, 2, -1, -1, -1, -1},
     {2, 10, 5, 2, 5, 3, 3, 5, 7, -1, -1, -1, -1, -1, -1, -1},
     {7, 9, 5, 7, 8, 9, 3, 11, 2, -1, -1, -1, -1, -1, -1, -1},
     {9, 5, 7, 9, 7, 2, 9, 2, 0, 2, 7, 11, -1, -1, -1, -1},
     {2, 3, 11, 0, 1, 8, 1, 7, 8, 1, 5, 7, -1, -1, -1, -1},
     {11, 2, 1, 11, 1, 7, 7, 1, 5, -1, -1, -1, -1, -1, -1, -1},
     {9, 5, 8, 8, 5, 7, 10, 1, 3, 10, 3, 11, -1, -1, -1, -1},
     {5, 7, 0, 5, 0, 9, 7, 11, 0, 1, 0, 10, 11, 10, 0, -1},
     {11, 10, 0, 11, 0, 3, 10, 5, 0, 8, 0, 7, 5, 7, 0, -1},
     {11, 10, 5, 7, 11, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {10, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {0, 8, 3, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {9, 0, 1, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {1, 8, 3, 1, 9, 8, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1},
     {1, 6, 5, 2, 6, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {1, 6, 5, 1, 2, 6, 3, 0, 8, -1, -1, -1, -1, -1, -1, -1},
     {9, 6, 5, 9, 0, 6, 0, 2, 6, -1, -1, -1, -1, -1, -1, -1},
     {5, 9, 8, 5, 8, 2, 5, 2, 6, 3, 2, 8, -1, -1, -1, -1},
     {2, 3, 11, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {11, 0, 8, 11, 2, 0, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1},
     {0, 1, 9, 2, 3, 11, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1},
     {5, 10, 6, 1, 9, 2, 9, 11, 2, 9, 8, 11, -1, -1, -1, -1},
     {6, 3, 11, 6, 5, 3, 5, 1, 3, -1, -1, -1, -1, -1, -1, -1},
     {0, 8, 11, 0, 11, 5, 0, 5, 1, 5, 11, 6, -1, -1, -1, -1},
     {3, 11, 6, 0, 3, 6, 0, 6, 5, 0, 5, 9, -1, -1, -1, -1},
     {6, 5, 9, 6, 9, 11, 11, 9, 8, -1, -1, -1, -1, -1, -1, -1},
     {5, 10, 6, 4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {4, 3, 0, 4, 7, 3, 6, 5, 10, -1, -1, -1, -1, -1, -1, -1},
     {1, 9, 0, 5, 10, 6, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1},
     {10, 6, 5, 1, 9, 7, 1, 7, 3, 7, 9, 4, -1, -1, -1, -1},
     {6, 1, 2, 6, 5, 1, 4, 7, 8, -1, -1, -1, -1, -1, -1, -1},
     {1, 2, 5, 5, 2, 6, 3, 0, 4, 3, 4, 7, -1, -1, -1, -1},
     {8, 4, 7, 9, 0, 5, 0, 6, 5, 0, 2, 6, -1, -1, -1, -1},
     {7, 3, 9, 7, 9, 4, 3, 2, 9, 5, 9, 6, 2, 6, 9, -1},
     {3, 11, 2, 7, 8, 4, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1},
     {5, 10, 6, 4, 7, 2, 4, 2, 0, 2, 7, 11, -1, -1, -1, -1},
     {0, 1, 9, 4, 7, 8, 2, 3, 11, 5, 10, 6, -1, -1, -1, -1},
     {9, 2, 1, 9, 11, 2, 9, 4, 11, 7, 11, 4, 5, 10, 6, -1},
     {8, 4, 7, 3, 11, 5, 3, 5, 1, 5, 11, 6, -1, -1, -1, -1},
     {5, 1, 11, 5, 11, 6, 1, 0, 11, 7, 11, 4, 0, 4, 11, -1},
     {0, 5, 9, 0, 6, 5, 0, 3, 6, 11, 6, 3, 8, 4, 7, -1},
     {6, 5, 9, 6, 9, 11, 4, 7, 9, 7, 11, 9, -1, -1, -1, -1},
     {10, 4, 9, 6, 4, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {4, 10, 6, 4, 9, 10, 0, 8, 3, -1, -1, -1, -1, -1, -1, -1},
     {10, 0, 1, 10, 6, 0, 6, 4, 0, -1, -1, -1, -1, -1, -1, -1},
     {8, 3, 1, 8, 1, 6, 8, 6, 4, 6, 1, 10, -1, -1, -1, -1},
     {1, 4, 9, 1, 2, 4, 2, 6, 4, -1, -1, -1, -1, -1, -1, -1},
     {3, 0, 8, 1, 2, 9, 2, 4, 9, 2, 6, 4, -1, -1, -1, -1},
     {0, 2, 4, 4, 2, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {8, 3, 2, 8, 2, 4, 4, 2, 6, -1, -1, -1, -1, -1, -1, -1},
     {10, 4, 9, 10, 6, 4, 11, 2, 3, -1, -1, -1, -1, -1, -1, -1},
     {0, 8, 2, 2, 8, 11, 4, 9, 10, 4, 10, 6, -1, -1, -1, -1},
     {3, 11, 2, 0, 1, 6, 0, 6, 4, 6, 1, 10, -1, -1, -1, -1},
     {6, 4, 1, 6, 1, 10, 4, 8, 1, 2, 1, 11, 8, 11, 1, -1},
     {9, 6, 4, 9, 3, 6, 9, 1, 3, 11, 6, 3, -1, -1, -1, -1},
     {8, 11, 1, 8, 1, 0, 11, 6, 1, 9, 1, 4, 6, 4, 1, -1},
     {3, 11, 6, 3, 6, 0, 0, 6, 4, -1, -1, -1, -1, -1, -1, -1},
     {6, 4, 8, 11, 6, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {7, 10, 6, 7, 8, 10, 8, 9, 10, -1, -1, -1, -1, -1, -1, -1},
     {0, 7, 3, 0, 10, 7, 0, 9, 10, 6, 7, 10, -1, -1, -1, -1},
     {10, 6, 7, 1, 10, 7, 1, 7, 8, 1, 8, 0, -1, -1, -1, -1},
     {10, 6, 7, 10, 7, 1, 1, 7, 3, -1, -1, -1, -1, -1, -1, -1},
     {1, 2, 6, 1, 6, 8, 1, 8, 9, 8, 6, 7, -1, -1, -1, -1},
     {2, 6, 9, 2, 9, 1, 6, 7, 9, 0, 9, 3, 7, 3, 9, -1},
     {7, 8, 0, 7, 0, 6, 6, 0, 2, -1, -1, -1, -1, -1, -1, -1},
     {7, 3, 2, 6, 7, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {2, 3, 11, 10, 6, 8, 10, 8, 9, 8, 6, 7, -1, -1, -1, -1},
     {2, 0, 7, 2, 7, 11, 0, 9, 7, 6, 7, 10, 9, 10, 7, -1},
     {1, 8, 0, 1, 7, 8, 1, 10, 7, 6, 7, 10, 2, 3, 11, -1},
     {11, 2, 1, 11, 1, 7, 10, 6, 1, 6, 7, 1, -1, -1, -1, -1},
     {8, 9, 6, 8, 6, 7, 9, 1, 6, 11, 6, 3, 1, 3, 6, -1},
     {0, 9, 1, 11, 6, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {7, 8, 0, 7, 0, 6, 3, 11, 0, 11, 6, 0, -1, -1, -1, -1},
     {7, 11, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {7, 6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {3, 0, 8, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {0, 1, 9, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {8, 1, 9, 8, 3, 1, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1},
     {10, 1, 2, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {1, 2, 10, 3, 0, 8, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1},
     {2, 9, 0, 2, 10, 9, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1},
     {6, 11, 7, 2, 10, 3, 10, 8, 3, 10, 9, 8, -1, -1, -1, -1},
     {7, 2, 3, 6, 2, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {7, 0, 8, 7, 6, 0, 6, 2, 0, -1, -1, -1, -1, -1, -1, -1},
     {2, 7, 6, 2, 3, 7, 0, 1, 9, -1, -1, -1, -1, -1, -1, -1},
     {1, 6, 2, 1, 8, 6, 1, 9, 8, 8, 7, 6, -1, -1, -1, -1},
     {10, 7, 6, 10, 1, 7, 1, 3, 7, -1, -1, -1, -1, -1, -1, -1},
     {10, 7, 6, 1, 7, 10, 1, 8, 7, 1, 0, 8, -1, -1, -1, -1},
     {0, 3, 7, 0, 7, 10, 0, 10, 9, 6, 10, 7, -1, -1, -1, -1},
     {7, 6, 10, 7, 10, 8, 8, 10, 9, -1, -1, -1, -1, -1, -1, -1},
     {6, 8, 4, 11, 8, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {3, 6, 11, 3, 0, 6, 0, 4, 6, -1, -1, -1, -1, -1, -1, -1},
     {8, 6, 11, 8, 4, 6, 9, 0, 1, -1, -1, -1, -1, -1, -1, -1},
     {9, 4, 6, 9, 6, 3, 9, 3, 1, 11, 3, 6, -1, -1, -1, -1},
     {6, 8, 4, 6, 11, 8, 2, 10, 1, -1, -1, -1, -1, -1, -1, -1},
     {1, 2, 10, 3, 0, 11, 0, 6, 11, 0, 4, 6, -1, -1, -1, -1},
     {4, 11, 8, 4, 6, 11, 0, 2, 9, 2, 10, 9, -1, -1, -1, -1},
     {10, 9, 3, 10, 3, 2, 9, 4, 3, 11, 3, 6, 4, 6, 3, -1},
     {8, 2, 3, 8, 4, 2, 4, 6, 2, -1, -1, -1, -1, -1, -1, -1},
     {0, 4, 2, 4, 6, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {1, 9, 0, 2, 3, 4, 2, 4, 6, 4, 3, 8, -1, -1, -1, -1},
     {1, 9, 4, 1, 4, 2, 2, 4, 6, -1, -1, -1, -1, -1, -1, -1},
     {8, 1, 3, 8, 6, 1, 8, 4, 6, 6, 10, 1, -1, -1, -1, -1},
     {10, 1, 0, 10, 0, 6, 6, 0, 4, -1, -1, -1, -1, -1, -1, -1},
     {4, 6, 3, 4, 3, 8, 6, 10, 3, 0, 3, 9, 10, 9, 3, -1},
     {10, 9, 4, 6, 10, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {4, 9, 5, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {0, 8, 3, 4, 9, 5, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1},
     {5, 0, 1, 5, 4, 0, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1},
     {11, 7, 6, 8, 3, 4, 3, 5, 4, 3, 1, 5, -1, -1, -1, -1},
     {9, 5, 4, 10, 1, 2, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1},
     {6, 11, 7, 1, 2, 10, 0, 8, 3, 4, 9, 5, -1, -1, -1, -1},
     {7, 6, 11, 5, 4, 10, 4, 2, 10, 4, 0, 2, -1, -1, -1, -1},
     {3, 4, 8, 3, 5, 4, 3, 2, 5, 10, 5, 2, 11, 7, 6, -1},
     {7, 2, 3, 7, 6, 2, 5, 4, 9, -1, -1, -1, -1, -1, -1, -1},
     {9, 5, 4, 0, 8, 6, 0, 6, 2, 6, 8, 7, -1, -1, -1, -1},
     {3, 6, 2, 3, 7, 6, 1, 5, 0, 5, 4, 0, -1, -1, -1, -1},
     {6, 2, 8, 6, 8, 7, 2, 1, 8, 4, 8, 5, 1, 5, 8, -1},
     {9, 5, 4, 10, 1, 6, 1, 7, 6, 1, 3, 7, -1, -1, -1, -1},
     {1, 6, 10, 1, 7, 6, 1, 0, 7, 8, 7, 0, 9, 5, 4, -1},
     {4, 0, 10, 4, 10, 5, 0, 3, 10, 6, 10, 7, 3, 7, 10, -1},
     {7, 6, 10, 7, 10, 8, 5, 4, 10, 4, 8, 10, -1, -1, -1, -1},
     {6, 9, 5, 6, 11, 9, 11, 8, 9, -1, -1, -1, -1, -1, -1, -1},
     {3, 6, 11, 0, 6, 3, 0, 5, 6, 0, 9, 5, -1, -1, -1, -1},
     {0, 11, 8, 0, 5, 11, 0, 1, 5, 5, 6, 11, -1, -1, -1, -1},
     {6, 11, 3, 6, 3, 5, 5, 3, 1, -1, -1, -1, -1, -1, -1, -1},
     {1, 2, 10, 9, 5, 11, 9, 11, 8, 11, 5, 6, -1, -1, -1, -1},
     {0, 11, 3, 0, 6, 11, 0, 9, 6, 5, 6, 9, 1, 2, 10, -1},
     {11, 8, 5, 11, 5, 6, 8, 0, 5, 10, 5, 2, 0, 2, 5, -1},
     {6, 11, 3, 6, 3, 5, 2, 10, 3, 10, 5, 3, -1, -1, -1, -1},
     {5, 8, 9, 5, 2, 8, 5, 6, 2, 3, 8, 2, -1, -1, -1, -1},
     {9, 5, 6, 9, 6, 0, 0, 6, 2, -1, -1, -1, -1, -1, -1, -1},
     {1, 5, 8, 1, 8, 0, 5, 6, 8, 3, 8, 2, 6, 2, 8, -1},
     {1, 5, 6, 2, 1, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {1, 3, 6, 1, 6, 10, 3, 8, 6, 5, 6, 9, 8, 9, 6, -1},
     {10, 1, 0, 10, 0, 6, 9, 5, 0, 5, 6, 0, -1, -1, -1, -1},
     {0, 3, 8, 5, 6, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {10, 5, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {11, 5, 10, 7, 5, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {11, 5, 10, 11, 7, 5, 8, 3, 0, -1, -1, -1, -1, -1, -1, -1},
     {5, 11, 7, 5, 10, 11, 1, 9, 0, -1, -1, -1, -1, -1, -1, -1},
     {10, 7, 5, 10, 11, 7, 9, 8, 1, 8, 3, 1, -1, -1, -1, -1},
     {11, 1, 2, 11, 7, 1, 7, 5, 1, -1, -1, -1, -1, -1, -1, -1},
     {0, 8, 3, 1, 2, 7, 1, 7, 5, 7, 2, 11, -1, -1, -1, -1},
     {9, 7, 5, 9, 2, 7, 9, 0, 2, 2, 11, 7, -1, -1, -1, -1},
     {7, 5, 2, 7, 2, 11, 5, 9, 2, 3, 2, 8, 9, 8, 2, -1},
     {2, 5, 10, 2, 3, 5, 3, 7, 5, -1, -1, -1, -1, -1, -1, -1},
     {8, 2, 0, 8, 5, 2, 8, 7, 5, 10, 2, 5, -1, -1, -1, -1},
     {9, 0, 1, 5, 10, 3, 5, 3, 7, 3, 10, 2, -1, -1, -1, -1},
     {9, 8, 2, 9, 2, 1, 8, 7, 2, 10, 2, 5, 7, 5, 2, -1},
     {1, 3, 5, 3, 7, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {0, 8, 7, 0, 7, 1, 1, 7, 5, -1, -1, -1, -1, -1, -1, -1},
     {9, 0, 3, 9, 3, 5, 5, 3, 7, -1, -1, -1, -1, -1, -1, -1},
     {9, 8, 7, 5, 9, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {5, 8, 4, 5, 10, 8, 10, 11, 8, -1, -1, -1, -1, -1, -1, -1},
     {5, 0, 4, 5, 11, 0, 5, 10, 11, 11, 3, 0, -1, -1, -1, -1},
     {0, 1, 9, 8, 4, 10, 8, 10, 11, 10, 4, 5, -1, -1, -1, -1},
     {10, 11, 4, 10, 4, 5, 11, 3, 4, 9, 4, 1, 3, 1, 4, -1},
     {2, 5, 1, 2, 8, 5, 2, 11, 8, 4, 5, 8, -1, -1, -1, -1},
     {0, 4, 11, 0, 11, 3, 4, 5, 11, 2, 11, 1, 5, 1, 11, -1},
     {0, 2, 5, 0, 5, 9, 2, 11, 5, 4, 5, 8, 11, 8, 5, -1},
     {9, 4, 5, 2, 11, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {2, 5, 10, 3, 5, 2, 3, 4, 5, 3, 8, 4, -1, -1, -1, -1},
     {5, 10, 2, 5, 2, 4, 4, 2, 0, -1, -1, -1, -1, -1, -1, -1},
     {3, 10, 2, 3, 5, 10, 3, 8, 5, 4, 5, 8, 0, 1, 9, -1},
     {5, 10, 2, 5, 2, 4, 1, 9, 2, 9, 4, 2, -1, -1, -1, -1},
     {8, 4, 5, 8, 5, 3, 3, 5, 1, -1, -1, -1, -1, -1, -1, -1},
     {0, 4, 5, 1, 0, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {8, 4, 5, 8, 5, 3, 9, 0, 5, 0, 3, 5, -1, -1, -1, -1},
     {9, 4, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {4, 11, 7, 4, 9, 11, 9, 10, 11, -1, -1, -1, -1, -1, -1, -1},
     {0, 8, 3, 4, 9, 7, 9, 11, 7, 9, 10, 11, -1, -1, -1, -1},
     {1, 10, 11, 1, 11, 4, 1, 4, 0, 7, 4, 11, -1, -1, -1, -1},
     {3, 1, 4, 3, 4, 8, 1, 10, 4, 7, 4, 11, 10, 11, 4, -1},
     {4, 11, 7, 9, 11, 4, 9, 2, 11, 9, 1, 2, -1, -1, -1, -1},
     {9, 7, 4, 9, 11, 7, 9, 1, 11, 2, 11, 1, 0, 8, 3, -1},
     {11, 7, 4, 11, 4, 2, 2, 4, 0, -1, -1, -1, -1, -1, -1, -1},
     {11, 7, 4, 11, 4, 2, 8, 3, 4, 3, 2, 4, -1, -1, -1, -1},
     {2, 9, 10, 2, 7, 9, 2, 3, 7, 7, 4, 9, -1, -1, -1, -1},
     {9, 10, 7, 9, 7, 4, 10, 2, 7, 8, 7, 0, 2, 0, 7, -1},
     {3, 7, 10, 3, 10, 2, 7, 4, 10, 1, 10, 0, 4, 0, 10, -1},
     {1, 10, 2, 8, 7, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {4, 9, 1, 4, 1, 7, 7, 1, 3, -1, -1, -1, -1, -1, -1, -1},
     {4, 9, 1, 4, 1, 7, 0, 8, 1, 8, 7, 1, -1, -1, -1, -1},
     {4, 0, 3, 7, 4, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {4, 8, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {9, 10, 8, 10, 11, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {3, 0, 9, 3, 9, 11, 11, 9, 10, -1, -1, -1, -1, -1, -1, -1},
     {0, 1, 10, 0, 10, 8, 8, 10, 11, -1, -1, -1, -1, -1, -1, -1},
     {3, 1, 10, 11, 3, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {1, 2, 11, 1, 11, 9, 9, 11, 8, -1, -1, -1, -1, -1, -1, -1},
     {3, 0, 9, 3, 9, 11, 1, 2, 9, 2, 11, 9, -1, -1, -1, -1},
     {0, 2, 11, 8, 0, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {3, 2, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {2, 3, 8, 2, 8, 10, 10, 8, 9, -1, -1, -1, -1, -1, -1, -1},
     {9, 10, 2, 0, 9, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {2, 3, 8, 2, 8, 10, 0, 1, 8, 1, 10, 8, -1, -1, -1, -1},
     {1, 10, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {1, 3, 8, 9, 1, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {0, 9, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {0, 3, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}};

  /*
    Determine the index into the edge table which
    tells us which vertices are inside of the surface
  */
  cubeindex = 0;
  if (p0d[isoComp] < isoVal) cubeindex |= 1;
  if (p1d[isoComp] < isoVal) cubeindex |= 2;
  if (p2d[isoComp] < isoVal) cubeindex |= 4;
  if (p3d[isoComp] < isoVal) cubeindex |= 8;
  if (p4d[isoComp] < isoVal) cubeindex |= 16;
  if (p5d[isoComp] < isoVal) cubeindex |= 32;
  if (p6d[isoComp] < isoVal) cubeindex |= 64;
  if (p7d[isoComp] < isoVal) cubeindex |= 128;

  /* Cube is entirely in/out of the surface */
  if (edgeTable[cubeindex] == 0)
    return triangles;

  /* Find the vertices where the surface intersects the cube */
  if (edgeTable[cubeindex] & 1)
    vertlist[0]  = VertexInterp(isoVal,isoComp,p0,p0d,p1,p1d,vertCache);
  if (edgeTable[cubeindex] & 2)
    vertlist[1]  = VertexInterp(isoVal,isoComp,p1,p1d,p2,p2d,vertCache);
  if (edgeTable[cubeindex] & 4)
    vertlist[2]  = VertexInterp(isoVal,isoComp,p2,p2d,p3,p3d,vertCache);
  if (edgeTable[cubeindex] & 8)
    vertlist[3]  = VertexInterp(isoVal,isoComp,p3,p3d,p0,p0d,vertCache);
  if (edgeTable[cubeindex] & 16)
    vertlist[4]  = VertexInterp(isoVal,isoComp,p4,p4d,p5,p5d,vertCache);
  if (edgeTable[cubeindex] & 32)
    vertlist[5]  = VertexInterp(isoVal,isoComp,p5,p5d,p6,p6d,vertCache);
  if (edgeTable[cubeindex] & 64)
    vertlist[6]  = VertexInterp(isoVal,isoComp,p6,p6d,p7,p7d,vertCache);
  if (edgeTable[cubeindex] & 128)
    vertlist[7]  = VertexInterp(isoVal,isoComp,p7,p7d,p4,p4d,vertCache);
  if (edgeTable[cubeindex] & 256)
    vertlist[8]  = VertexInterp(isoVal,isoComp,p0,p0d,p4,p4d,vertCache);
  if (edgeTable[cubeindex] & 512)
    vertlist[9]  = VertexInterp(isoVal,isoComp,p1,p1d,p5,p5d,vertCache);
  if (edgeTable[cubeindex] & 1024)
    vertlist[10]  = VertexInterp(isoVal,isoComp,p2,p2d,p6,p6d,vertCache);
  if (edgeTable[cubeindex] & 2048)
    vertlist[11]  = VertexInterp(isoVal,isoComp,p3,p3d,p7,p7d,vertCache);

  /* Create the triangles */
  int nTriang = 0;
  for (int i=0;triTable[cubeindex][i]!=-1;i+=3)
    nTriang++;

  triangles.resize(nTriang);
  for (int j=0; j<nTriang; ++j)
  {
    int j3 = 3*j;
    triangles[j][0] = vertlist[triTable[cubeindex][j3  ]];
    triangles[j][1] = vertlist[triTable[cubeindex][j3+1]];
    triangles[j][2] = vertlist[triTable[cubeindex][j3+2]];
  }

  return triangles;
}
#endif

struct Node
{
  // Base constructor
  Node() { m_vec = nullptr; }
  Node(const vector<Real>& vec, long idx) : m_idx(idx)
  {
    m_size = vec.size();
    m_vec = new Real[m_size];
    for (int i=0; i<m_size; ++i) {
      m_vec[i] = vec[i];
    }
  }
  // Copy constructor
  Node(const Node& rhs)
  {
    m_size = rhs.m_size;
    m_idx = rhs.m_idx;
    m_vec = new Real[m_size];
    for (int i=0; i<m_size; ++i)
      m_vec[i] = rhs.m_vec[i];
  }

  // Destructor
  ~Node() { delete [] m_vec; }

  // Accessor to spatial coordinates
  Real operator[] (int n) const { return m_vec[n]; }

  // Relational operator to enable std::set
  inline bool operator< (const Node& rhs) const {
    Real sum = 0.0;
    for (int i=0; i<AMREX_SPACEDIM; ++i) {
      sum += (m_vec[i]-rhs[i])*(m_vec[i]-rhs[i]);
    }

    if (myLessVerbose)
    {
      std::cerr << "In LT" << std::endl;
      std::cerr << "  lhs: " << m_vec[0] << " " << m_vec[1] << " " << m_vec[2] << std::endl;
      std::cerr << "  rhs: " << rhs[0] << " " << rhs[1] << " " << rhs[2] << std::endl;
      if (std::sqrt(sum) < epsilon_DEF)
      {
        std::cerr << "under eps tol" << std::endl;
      }
      else
      {
        if (m_vec[0]==rhs[0]) {
          if (m_vec[1]==rhs[1]) {
            std::cerr << "          based on z" << std::endl;
          } else {
            std::cerr << "          based on y" << std::endl;
          }
        } else {
          std::cerr << "          based on x" << std::endl;
        }
      }
    }
    if (std::sqrt(sum) < epsilon_DEF) return false;

    if (m_vec[0]==rhs[0]) {
      if (m_vec[1]==rhs[1]) {
        return m_vec[2] < rhs[2];
      } else {
        return m_vec[1] < rhs[1];
      }
    } else {
      return m_vec[0] < rhs[0];
    }
  }

  // Data
  int m_size;
  int m_idx;
  Real* m_vec;

  static bool myLessVerbose;
};

bool Node::myLessVerbose=false;

// Array of node indices forming an element,
// segment or triangle in 2D and 3D, resp.
struct Element
{
  Element (const vector<int>& vec)
  {
    AMREX_D_TERM(m_vec[0] = vec[0];,
                 m_vec[1] = vec[1];,
                 m_vec[2] = vec[2]);
    //
    // Rotate elements so that first one is the smallest.
    //
    int smallest = 0;
    for (int i=1; i<AMREX_SPACEDIM; ++i) {
      if (m_vec[smallest]>m_vec[i]) smallest=i;
    }
    std::rotate(&m_vec[0],&m_vec[smallest],&m_vec[0]+AMREX_SPACEDIM);
  }

  int size () const {return AMREX_SPACEDIM;}

  // Accessor to node indices
  int operator[] (int n) const
  {
    AMREX_ASSERT(n >= 0 && n < AMREX_SPACEDIM);
    return m_vec[n];
  }
  bool operator< (const Element& rhs) const
  {
    if (m_vec[0]==rhs[0]) {
#if AMREX_SPACEDIM==3
      if (m_vec[1]==rhs[1]) {
        return m_vec[2] < rhs[2];
      } else
#endif
      {
        return m_vec[1] < rhs[1];
      }
    } else {
      return m_vec[0] < rhs[0];
    }
  }

  // Data
  int m_vec[AMREX_SPACEDIM];
};

void
Collate(Vector<Real>& NodeRaw,
        Vector<int>&  EltRaw,
        int          nCompPerNode)
{

#if AMREX_USE_MPI

  const int nProcs = ParallelDescriptor::NProcs();

  if (nProcs < 2) return;

  const int IOProc = ParallelDescriptor::IOProcessorNumber();

  AMREX_ASSERT(IOProc==0);

  Vector<int> nmdataR(nProcs,0);
  Vector<int> offsetR(nProcs,0);
  //
  // Tell root CPU how many Real data elements each CPU will be sending.
  //
  int countR = NodeRaw.size();
  MPI_Gather(&countR,
             1,
             ParallelDescriptor::Mpi_typemap<int>::type(),
             nmdataR.dataPtr(),
             1,
             ParallelDescriptor::Mpi_typemap<int>::type(),
             IOProc,
             ParallelDescriptor::Communicator());

  Vector<Real> nrRecv;
  if (ParallelDescriptor::IOProcessor())
  {
    for (int i = 1; i < nProcs; i++)
      offsetR[i] = offsetR[i-1] + nmdataR[i-1];

    nrRecv.resize(offsetR[nProcs-1] + nmdataR[nProcs-1]);
  }
  //
  // Gather all the Real data to IOProc into NodeRaw.
  //
  MPI_Gatherv(NodeRaw.dataPtr(),
              countR,
              ParallelDescriptor::Mpi_typemap<Real>::type(),
              nrRecv.dataPtr(),
              nmdataR.dataPtr(),
              offsetR.dataPtr(),
              ParallelDescriptor::Mpi_typemap<Real>::type(),
              IOProc,
              ParallelDescriptor::Communicator());

  NodeRaw.swap(nrRecv);
  nrRecv.clear();

  // Now communicate the element info
  Vector<int> nmdataI(nProcs,0);
  Vector<int> offsetI(nProcs,0);
  int countI = EltRaw.size();
  MPI_Gather(&countI,
             1,
             ParallelDescriptor::Mpi_typemap<int>::type(),
             nmdataI.dataPtr(),
             1,
             ParallelDescriptor::Mpi_typemap<int>::type(),
             IOProc,
             ParallelDescriptor::Communicator());

  Vector<int> erRecv;
  if (ParallelDescriptor::IOProcessor())
  {
    for (int i = 1; i < nProcs; i++)
      offsetI[i] = offsetI[i-1] + nmdataI[i-1];

    erRecv.resize(offsetI[nProcs-1] + nmdataI[nProcs-1]);
  }
  //
  // Gather all the data to IOProc into EltRaw
  //
  MPI_Gatherv(EltRaw.dataPtr(),
              countI,
              ParallelDescriptor::Mpi_typemap<int>::type(),
              erRecv.dataPtr(),
              nmdataI.dataPtr(),
              offsetI.dataPtr(),
              ParallelDescriptor::Mpi_typemap<int>::type(),
              IOProc,
              ParallelDescriptor::Communicator());

  EltRaw.swap(erRecv);
  erRecv.clear();

  //
  // Shift nodeIDs in element definitions
  //
  if (ParallelDescriptor::IOProcessor())
  {
    for (int i = 1; i < nProcs; i++)
    {
      const int nodeOffset = offsetR[i]/nCompPerNode;
      for (int j = 0; j < nmdataI[i]; j++)
        EltRaw[offsetI[i]+j] += nodeOffset;
    }
  }
#endif
}

class FABdata
{
public:
  FABdata(size_t i, int n)
    {fab.resize(Box(IntVect::TheZeroVector(),
                    IntVect(AMREX_D_DECL(i-1,0,0))),n);}
  Real& operator[](size_t i) {return fab.dataPtr()[i];}
  const FArrayBox& getFab() const {return fab;}
  FArrayBox fab;
  size_t boxSize;
};



struct Seg
{
    Seg() : mL(-1), mR(-1) {}
    Seg(int L, int R) : mL(L), mR(R) {}
    int ID_l() const
        { return mL; }
    int ID_r() const
        { return mR; }
    void flip()
        { int t=mL; mL=mR; mR=t;}
    int mL, mR;
};

typedef list<Seg> Line;

Line::iterator FindMySeg(Line& segs, int idx)
{
    for (Line::iterator it=segs.begin(); it!=segs.end(); ++it)
    {
        if ( ((*it).ID_l() == idx) || ((*it).ID_r() == idx) )
            return it;
    }
    return segs.end();
}

bool operator==(const Seg& lhs, const Seg& rhs)
{
    return lhs.ID_l() == rhs.ID_l() && lhs.ID_r() == rhs.ID_r();
}

bool operator!=(const Seg& lhs, const Seg& rhs)
{
    return !operator==(lhs,rhs);
}

typedef Vector<Real> Point;

std::pair<bool,Point>
nodesAreClose(const FArrayBox& nodes,
              int              lhs,
              int              rhs,
              Real             eps)
{
    Real sep = 0;
    int nComp = nodes.nComp();
    Point ret(nComp,0);
    for (int i=0; i<ret.size(); ++i)
    {
        Real xL = nodes.dataPtr(i)[lhs];
        Real xR = nodes.dataPtr(i)[rhs];
        sep += (xL-xR)*(xL-xR);
        ret[i] = 0.5*(xL+xR);
    }
    return std::pair<bool,Point>(std::sqrt(sep) <= eps, ret);
}

void
join(Line&     l1,
     Line&     l2,
     FArrayBox&   nodes,
     Point&       pt,
     list<Point>& newPts,
     list<int>&   deadPts)
{
    int newPtID = nodes.box().numPts() + newPts.size();
    newPts.push_back(pt);

    int rhs = l1.back().ID_l();
    deadPts.push_back(l1.back().ID_r());
    AMREX_ASSERT(l1.size()>0);
    l1.pop_back();
    l1.push_back(Seg(rhs,newPtID));

    int lhs = l2.front().ID_r();
    deadPts.push_back(l2.front().ID_l());
    AMREX_ASSERT(l2.size()>0);
    l2.pop_front();
    l2.push_front(Seg(newPtID,lhs));

    l1.splice(l1.end(),l2);
}

// Provides an ordering to uniquely sort segments
struct SegLT
{
    bool operator()(const Seg& lhs,
                    const Seg& rhs)
        {
            int smL = std::min(lhs.ID_l(),lhs.ID_r());
            int smR = std::min(rhs.ID_l(),rhs.ID_r());
            if (smL < smR)
            {
                return true;
            }
            else if (smL == smR)
            {
                int lgL = std::max(lhs.ID_l(),lhs.ID_r());
                int lgR = std::max(rhs.ID_l(),rhs.ID_r());
                return lgL < lgR;
            }
            return false;
        }
};

using std::list;

list<Line>
MakeCLines(const FArrayBox&  nodes,
           Vector<int>& faceData,
           int         nodesPerElt)
{
    // We want to clean things up a little...the contour is now in a form
    // of a huge list of itty-bitty segments.  We want to connect up the segments
    // to form a minimal number of polylines.

    // First build up a Line
    int nElts = faceData.size() / nodesPerElt;
    AMREX_ASSERT(nElts * nodesPerElt == faceData.size());
    Line segList;
    for (int i=0; i<nElts; ++i)
    {
        int offset = i*nodesPerElt;
        segList.push_back(Seg(faceData[offset]-1,faceData[offset+1]-1)); // make 0-based
    }

    // Find a segment with the specified vertex as one of its endpoints,
    // then assemble the list of segments to form the contour line.  If
    // we finish, and segments remain, start a new line.
    int idx = segList.front().ID_r();
    segList.pop_front();
    list<Line> cLines;
    cLines.push_back(Line());

    while (segList.begin() != segList.end())
    {
        Line::iterator segIt = FindMySeg(segList,idx);
        if (segIt != segList.end())
        {
            int idx_l = (*segIt).ID_l();
            int idx_r = (*segIt).ID_r();

            if ( idx_l == idx )
            {
                idx = idx_r;
                cLines.back().push_back(*segIt);
            }
            else
            {
                idx = idx_l;
                Seg nseg(*segIt); nseg.flip();
                cLines.back().push_back(nseg);
            }

            segList.erase(segIt);
        }
        else
        {
            cLines.push_back(Line());
            idx = segList.front().ID_r();
            segList.pop_front();
        }
    }

    // Connect up the line fragments as much as possible
    bool changed;
    do
    {
        changed = false;
        for (std::list<Line>::iterator it = cLines.begin(); it!=cLines.end(); ++it)
        {
            if (!(*it).empty())
            {
                const int idx_l = (*it).front().ID_l();
                const int idx_r = (*it).back().ID_r();
                for (std::list<Line>::iterator it1 = cLines.begin(); it1!=cLines.end(); ++it1)
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
                            for (Line::iterator it2=(*it1).begin(); it2!=(*it1).end(); ++it2)
                                (*it2).flip();
                            (*it).splice((*it).end(),*it1);
                            changed = true;
                        }
                        else if (idx_l == (*it1).front().ID_l())
                        {
                            (*it1).reverse();
                            for (Line::iterator it2=(*it1).begin(); it2!=(*it1).end(); ++it2)
                                (*it2).flip();
                            (*it).splice((*it).begin(),*it1);
                            changed = true;
                        }
                    }
                }
            }
        }
    } while(changed);

    // Clear out empty placeholders for lines we connected up to others.
    for (std::list<Line>::iterator it = cLines.begin(); it!=cLines.end();)
    {
        if (it->empty())
            cLines.erase(it++);
        else
            it++;
    }

    // Dump out the final number of contours
    cerr << "  number of contour lines: " << cLines.size() << endl;

    return cLines;
}

static
void CheckSurfaceNormal(const Vector<Triangle>& eltTris,const FArrayBox& sfab,const IntVect& iv)
{
}

int
main (int   argc,
      char* argv[])
{
  amrex::Initialize(argc,argv);
  {
    if (argc < 2) {
      print_usage(argc,argv);
    }

    ParmParse pp;

    if (pp.contains("help")) {
      print_usage(argc,argv);
    }

    int verbose=0;
    pp.query("verbose",verbose);

    int doCollate=1;
    pp.query("collate",doCollate);

    std::string infile("");
    pp.get("infile", infile);
    if (infile.empty()) {
       Abort("Plotfile not specified, Use infile=");
    }

    // Load pltfile
    PlotFileData pf(infile);

    // Iso-surface definition
    Real isoVal = isoVal_DEF;
    pp.query("isoVal",isoVal);
    string isoCompName = isoCompName_DEF;
    pp.query("isoCompName",isoCompName);

    // Plotfile fields to read in
    Vector<int> pltComps;
    if (int nc = pp.countval("comps")) {
      pltComps.resize(nc);
      pp.getarr("comps",pltComps,0,nc);
    } else {
      int pltsComp = 0;
      pp.query("sComp",pltsComp);
      int pltnComp = 1;
      pp.query("nComp",pltnComp);
      AMREX_ASSERT(pltsComp+pltnComp <= pf.nComp());
      pltComps.resize(pltnComp);
      for (int i=0; i<pltnComp; ++i) {
        pltComps[i] = pltsComp + i;
      }
    }

    int isoComp = -1;
    const Vector<std::string>& pltNames = pf.varNames();
    Vector<string> varnames(pltComps.size()); // names of variables to get
    for (int i=0; i<varnames.size(); ++i) {
      if (pltComps[i]>=pltNames.size()) {
        Abort("At least one of the components requested is not in pltfile");
      }
      varnames[i] = pltNames[pltComps[i]];
      if (varnames[i]==isoCompName) isoComp = AMREX_SPACEDIM + i;
    }
    if (isoComp<AMREX_SPACEDIM) {
      Abort("isoCompName not in list of variables to read in");
    }

    const int nComp = pltComps.size();

    int finestLevel = pf.finestLevel();
    pp.query("finestLevel",finestLevel);
    AMREX_ASSERT(finestLevel <= pf.finestLevel());
    int Nlev = finestLevel + 1;
    Vector<BoxArray> grids(Nlev);
    for (int lev=0; lev<Nlev; ++lev)
    {
      grids[lev] = pf.boxArray(lev);
    }

    const int nodesPerElt = AMREX_SPACEDIM;
    bool rm_external_elements = true;
    pp.query("rm_external_elements",rm_external_elements);
    bool build_distance_function = false;
    pp.query("build_distance_function",build_distance_function);
    Vector<std::unique_ptr<MultiFab>> distance(Nlev);
    if (build_distance_function && AMREX_SPACEDIM==2) {
      Abort("Distance function not worked out for 2D yet");
    }

    Vector<int> nGrow(Nlev,1);
    Real dmax = pf.probSize()[0] / pf.probDomain(0).length(0);
    if (build_distance_function) {
      pp.query("dmax",dmax);
      Print() << "dmax: " << dmax << std::endl;
      for (int lev=0; lev<Nlev; ++lev) {
        Real dxL = pf.probSize()[lev] / pf.probDomain(lev).length(0);
        nGrow[lev] = int(dmax*(1.0000001) / dxL);
      }
    } else {
      nGrow[0] = 1;
      pp.query("nGrow",nGrow[0]);
      for (int lev=1; lev<Nlev; ++lev) {
        nGrow[lev] = nGrow[0];
      }
    }

    //-----------------------------------------------------------------
    // IO
    //-----------------------------------------------------------------
    const Real strt_time_io = ParallelDescriptor::second();

    Vector<Vector<MultiFab>> pfdata(Nlev);
    Vector<Geometry> geoms(Nlev);
    RealBox rb(&(pf.probLo()[0]), &(pf.probHi()[0]));
    int coord = pf.coordSys();
    Vector<int> is_per(AMREX_SPACEDIM,0);
    pp.queryarr("is_per",is_per,0,AMREX_SPACEDIM);
    Print() << "Periodicity assumed for this case: ";
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
      Print() << is_per[idim] << " ";
    }
    Print() << std::endl;

    for (int lev=0; lev<Nlev; ++lev) {
      pfdata[lev].resize(nComp);
      Print() << "Reading the plotfile data at level " << lev
              << "... (nComp,nGrow) = (" << nComp
              << "," << nGrow[lev] << ")"<< std::endl;
      for (int n=0; n<nComp; ++n) {
        Print() << "   var = " << varnames[n] << "..." << std::endl;
        pfdata[lev][n] = pf.get(lev,varnames[n]);
      }
      Print() << "...done reading the plotfile data at level " << lev << "..." << std::endl;
      geoms[lev] = Geometry(pf.probDomain(lev),&rb,coord,&(is_per[0]));
    }
    const Real end_time_io = ParallelDescriptor::second();
    Real io_time = end_time_io - strt_time_io;

    Real time = 0;
    PhysBCFunctNoOp bcFunc;
    PCInterp pci;

    //-----------------------------------------------------------------
    // Isosurface
    //-----------------------------------------------------------------
    const Real strt_time_surf = ParallelDescriptor::second();

    // Global isosurface data containers
    std::set<Node> nodeSet;
    std::set<Element> eltSet;
    int nodeCtr = 0;

    Vector<MultiFab> states(Nlev);
    Vector<DistributionMapping> dmaps(Nlev);

    for (int lev=0; lev<Nlev; ++lev) {

      // Get a domain box grown in periodic dir
      const Box& domain = geoms[lev].Domain();
      const Box& gpdomain = geoms[lev].growPeriodicDomain(nGrow[lev]);

      // Define data holders
      dmaps[lev] = DistributionMapping(grids[lev]);
      BoxArray gba = BoxArray(grids[lev]).grow(nGrow[lev]);

      states[lev].define(grids[lev],dmaps[lev],nComp+AMREX_SPACEDIM,nGrow[lev]);
      if (build_distance_function) {
        distance[lev].reset(new MultiFab(grids[lev],dmaps[lev],1,nGrow[lev]));
      }

      // Get cell center coordinates
      const auto geomdata = geoms[lev].data();
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (MFIter mfi(states[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
      {
        Box const& gbx = mfi.growntilebox();
        const auto& sta = states[lev].array(mfi);
        AMREX_PARALLEL_FOR_3D ( gbx, i, j, k,
        {
          const amrex::Real* plo = geomdata.ProbLo();
          const amrex::Real* dx  = geomdata.CellSize();
          AMREX_D_TERM(sta(i,j,k,0) = (i + 0.5)*dx[0] + plo[0];,
                       sta(i,j,k,1) = (j + 0.5)*dx[1] + plo[1];,
                       sta(i,j,k,2) = (k + 0.5)*dx[2] + plo[2]);
        });
      }

      states[lev].FillBoundary(0,AMREX_SPACEDIM,geoms[lev].periodicity()); // After this, bad data in periodic directions however

      if (lev > 0) {
        Vector<BCRec> bc(AMREX_SPACEDIM);
        int r = pf.refRatio(lev-1);
        IntVect ratio(AMREX_D_DECL(r,r,r));
        FillPatchTwoLevels(states[lev], time,
                           {&states[lev-1]},{time},
                           {&states[lev]},{time},
                           0, 0, AMREX_SPACEDIM, geoms[lev-1], geoms[lev],
                           bcFunc, 0, bcFunc, 0, ratio, &pci, bc, 0);
      }

      // In periodic direction, manually shift the coordinate back to its original location
      if (geoms[lev].isAnyPeriodic()) {
        Vector<IntVect> pshifts(27);

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(states[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
          Box const& bx = mfi.tilebox();
          auto& sfab = states[lev][mfi];
          geoms[lev].periodicShift(domain,bx,pshifts);
          for (const auto& iv : pshifts) {
            const auto shbox = bx + iv;
            auto pisect = domain & shbox;
            if (pisect.ok()) {
              pisect -= iv;
              for (int d=0; d<AMREX_SPACEDIM; ++d) {
                auto L = geoms[lev].ProbLength(d);
                if (iv[d] != 0) {
                  sfab.plus<amrex::RunOn::Device>((iv[d]>0 ? -L : L),pisect,d,1); // shift the location back
                }
              }
            }
          }
        }
      }

      // Use FillPatch to get properly filled ghost cells
      Print() << "FillPatching the grown structures at level " << lev << "..." << std::endl;
      MultiFab gstate(grids[lev],dmaps[lev],1,nGrow[lev]);
      gstate.setVal(-666);
      for (int n=0; n<nComp; ++n) {
        if (lev==0) {
          FillPatchSingleLevel(gstate,time,{&pfdata[lev][n]},{time},0,0,1,geoms[0],bcFunc,0);
        } else {
          BCRec bc;
          int r = pf.refRatio(lev-1);
          IntVect ratio(AMREX_D_DECL(r,r,r));
          FillPatchTwoLevels(gstate, time,
                             {&pfdata[lev-1][n]},{time},
                             {&pfdata[lev][n]},{time},
                             0, 0, 1, geoms[lev-1], geoms[lev],
                             bcFunc, 0, bcFunc, 0, ratio, &pci, {bc}, 0);
        }
        states[lev].ParallelCopy(gstate, 0, AMREX_SPACEDIM+n, 1, nGrow[lev], nGrow[lev], geoms[lev].periodicity());
      }
      Print() << "...done FillPatching the grown structures at level " << lev << "..." << std::endl;

      // Populate the list of elements of the isosurface
      for (MFIter mfi(states[lev]); mfi.isValid(); ++mfi) {

        const auto& sfab = states[lev][mfi];

        // Build a fine-covered mask using gstate
        auto& mask = gstate[mfi];
        const auto& gbox = mask.box();
        const auto g1box = grow(mfi.validbox(),1);

        mask.setVal<amrex::RunOn::Device>(1.0);

        if (lev<finestLevel && !build_distance_function) {
          const auto fratio = pf.refRatio(lev);
          const auto& fineBoxes = grids[lev+1];
          for (int i=0; i<fineBoxes.size(); ++i) {
            const auto cgFineBox = Box(fineBoxes[i]).coarsen(fratio);
            const auto isect = gbox & cgFineBox;
            if (isect.ok()) {
              mask.setVal<amrex::RunOn::Device>(-1.0,isect,0);
            }
            if (geoms[lev].isAnyPeriodic()) {
              Vector<IntVect> pshifts(27);
              geoms[lev].periodicShift(gbox,cgFineBox,pshifts);
              for (const auto& iv : pshifts) {
                auto shbox = cgFineBox + iv;
                const auto pisect = gbox & shbox;
                if (pisect.ok()) {
                  mask.setVal<amrex::RunOn::Device>(-1.0,pisect,0);
                }
              }
            }
          }
        }

        // For looping over quads/bricks, make set of "base" points
        Box loopBox(gbox & gpdomain);
        for (int i=0; i<AMREX_SPACEDIM; ++i) {
          loopBox.growHi(i,-1);
        }

        // Build list of <Edge,Point>, collecting over all the cells in the current box
        PMap vertCache;

#if AMREX_SPACEDIM==2
        SegList elements;
        for (IntVect iv=loopBox.smallEnd(); iv<=loopBox.bigEnd(); loopBox.next(iv))
        {
          auto eltSegs = Segmentise(sfab,mask,vertCache,iv,isoVal,isoComp);
          for (int i = 0; i < eltSegs.size(); i++) {
            elements.push_back(eltSegs[i]);
          }
        }
#else
        TriList elements;
        for (IntVect iv=loopBox.smallEnd(); iv<=loopBox.bigEnd(); loopBox.next(iv))
        {
          auto eltTris = Polygonise(sfab,mask,vertCache,iv,isoVal,isoComp);
          for (int i = 0; i < eltTris.size(); i++) {
            elements.push_back(eltTris[i]);
          }
          CheckSurfaceNormal(eltTris,sfab,iv);
        }
#endif

        if (build_distance_function) {
          if (elements.size()>0) {
#if AMREX_SPACEDIM==3
            std::vector<Vec3f> vertList;
            std::vector<Vec3ui> faceList;

            std::map<PMapIt,unsigned int,PMapItCompare> ptID;
            unsigned int id = 0;
            for (PMapIt it=vertCache.begin(); it!=vertCache.end(); ++it) {
              const auto& loc = it->second;
              vertList.push_back(Vec3f(loc[0],loc[1],loc[2]));
              ptID[it] = id++;
            }

            for (const auto& elt : elements) {
              faceList.push_back(Vec3ui(ptID[elt[0]],ptID[elt[1]],ptID[elt[2]]));
            }

            //const Box& vbox = grids[lev][mfi.index()];
            Array<Real,AMREX_SPACEDIM> plo = pf.probLo();
            GpuArray<Real,AMREX_SPACEDIM> dxf;
            for (int i=0; i<AMREX_SPACEDIM; ++i) {
              dxf[i] = pf.probSize()[i] / pf.probDomain(lev).length(i);
            }
            const Box& vbox = (*distance[lev])[mfi].box();
            Vec3f local_origin(plo[0] + vbox.smallEnd()[0]*dxf[0],
                               plo[1] + vbox.smallEnd()[1]*dxf[1],
                               plo[2] + vbox.smallEnd()[2]*dxf[2]);
            Array3f phi_grid;
            float dx = float(dxf[0]);
            make_level_set3(faceList, vertList, local_origin, dx,
                            vbox.length(0),vbox.length(1),vbox.length(2), phi_grid);

            vertList.clear();
            faceList.clear();
            ptID.clear();

            const auto& d = distance[lev]->array(mfi);
            const auto& sfaba = states[lev].array(mfi);
            const int* lo = vbox.loVect();
            const int* hi = vbox.hiVect();

            for (int k=lo[2]; k<=hi[2]; ++k) {
              int kL=k-lo[2];
              for (int j=lo[1]; j<=hi[1]; ++j) {
                int jL=j-lo[1];
                for (int i=lo[0]; i<=hi[0]; ++i) {
                  int iL=i-lo[0];
                  Real abs_d = std::min(dmax,Real(phi_grid(iL,jL,kL)));
                  int sgn = sfaba(i,j,k,isoComp) < isoVal ? -1 : +1;
                  d(i,j,k) = sgn * abs_d;
                }
              }
            }
#endif
          } else {
            const auto& vb = grids[lev][mfi.index()];
            const auto& iv = vb.smallEnd();
            (*distance[lev])[mfi].setVal<amrex::RunOn::Device>(sfab(iv,isoComp) < isoVal ? -dmax : + dmax);
          }
        }

        // If all nodes of element not in g1box, remove before merging set with master list
        if (rm_external_elements) {
          std::set<PMapIt,PMapItCompare> ptsToRm;
#if AMREX_SPACEDIM==2
          SegList eltsInside;
#else
          TriList eltsInside;
#endif
          for (const auto& it : elements) {
            bool eltValid = true;
            // Loop on PMap defining elements
            for (int i = 0; i < nodesPerElt; ++i) {
              const auto& edge = it[i]->first;
              if ( !(g1box.contains(edge.IV_l) && g1box.contains(edge.IV_r)) ) {
                ptsToRm.insert(it[i]);
                eltValid = false;
              }
            }
            if (eltValid) eltsInside.push_back(it);
          }
          for (const auto& it : ptsToRm) {
            vertCache.erase(it);
          }
          std::swap(elements,eltsInside);
          eltsInside.clear();
        }

        // Loop through the vertices in cache and append to the global list of nodes
        // while keeping track of the correspondance between the vertices in vertCache
        // and the global node indices
        std::map<PMapIt,std::set<Node>::iterator,PMapItCompare> PMI_N_map;

        for (PMapIt it=vertCache.begin(); it!=vertCache.end(); ++it) {
          // Build a Node from the vertex
          const Node n(it->second,nodeSet.size());
          // Add to global list if not already contained
          std::pair<std::set<Node>::iterator,bool> nsit = nodeSet.insert(n);
          // If the node was added, store the key just added to nodeSet
          // otherwise find which key this node was already at
          if (nsit.second) {
            PMI_N_map[it] = nsit.first;
          } else {
            PMI_N_map[it] = nodeSet.find(n);
          }
        }

        // Loop through the elements created above and add to global list
        // associating element vertices with node indices from the nodeSet.
        // (unless element are degenerate)
        std::vector<int> v(AMREX_SPACEDIM);
#if AMREX_SPACEDIM==2
        for (SegList::const_iterator it=elements.begin(); it!=elements.end(); ++it)
        {
          // For each vertices of the element, get the global node index
          for (int k=0; k<AMREX_SPACEDIM; ++k) {
            v[k] = PMI_N_map[(*it)[k]]->m_idx;
          }
          if (v[0] != v[1]) eltSet.insert(Element(v));
        }
#else
        for (TriList::const_iterator it=elements.begin(); it!=elements.end(); ++it)
        {
          // For each vertices of the element, get the global node index
          for (int k=0; k<AMREX_SPACEDIM; ++k) {
            v[k] = PMI_N_map[(*it)[k]]->m_idx;
          }
          const bool degenerate = (v[0]==v[1] || v[1]==v[2] || v[0]==v[2]);
          if (!degenerate) eltSet.insert(Element(v));
        }
#endif
      } // MFIter
    } // Loop on levels

    // Write distance function to disk if required
    ParallelDescriptor::Barrier();
    if (build_distance_function) {
      std::string outfile("distance");
      pp.query("outfile",outfile);
      Vector<int> levelSteps(Nlev);
      Vector<IntVect> refRatio(Nlev-1);
      Vector<const MultiFab*> ptrs(Nlev);
      for (int lev=0; lev<Nlev; ++lev) {
        ptrs[lev] = distance[lev].get();
        levelSteps[lev] = pf.levelStep(lev);
        if (lev<finestLevel) {
          int ir = pf.refRatio(lev);
          refRatio[lev] = IntVect(AMREX_D_DECL(ir,ir,ir));
        }
      }
      WriteMultiLevelPlotfile(outfile,Nlev,ptrs,{"distance"},geoms,time,levelSteps,refRatio);
      distance.clear();
    }

    // Simulate a sort of the node set to node ordering to be consistent with element pointers
    std::vector<std::set<Node>::iterator> sortedNodes(nodeSet.size());
    for (std::set<Node>::iterator it=nodeSet.begin(); it!=nodeSet.end(); ++it) {
      sortedNodes[it->m_idx] = it;
    }

    const Real end_time_surf = ParallelDescriptor::second();
    const Real surf_time = end_time_surf - strt_time_surf;

    Real surf_time_max, surf_time_min; surf_time_max = surf_time_min = surf_time;
    Real io_time_max, io_time_min; io_time_max = io_time_min = io_time;

    const int IOProc   = ParallelDescriptor::IOProcessorNumber();
    ParallelDescriptor::ReduceRealMax(surf_time_max,IOProc);
    ParallelDescriptor::ReduceRealMax(io_time_max,IOProc);
    ParallelDescriptor::ReduceRealMin(surf_time_min,IOProc);
    ParallelDescriptor::ReduceRealMin(io_time_min,IOProc);
    if (ParallelDescriptor::IOProcessor()) {
      std::cout << "Max Compute Surface time: " << surf_time_max << '\n';
      std::cout << "Min Compute Surface time: " << surf_time_min << '\n';
      std::cout << "Max I/O time: " << io_time_max << '\n';
      std::cout << "Min I/O time: " << io_time_min << '\n';
    }

    // Prepare floating point and integer data for MPI communication (make two arrays to pass around)
    const int nReal = (pltComps.size()+AMREX_SPACEDIM)*nodeSet.size();
    Vector<Real> nodeRaw(nReal);
    for (long i = 0; i < sortedNodes.size(); ++i)
    {
      const Real* vec = sortedNodes[i]->m_vec;
      const int N = sortedNodes[i]->m_size;
      for (int j = 0; j < N; ++j) {
        nodeRaw[i*N+j] = vec[j];
      }
    }

    int cnt = 0;
    Vector<int> eltRaw(nodesPerElt*eltSet.size());

    // Checks
    for (std::set<Element>::const_iterator it = eltSet.begin(); it != eltSet.end(); ++it) {
      const Element& e = *it;
      if ( AMREX_D_TERM(   e[0]>=sortedNodes.size() || e[0]<0,
                        || e[1]>=sortedNodes.size() || e[1]<0,
                        || e[2]>=sortedNodes.size() || e[2]<0)
        )
      {
        std::cerr << "***** " << ParallelDescriptor::MyProc() << " *** bad element" << '\n';
        std::cerr << AMREX_D_TERM(e[0], << ' ' << e[1], << ' ' << e[2]) << '\n';
        Abort("Bailing...");
      }
      const int N = e.size();
      AMREX_ASSERT(N==nodesPerElt);
      for (int j=0; j<N; ++j) {
        eltRaw[cnt*N+j] = e[j];
      }
      ++cnt;
    }

    // All relevant data now in "raw" arrays
    nodeSet.clear();
    eltSet.clear();
    sortedNodes.clear();

    // Communicate node and element info from all procs to IOProc
    if (doCollate) {
       Collate(nodeRaw,eltRaw,nComp+AMREX_SPACEDIM);
    }

    if (ParallelDescriptor::IOProcessor()) {
      // Uniquify nodes, and make elements consistent
      const Real strt_time_uniq = ParallelDescriptor::second();

      Vector<Real> newData(nComp+AMREX_SPACEDIM);
      int nNodes = nodeRaw.size()/(nComp+AMREX_SPACEDIM);

      nodeCtr = 0;
      Vector<std::set<Node>::iterator> nodeVec(nNodes);

      // Populate nodeSet and nodeVec
      AMREX_ASSERT(nodeSet.size()==0);

      if (verbose) {
        Print() << "Number of collated nodes: " << nNodes << std::endl;
        Print() << "  size of nodeSet iterator: " << sizeof(std::set<Node>::iterator) << std::endl;
      }

      // Re-use nodeSet to build final list of nodes
      for (int j=0; j<nNodes; ++j) {
        for (int k=0; k<nComp+AMREX_SPACEDIM; ++k) {
          newData[k] = nodeRaw[j*(nComp+AMREX_SPACEDIM)+k];
        }

        Node n(newData,nodeSet.size());
        std::pair<std::set<Node>::iterator,bool> nsit = nodeSet.insert(n);
        if (nsit.second) {
          nodeVec[j] = nsit.first;
        } else {
          nodeVec[j] = nodeSet.find(n);
        }
      }
      // Free up some memory
      nodeRaw.clear();

      if (verbose) {
        Print() << "  final node merge complete, total number nodes now " << nodeSet.size() << std::endl;
      }

      // simulate a sort of the node set to node ordering to be consistent with element pointers
      sortedNodes.resize(nodeSet.size());
      for (std::set<Node>::iterator it=nodeSet.begin(); it!=nodeSet.end(); ++it) {
        sortedNodes[it->m_idx] = it;
      }

      Vector<int> eltData(AMREX_SPACEDIM);
      int nElts = eltRaw.size()/AMREX_SPACEDIM;

      if (verbose) {
        Print() << "Number of collated elements: " << nElts << std::endl;
        Print() << "  size of Element: " << sizeof(Element) << std::endl;
      }

      AMREX_ASSERT(eltSet.size()==0);
      for (int j=0; j<nElts; ++j) {
        for (int k=0; k<AMREX_SPACEDIM; ++k) {
          eltData[k] = (nodeVec[eltRaw[j*AMREX_SPACEDIM+k]])->m_idx;
        }
        eltSet.insert(eltData);
      }
      nElts = eltSet.size();

      // Clear up some memory
      nodeVec.clear();
      eltRaw.clear();

      if (verbose) {
        Print() << "  final element renumbering complete, total number elements now " << eltSet.size() << std::endl;
      }
      const Real end_time_uniq = ParallelDescriptor::second();
      const Real uniq_time = end_time_uniq - strt_time_uniq;
      Print() << "Uniquify time: " << uniq_time << '\n';

      const Real strt_time_sout = ParallelDescriptor::second();
      bool writeSurf = true;
      pp.query("writeSurf",writeSurf);
      std::string surfFormat = "MEF";
      pp.query("surfFormat",surfFormat);

      if (writeSurf) {

        if (surfFormat == "MEF") {
          Print() << "...write surface in mef format (mef = Marcs element format)" << std::endl;

          Print() << "      (Nelts,Nnodes):(" << nElts << ", " << sortedNodes.size() << ")" << std::endl;

          // eltRaw presently contains elts prior to uniquefying, shrink and reload correctly
          eltRaw.resize(nElts*nodesPerElt);
          size_t icnt = 0;
          for (std::set<Element>::const_iterator eit=eltSet.begin(); eit!=eltSet.end(); ++eit) {
            for (int i=0; i<nodesPerElt; ++i) {
              eltRaw[icnt++] = (*eit)[i] + 1;
            }
          }
          // Get back some memory
          eltSet.clear();

          int nNodeSize = sortedNodes[0]->m_size;

          // If the surface is large, write data to disk/clear mem/read up into a single fab
          bool surface_is_large = false;
          pp.query("surface_is_large",surface_is_large);
          int chunk_size = 32768;
          pp.query("chunk_size",chunk_size);
          FABdata* tmpDataP;
          if (surface_is_large) {
            std::string tmpFile="isoTEMPFILE";
            pp.query("tmpFile",tmpFile);
            std::ofstream ost;
            ost.open(tmpFile.c_str());

            size_t Npts = sortedNodes.size();
            Box bigBox(IntVect::TheZeroVector(),IntVect(AMREX_D_DECL(Npts-1,0,0)));
            BoxArray little_boxes(bigBox); little_boxes.maxSize(chunk_size);
            FArrayBox little_fab;
            if (verbose) {
              Print() << "  staging vertex data to disk in " << little_boxes.size() << " chunks..." << std::endl;
            }

            for (int j=0; j<little_boxes.size(); ++j) {
              const Box& little_box = little_boxes[j];
              little_fab.resize(little_box,nNodeSize);
              Real* fdat = little_fab.dataPtr();
              icnt = 0;
              for (int i=little_box.smallEnd(0); i<=little_box.bigEnd(0); ++i) {
                const Real* vec = sortedNodes[i]->m_vec;
                for (int k=0; k<nNodeSize; ++k) {
                  fdat[icnt++] = vec[k];
                }
              }
              little_fab.writeOn(ost);
            }
            little_fab.clear();

            ost.close();
            if (verbose) {
              Print() << "  ... data staged." << std::endl;
            }

            // Clear out some memory now
            sortedNodes.clear();
            nodeSet.clear();

            // Now allocate final fabdata, and populate with file data
            if (verbose) {
              Print() << "  allocating final fab data structure (" << Npts*nNodeSize << " data elements) ..." << std::endl;
            }

            tmpDataP = new FABdata(Npts,nNodeSize);
            if (verbose) {
              Print() << "  .... final fab structure allocated" << std::endl;
            }
            FABdata& tmpData = *tmpDataP;
            std::ifstream ist;
            ist.open(tmpFile.c_str());

            if (verbose) {
              Print() << "  retrieving vertex data back from disk...";
            }
            Real* fdat = tmpData.fab.dataPtr();
            for (int j=0; j<little_boxes.size(); ++j)
            {
              little_fab.readFrom(ist);
              const Box& little_box = little_fab.box();
              size_t icntL=0;
              icnt=little_fab.box().smallEnd(0) * nNodeSize;
              const Real* fdatL = little_fab.dataPtr();
              for (int i=little_box.smallEnd(0); i<=little_box.bigEnd(0); ++i)
              {
                for (int k=0; k<nNodeSize; ++k)
                {
                  fdat[icnt++] = fdatL[icntL++];
                }
              }
              AMREX_ASSERT(icntL==little_boxes[j].numPts() * nNodeSize);
            }
            ist.close();
            if (verbose) {
              Print() << "  ... data retrieved" << std::endl;
            }
          } else {
            tmpDataP = new FABdata(sortedNodes.size(),nNodeSize);
            FABdata& tmpData = *tmpDataP;
            icnt = 0;
            Real* fdat = tmpData.fab.dataPtr();
            for (size_t i=0; i<sortedNodes.size(); ++i) {
              const Real* vec = sortedNodes[i]->m_vec;
              for (int k=0; k<nNodeSize; ++k) {
                fdat[icnt++] = vec[k];
              }
            }
          }

          const FABdata& tmpData = *tmpDataP;

          /*
            In 2D we can connect up the line segments in order to get a set of
            lines.  We do this and hack in an output function that creates a
            file that ParaView can read.  May need to split out one zone per
            file...not sure.  Also, ParaView doesn't seem to like variables
            with names like Y(XX), so we translate them to XX.
          */
#if AMREX_SPACEDIM==2
          list<Line> contours = MakeCLines(tmpData.getFab(),eltRaw,nodesPerElt);

#if 0
          std::string fileName = "MARC.dat";
          std::ofstream ofm(fileName.c_str());

          // Build variable name string
          std::string mvars;
          for (int j=0; j<varnames.size(); ++j) {
            mvars += " " + varnames[j];
          }
          mvars.erase(std::remove(mvars.begin(), mvars.end(), 'X'), mvars.end());
          mvars.erase(std::remove(mvars.begin(), mvars.end(), 'Y'), mvars.end());
          mvars.erase(std::remove(mvars.begin(), mvars.end(), '('), mvars.end());
          mvars.erase(std::remove(mvars.begin(), mvars.end(), ')'), mvars.end());
          mvars = "X Y Z" + mvars;
          ofm << "VARIABLES = " << mvars << '\n';

          for (std::list<Line>::iterator it = contours.begin(); it!=contours.end(); ++it)
          {
            const auto& iLine = *it;
            ofm << "ZONE I=1 J=" << iLine.size()+1 << " k=1 FORMAT=POINT\n";
            for (Line::const_iterator it1 = iLine.begin(); it1!=iLine.end(); ++it1)
            {
              int nodeIdx = it1->ID_l();
              const Real* vec = sortedNodes[nodeIdx]->m_vec;
              for (int n=0; n<nNodeSize; ++n)
              {
                ofm << vec[n] << " ";
                if (n==1) ofm << "0 ";
              }
              ofm << '\n';
            }
            const Real* vec = sortedNodes[iLine.back().ID_r()]->m_vec;
            for (int n=0; n<nNodeSize; ++n)
            {
              ofm << vec[n] << " ";
              if (n==1) ofm << "0 ";
            }
            ofm << '\n';
          }
#endif

          const auto& xComp = Segment::xComp;
          const auto& yComp = Segment::yComp;
          const int pComp = yComp+1; // Guessing here
          for (std::list<Line>::iterator it = contours.begin(); it!=contours.end(); ++it) {
            Array<Real,2> integral = {0, 0};
            const auto& iLine = *it;
            for (Line::const_iterator it1 = iLine.begin(); it1!=iLine.end(); ++it1) {
              const auto seg = *it1;
              const Real* p0 = sortedNodes[seg.ID_l()]->m_vec;
              const Real* p1 = sortedNodes[seg.ID_r()]->m_vec;

              const Real x0 = p0[xComp];
              const Real x1 = p1[xComp];
              const Real y0 = p0[yComp];
              const Real y1 = p1[yComp];
              const Real avgp0 = p0[pComp];
              const Real avgp1 = p1[pComp];

              Real len = std::sqrt(((x1 - x0)*(x1 - x0)) + (y1-y0)*(y1-y0));
              Array<Real,2> normal = {0, 0};
              if (len > 0) {
                normal = {(y0 - y1)/len, (x1 - x0)/len};
              }
              for (int i=0; i<2; ++i) {
                integral[i] += normal[i] * 0.5 * (avgp0 + avgp1) * len;
              }
            }
            Print() << "Integral: " << integral[0] << " " << integral[1] << std::endl;
          }
#endif


          // Build variable name string
          std::string vars = "X Y";
#if AMREX_SPACEDIM==3
          vars += " Z";
#endif
          for (int j=0; j<varnames.size(); ++j) {
            vars += " " + varnames[j];
          }

          // Build a label and a filename
          char buf[72];
          sprintf(buf,"%g",pf.time());
          string label(buf);

          sprintf(buf, "%g", isoVal);
          std::string outfile_base = infile+"_"+isoCompName+"_"+string(buf);
          pp.query("outfile_base",outfile_base);
          std::string outfile(outfile_base+".mef");

          Print() << "  Writing the file..." << std::endl;
          std::ofstream ofs;
          std::ostream* os =
            (outfile=="-" ? (std::ostream*)(&std::cout) : (std::ostream*)(&ofs) );
          if (outfile!="-")
            ofs.open(outfile.c_str(),std::ios::out|std::ios::trunc|std::ios::binary);
          (*os) << label << std::endl;
          (*os) << vars << std::endl;
          (*os) << nElts << " " << nodesPerElt << std::endl;
          tmpData.fab.writeOn(*os);
          if (verbose) {
            Print() << "  clearing temporary fab data" << std::endl;
          }
          delete tmpDataP;
          (*os).write((char*)eltRaw.dataPtr(),sizeof(int)*eltRaw.size());
          if (outfile!="-") {
            ofs.close();
          }
          Print() << "            ...done" << std::endl;
        } else if (surfFormat == "XDMF") {
          // Setup eltRaw
          eltRaw.resize(nElts*nodesPerElt);
          size_t icnt = 0;
          for (std::set<Element>::const_iterator eit=eltSet.begin(); eit!=eltSet.end(); ++eit) {
            for (int i=0; i<nodesPerElt; ++i) {
              eltRaw[icnt++] = (*eit)[i];
            }
          }
          // Get back some memory
          eltSet.clear();

          size_t binEltSize = sizeof(int)*eltRaw.size();

          // Output file
          char buf[72];
          sprintf(buf,"%g",pf.time());
          string label(buf);
          std::string outfile_base = infile+"_"+isoCompName+"_"+string(buf);
          pp.query("outfile_base",outfile_base);
          std::string outfile(outfile_base+".mesh");

          // Write XDMF file
          std::ofstream outfile_data(outfile_base+".xmf");
          if(outfile_data.is_open()) {
             outfile_data.precision(8);
             outfile_data << "<?xml version=\"1.0\"?>\n";
             outfile_data << "<Xdmf Version=\"3.0\" xmlns:xi=\"http://www.w3.org/2001/XInclude\">\n";
             outfile_data << "   <Domain>\n";
             outfile_data << "      <Grid Name=\"isoSurface\">\n";
             outfile_data << "      <Information Name=\"Variable\" Value=\""<< isoCompName <<"\"/>\n";
             outfile_data << "      <Information Name=\"IsoValue\" Value=\""<< isoVal <<"\"/>\n";
             outfile_data << "      <Time Value=\"" << pf.time() << "\"/>\n";
             // Mesh
             if (AMREX_SPACEDIM == 2) {
                outfile_data << "         <Topology TopologyType=\"Polyline\" NodesPerElement=\"2\" NumberOfElements=\""<<nElts<< "\">\n";
             } else {
                outfile_data << "         <Topology TopologyType=\"Triangle\" NumberOfElements=\""<<nElts<< "\">\n";
             }
             outfile_data << "            <DataItem Name=\"Conn\" Format=\"Binary\" DataType=\"Int\" Dimensions=\""
                          << AMREX_SPACEDIM*nElts << "\">\n";
             outfile_data << "               "<<outfile<<"\n";
             outfile_data << "            </DataItem>\n";
             outfile_data << "         </Topology>\n";
             if (AMREX_SPACEDIM == 2) {
                outfile_data << "         <Geometry GeometryType=\"XY\">\n";
             } else {
                outfile_data << "         <Geometry GeometryType=\"XYZ\">\n";
             }
             outfile_data << "            <DataItem Name=\"Coord\" Format=\"Binary\" Precision=\"8\" DataType=\"Float\" Seek=\""<< binEltSize << "\" Dimensions=\""
                          << AMREX_SPACEDIM*sortedNodes.size() <<"\">\n";
             outfile_data << "               "<<outfile<<"\n";
             outfile_data << "            </DataItem>\n";
             outfile_data << "         </Geometry>\n";
             binEltSize += AMREX_SPACEDIM*sortedNodes.size()*sizeof(Real);
             // Data
             for (size_t comp = 0; comp < nComp; ++comp) {
                outfile_data << "         <Attribute Name=\""<<varnames[comp]<<"\" AttributeType=\"Scalar\" Center=\"Node\">\n";
                outfile_data << "            <DataItem Format=\"Binary\" Precision=\"8\" DataType=\"Float\" Seek=\""<< binEltSize << "\" Dimensions=\""
                             << sortedNodes.size() <<"\">\n";
                outfile_data << "               "<<outfile<<"\n";
                outfile_data << "            </DataItem>\n";
                outfile_data << "         </Attribute>\n";
                binEltSize += sortedNodes.size()*sizeof(Real);
             }
             outfile_data << "      </Grid>\n";
             outfile_data << "   </Domain>\n";
             outfile_data << "</Xdmf>\n";
             outfile_data.close();

             // Mesh data
             std::ofstream ofs;
             std::ostream* os = (std::ostream*)(&ofs);
             ofs.open(outfile.c_str(),std::ios::out|std::ios::trunc|std::ios::binary);

             // Connectivity
             (*os).write((char*)eltRaw.dataPtr(),sizeof(int)*eltRaw.size());

             // Coordinates
             for (size_t i = 0; i < sortedNodes.size(); ++i) {
                const Real* vec = sortedNodes[i]->m_vec;
                (*os).write((char*)vec,sizeof(Real)*AMREX_SPACEDIM);
             }

             // Data
             for (size_t comp = 0; comp < nComp; ++comp) {
                for (size_t i = 0; i < sortedNodes.size(); ++i) {
                   const Real* vec = sortedNodes[i]->m_vec;
                   int offset = (AMREX_SPACEDIM+comp)*sizeof(Real);
                   (*os).write((char*)vec+offset,sizeof(Real));
                }
             }
             ofs.close();
          }
        }
      }

      const Real end_time_sout = ParallelDescriptor::second();
      const Real sout_time = end_time_sout - strt_time_sout;
      std::cout << "Surface output time: " << sout_time << '\n';

      // Compute area of isosurface
      bool computeArea = false;
      pp.query("computeArea",computeArea);
      if (computeArea && (AMREX_SPACEDIM==3))  {
        Real Area = 0;
        for (std::set<Element>::const_iterator it = eltSet.begin(); it != eltSet.end(); ++it) {
          const Element& elt = *it;
          if (elt.size()==3) {
            if (elt[0]>=sortedNodes.size() || elt[1]>=sortedNodes.size() || elt[2]>=sortedNodes.size()) {
              std::cerr << "Accessing node past end: " << elt[0] << ", " << elt[1] << ", " << elt[2] << std::endl;
            }

            const Real* p0 = sortedNodes[elt[0]]->m_vec;
            const Real* p1 = sortedNodes[elt[1]]->m_vec;
            const Real* p2 = sortedNodes[elt[2]]->m_vec;

            Area += 0.5*sqrt(
              pow(( p1[1] - p0[1])*(p2[2]-p0[2])
                  -(p1[2] - p0[2])*(p2[1]-p0[1]), 2)

              + pow(( p1[2] - p0[2])*(p2[0]-p0[0])
                    -(p1[0] - p0[0])*(p2[2]-p0[2]), 2)

              + pow(( p1[0] - p0[0])*(p2[1]-p0[1])
                    -(p1[1] - p0[1])*(p2[0]-p0[0]), 2) );
          }
        }
        Print() << "Total area = " << Area << '\n';
      }
    } // IOProc
  }
  amrex::Finalize();
  return 0;
}
