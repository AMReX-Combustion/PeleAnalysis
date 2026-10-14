#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <set>
#include <list>
#include <map>
#include <cmath>
#include <algorithm>
using std::cout;
using std::cerr;
using std::endl;
using std::string;
using std::vector;
using std::set;
using std::map;
using std::list;
#ifndef WIN32
using std::pow;
#endif
#include <AMReX_REAL.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Array.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_Utility.H>

using namespace amrex;

static
std::vector<std::string> parseVarNames(std::istream& is)
{
    std::string line;
    std::getline(is,line);
    return Tokenize(line,std::string(" "));
}

static std::string parseTitle(std::istream& is)
{
    std::string line;
    std::getline(is,line);
    return line;
}

static
Real triangleArea(const vector<Real>& p0,
                  const vector<Real>& p1,
                  const vector<Real>& p2)
{
    // Note: assumes (x,y,z) are first 3 components
    return 0.5*sqrt(
        pow(  ( p1[1] - p0[1])*(p2[2]-p0[2])
              -(p1[2] - p0[2])*(p2[1]-p0[1]), 2)

        + pow(( p1[2] - p0[2])*(p2[0]-p0[0])
              -(p1[0] - p0[0])*(p2[2]-p0[2]), 2)

        + pow(( p1[0] - p0[0])*(p2[1]-p0[1])
              -(p1[1] - p0[1])*(p2[0]-p0[0]), 2));
}


// Length of a 2-node element. The 2D counterpart of triangleArea: for a contour
// written by isosurface in 2D, elements are segments and the measure that
// partitions among bins is arc length rather than area.
//
// Uses the leading AMREX_SPACEDIM components as the coordinates, which is how
// isosurface lays out a MEF node ("X Y [Z] <fields...>"). Note triangleArea above
// hard-codes components 0,1,2, so it must not be used on a 2D MEF -- component 2
// there is the first interpolated field, not z.
static
Real segmentLength(const vector<Real>& p0,
                   const vector<Real>& p1)
{
    Real sumSq = 0;
    for (int d=0; d<AMREX_SPACEDIM; ++d)
    {
        sumSq += (p1[d] - p0[d])*(p1[d] - p0[d]);
    }
    return std::sqrt(sumSq);
}


static
void orderNodes(vector<Real>& A, vector<int>& Abin,
                vector<Real>& B, vector<int>& Bbin,
                vector<Real>& C, vector<int>& Cbin, int binID)
{
    // Order big to small.
    if (Bbin[binID] > Abin[binID])
    {
        vector<Real> t = A;
        vector<int> tbin = Abin;
        A = B; Abin = Bbin;
        B = t; Bbin = tbin;
    }
    if (Cbin[binID] > Bbin[binID])
    {
        vector<Real> t = B;
        vector<int> tbin = Bbin;
        B = C; Bbin = Cbin;
        C = t; Cbin = tbin;
    }
    if (Bbin[binID] > Abin[binID])
    {
        vector<Real> t = A;
        vector<int> tbin = Abin;
        A = B; Abin = Bbin;
        B = t; Bbin = tbin;
    }
}

// Order a 2-node element so that A is in the higher bin for coordinate binID.
// The 2D counterpart of orderNodes.
static
void orderNodesSegment(vector<Real>& A, vector<int>& Abin,
                       vector<Real>& B, vector<int>& Bbin, int binID)
{
    if (Bbin[binID] > Abin[binID])
    {
        vector<Real> t = A;
        vector<int> tbin = Abin;
        A = B; Abin = Bbin;
        B = t; Bbin = tbin;
    }
}

// Find the single point D where segment AB crosses the bin boundary below A, and
// interpolate all node states to it.
//
// The 2D counterpart of findFG. A triangle straddling a boundary needs two cut
// points (and yields three sub-triangles); a segment needs exactly one, and
// yields two sub-segments, which is why this is so much simpler than the
// triangle case.
static
void findD(const vector<Real>& A,
           const vector<Real>& B,
           vector<Real>&       D,
           const vector<Real>& binLO,
           Real                binMax,
           int                 abin,
           int                 comp)
{
    // Assumes A[comp] > B[comp] (i.e. orderNodesSegment has been applied) and
    // that A and B share a component layout.
    Real fAB;
    if (abin >= binLO.size()) // A above the upper bin bound
    {
        fAB = (A[comp] - binMax)/(A[comp] - B[comp]);
    }
    else
    {
        // A is the higher node, so its bin index cannot be below the range here:
        // that would require B's to be lower still, and getBin() floors at -1, so
        // both would be -1 and the caller would have taken the "same bin" branch.
        AMREX_ALWAYS_ASSERT(abin >= 0);
        fAB = (A[comp] - binLO[abin])/(A[comp] - B[comp]);
    }

    AMREX_ALWAYS_ASSERT(fAB>=0 && fAB<=1);

    // Interpolate all states to the interface
    for (int i=0; i<A.size(); ++i)
    {
        D[i] = A[i] - fAB*(A[i] - B[i]);
    }
}

static
void findDE(const vector<Real>& A,
            const vector<Real>& B,
            const vector<Real>& C,
            vector<Real>&       D,
            vector<Real>&       E,
            const vector<Real>& binLO,
            Real                binMax,
            int                 abin,
            int                 comp)
{
    // Assume that A[comp]=B[comp] and A[comp]>C[comp], and that ABC all have the same component layout
    Real fAC, fBC;
    if (abin < 0) // then A and B below lower bin bound
    {
        fAC = (binLO[0] - A[comp])/(C[comp] - A[comp]);
        fBC = (binLO[0] - B[comp])/(C[comp] - B[comp]);
    }
    else if (abin >= binLO.size()) // then A and B above upper  bin bound
    {
        fAC = (A[comp] - binMax)/(A[comp]-C[comp]);
        fBC = (B[comp] - binMax)/(B[comp]-C[comp]);
    }
    else
    {
        fAC = (A[comp]-binLO[abin])/(A[comp]-C[comp]);
        fBC = (B[comp]-binLO[abin])/(B[comp]-C[comp]);
    }

    AMREX_ALWAYS_ASSERT(fAC>=0 && fAC<=1 && fBC>=0 && fBC<=1);

    // Interpolate all states to the interface
    for (int i=0; i<A.size(); ++i)
    {
        D[i] = A[i] - fAC*(A[i] - C[i]);
        E[i] = B[i] - fBC*(B[i] - C[i]);
    }
}

static
void findFG(const vector<Real>& A,
            const vector<Real>& B,
            const vector<Real>& C,
            vector<Real>&       F,
            vector<Real>&       G,
            const vector<Real>& binLO,
            Real                binMax,
            int                 abin,
            int                 comp)
{
    // Assume that A[comp]>B[comp] and A[comp]>C[comp], and that ABC all have the same component layout
    Real fAB, fAC;
    if (abin < 0) // then A below lower bin bound
    {
        fAB = (binLO[0] - A[comp])/(A[comp] - B[comp]);
        fAC = (binLO[0] - C[comp])/(A[comp] - C[comp]);
    }
    else if (abin >= binLO.size()) // then A above upper bin bound
    {
        fAB = (A[comp] - binMax)/(A[comp]-B[comp]);
        fAC = (A[comp] - binMax)/(A[comp]-C[comp]);
    }
    else
    {
        fAB = (A[comp]-binLO[abin])/(A[comp]-B[comp]);
        fAC = (A[comp]-binLO[abin])/(A[comp]-C[comp]);
    }
    AMREX_ALWAYS_ASSERT(fAB>=0 && fAB<=1 && fAC>=0 && fAC<=1);
    for (int i=0; i<A.size(); ++i)
    {
        F[i] = A[i] - fAB*(A[i] - B[i]);
        G[i] = A[i] - fAC*(A[i] - C[i]);
    }
}

static
vector<int>
getBin (const vector<Real>&          val,
        const vector<int>&           binComps,
        const vector<vector<Real> >& binLO,
        const vector<Real>&          binMax)
{
    const int nc = binComps.size();

    vector<int> retVal(nc,0);

    for (int j = 0; j < nc; ++j)
    {
        if (val[binComps[j]]<binLO[j][0])
        {
            retVal[j] = -1;
        }
        else if (val[binComps[j]]>binMax[j])
        {
            retVal[j] = binLO[j].size();
        }
        else
        {
            std::vector<Real>::const_iterator it =
                std::upper_bound(binLO[j].begin(), binLO[j].end(), val[binComps[j]]);

            --it;

            retVal[j] = it - binLO[j].begin();
        }
    }

    return retVal;
}

// This is the running sum of triangles not added because
//  they didn't satisfy the condition
static Real areaOutsideCondition = 0.;
static
bool satisfyCondition(const vector<Real>& A, const vector<Real>& B, const vector<Real>& C,
                      int condComp, Real condVal, int condSgn)
{
    if (condSgn > 0)
    {
        if (A[condComp] > condVal  &&  B[condComp] > condVal  &&  C[condComp] > condVal)
            return true;
    }
    else if (condSgn < 0)
    {
        if (A[condComp] < condVal  &&  B[condComp] < condVal  &&  C[condComp] < condVal)
            return true;
    }
    else
    {
        if (A[condComp] == condVal  &&  B[condComp] == condVal  &&  C[condComp] == condVal)
            return true;
    }

    return false;
}

// The 2D counterpart of satisfyCondition: keep a segment only if BOTH its nodes
// satisfy the condition.
static
bool satisfyConditionSegment(const vector<Real>& A, const vector<Real>& B,
                             int condComp, Real condVal, int condSgn)
{
    if (condSgn > 0)
    {
        if (A[condComp] > condVal  &&  B[condComp] > condVal)
            return true;
    }
    else if (condSgn < 0)
    {
        if (A[condComp] < condVal  &&  B[condComp] < condVal)
            return true;
    }
    else
    {
        if (A[condComp] == condVal  &&  B[condComp] == condVal)
            return true;
    }

    return false;
}

static long NmyTriangles = 0;
static long NmySegments = 0;

// Elements this rank actually binned: triangles in 3D, segments in 2D.
static long NmyElements() { return NmyTriangles + NmySegments; }

static
void processTriangle(const vector<Real>& Ai, const vector<int>& AbinI,
                     const vector<Real>& Bi, const vector<int>& BbinI,
                     const vector<Real>& Ci, const vector<int>& CbinI,
                     map<vector<int>,Real >& bins,
                     const vector<vector<Real> >& binLO,
                     const vector<Real>& binMax,
                     Real areaEps,
                     const vector<int>& binComps,
                     int condComp, Real condVal, int condSgn, bool condApply,
                     int binID=0)
{
    const Real area = triangleArea(Ai,Bi,Ci);

    if (area < areaEps) {
      return;
    }

    if (binID >= AbinI.size())
    {
        bool in_range = true;
        for (int i=0; i<AbinI.size(); ++i)
            if (AbinI[i]<0 || AbinI[i]>=binLO[i].size())
                in_range = false;

        if (in_range)
        {
            NmyTriangles++;

            if ( !(condApply) || satisfyCondition(Ai,Bi,Ci,condComp,condVal,condSgn))
            {
                bins[AbinI] += area;
            }
            else
            {
                areaOutsideCondition += area;
            }
        }
        return;
    }
    else if ( (AbinI[binID]==BbinI[binID]) && (BbinI[binID]==CbinI[binID]) )
    {
        processTriangle(Ai,AbinI,Bi,BbinI,Ci,CbinI,bins,binLO,binMax,areaEps,binComps,
                        condComp,condVal,condSgn,condApply,binID+1);
    }
    else
    {
        vector<int> Abin = AbinI;
        vector<int> Bbin = BbinI;
        vector<int> Cbin = CbinI;
        vector<Real> A = Ai;
        vector<Real> B = Bi;
        vector<Real> C = Ci;
        orderNodes(A,Abin,B,Bbin,C,Cbin,binID);
        if (Abin[binID] == Bbin[binID])
        {
            int nComp = A.size();
            vector<Real> D(nComp), E(nComp);
            findDE(A,B,C,D,E,binLO[binID],binMax[binID],Abin[binID],binComps[binID]);
            vector<int> Dbin = getBin(D,binComps,binLO,binMax);
            vector<int> Ebin = getBin(E,binComps,binLO,binMax);
            for (int i=0; i<=binID; ++i)
            {
                Dbin[i] = Abin[i];
                Ebin[i] = Dbin[i];
            }
            processTriangle(A,Abin,B,Bbin,E,Ebin,bins,binLO,binMax,areaEps,binComps,
                            condComp,condVal,condSgn,condApply,binID+1);
            processTriangle(A,Abin,E,Ebin,D,Dbin,bins,binLO,binMax,areaEps,binComps,
                            condComp,condVal,condSgn,condApply,binID+1);

            // Decrement bin for triangles on other side of interp line
            Dbin[binID] = Abin[binID] - 1;
            Ebin[binID] = Ebin[binID] - 1;
            processTriangle(D,Dbin,C,Cbin,E,Ebin,bins,binLO,binMax,areaEps,binComps,
                            condComp,condVal,condSgn,condApply,binID);
        }
        else
        {
            int nComp = A.size();
            vector<Real> F(nComp), G(nComp);
            findFG(A,B,C,F,G,binLO[binID],binMax[binID],Abin[binID],binComps[binID]);
            vector<int> Fbin = getBin(F,binComps,binLO,binMax);
            vector<int> Gbin = getBin(G,binComps,binLO,binMax);
            for (int i=0; i<=binID; ++i)
            {
                Fbin[i] = Abin[i];
                Gbin[i] = Fbin[i];
            }
            processTriangle(A,Abin,F,Fbin,G,Gbin,bins,binLO,binMax,areaEps,binComps,
                            condComp,condVal,condSgn,condApply,binID+1);

            // Decrement bin for triangles on other side of interp line
            Fbin[binID] = Abin[binID] - 1;
            Gbin[binID] = Abin[binID] - 1;
            processTriangle(F,Fbin,B,Bbin,C,Cbin,bins,binLO,binMax,areaEps,binComps,
                            condComp,condVal,condSgn,condApply,binID);
            processTriangle(F,Fbin,C,Cbin,G,Gbin,bins,binLO,binMax,areaEps,binComps,
                            condComp,condVal,condSgn,condApply,binID);
        }
    }
}

// Distribute a 2-node element's arc length into bins, clipping it at every bin
// boundary. The 2D counterpart of processTriangle, and it follows the same shape:
//
//   * recurse over the binning coordinates via binID;
//   * when both nodes share a bin for the current coordinate, move to the next;
//   * otherwise cut at the boundary and process each fragment.
//
// The cutting step is where it simplifies. A triangle straddling a boundary has
// two topologies (one node separated, or two) needing two cut points and yielding
// three sub-triangles; a segment has only one topology, one cut point D, and two
// sub-segments AD and DB. So there is no orderNodes-driven case split beyond
// putting A on the high side.
//
// As in the triangle version, the fragment on the far side of the cut is
// re-processed at the SAME binID with its bin index decremented, so that a
// segment spanning many bins is clipped repeatedly until each piece lies in one
// bin. Length is therefore partitioned exactly, not apportioned by an
// extent-overlap approximation.
static
void processSegment(const vector<Real>& Ai, const vector<int>& AbinI,
                    const vector<Real>& Bi, const vector<int>& BbinI,
                    map<vector<int>,Real >& bins,
                    const vector<vector<Real> >& binLO,
                    const vector<Real>& binMax,
                    Real lengthEps,
                    const vector<int>& binComps,
                    int condComp, Real condVal, int condSgn, bool condApply,
                    int binID=0)
{
    const Real length = segmentLength(Ai,Bi);

    if (length < lengthEps) {
      return;
    }

    if (binID >= AbinI.size())
    {
        bool in_range = true;
        for (int i=0; i<AbinI.size(); ++i)
            if (AbinI[i]<0 || AbinI[i]>=binLO[i].size())
                in_range = false;

        if (in_range)
        {
            NmySegments++;

            if ( !(condApply) || satisfyConditionSegment(Ai,Bi,condComp,condVal,condSgn))
            {
                bins[AbinI] += length;
            }
            else
            {
                areaOutsideCondition += length;
            }
        }
        return;
    }
    else if (AbinI[binID]==BbinI[binID])
    {
        processSegment(Ai,AbinI,Bi,BbinI,bins,binLO,binMax,lengthEps,binComps,
                       condComp,condVal,condSgn,condApply,binID+1);
    }
    else
    {
        vector<int> Abin = AbinI;
        vector<int> Bbin = BbinI;
        vector<Real> A = Ai;
        vector<Real> B = Bi;
        orderNodesSegment(A,Abin,B,Bbin,binID);

        int nComp = A.size();
        vector<Real> D(nComp);
        findD(A,B,D,binLO[binID],binMax[binID],Abin[binID],binComps[binID]);
        vector<int> Dbin = getBin(D,binComps,binLO,binMax);
        for (int i=0; i<=binID; ++i)
        {
            Dbin[i] = Abin[i];
        }

        // A-side fragment: lies wholly inside A's bin for this coordinate, so
        // move on to the next binning coordinate.
        processSegment(A,Abin,D,Dbin,bins,binLO,binMax,lengthEps,binComps,
                       condComp,condVal,condSgn,condApply,binID+1);

        // Far-side fragment: D sits exactly on the boundary, so getBin would put
        // it back in A's bin. Decrement it and re-process at this same binID, in
        // case D..B crosses further boundaries.
        Dbin[binID] = Abin[binID] - 1;
        processSegment(D,Dbin,B,Bbin,bins,binLO,binMax,lengthEps,binComps,
                       condComp,condVal,condSgn,condApply,binID);
    }
}

int
main (int   argc,
      char* argv[])
{
  Initialize(argc,argv);
  {
    ParmParse pp;

    bool dumpFab = false; pp.query("dumpFab",dumpFab);
    std::string fabFileBase="bin"; pp.query("fabFileBase",fabFileBase);

    vector<vector<Real> > nodeVec;
    vector<vector<int> > eltVec;

#define DEBUG_TEST
#undef DEBUG_TEST
#ifdef DEBUG_TEST
    int nPts = 4;
    int nElts = 2;
    int nComp = 3;
    nodeVec.resize(nPts);
    nodeVec[0] = {-.1, -.1, 0};
    nodeVec[1] = {1.1, -.1, 0};
    nodeVec[2] = {1.1, 1.1, 0};
    nodeVec[3] = {-.1, 1.1, 0};
    eltVec.resize(nElts);
    eltVec[0] = {1, 2, 3};
    eltVec[1] = {1, 3, 4};
    // run with: binComps=0 1 binMin=0 0  binMax=1 1 nBins=1 1
    // Total area of this surface: 1.44 (sum of bins: 1)
#else
    std::string infile; pp.get("infile",infile);
    std::ifstream ifs;
    ifs.open(infile.c_str());

    const std::string title = parseTitle(ifs);
    const std::vector<std::string> names = parseVarNames(ifs);
    int nComp = names.size();

    size_t nElts;
    int MYLEN;
    ifs >> nElts;
    ifs >> MYLEN;

    // 2 nodes per element is a 2D contour (arc length), 3 is a 3D surface
    // (area). isosurface writes nodesPerElt = AMREX_SPACEDIM, so a mismatch means
    // the MEF was produced by a build of a different dimensionality -- in which
    // case the coordinate columns would also be miscounted, silently, since the
    // leading AMREX_SPACEDIM components are taken as (x,y[,z]).
    if (MYLEN != 2 && MYLEN != 3)
      Abort("binMEF supports 2 nodes per element (2D contours) or 3 (3D surfaces)");
    if (MYLEN != AMREX_SPACEDIM)
      Abort("This MEF has " + std::to_string(MYLEN) + " nodes per element but "
            "binMEF was built for DIM=" + std::to_string(AMREX_SPACEDIM)
            + ". Rebuild with DIM=" + std::to_string(MYLEN) + ".");

    if (ParallelDescriptor::IOProcessor())
      cerr << "...finished reading data header" << endl;

    // Read MEF data
    FArrayBox nodeFab;
    nodeFab.readFrom(ifs);
    Real* nodeData = nodeFab.dataPtr();
    int nPts = nodeFab.box().numPts();
    if (ParallelDescriptor::IOProcessor())
      cerr << "..." << nPts << " nodes read from data file (nComp=" << nComp << ")" << endl;

    // Rotate data to be accessible pointwise for triangle stuff below
    nodeVec.resize(nPts);
    for (int i=0; i<nPts; ++i)
    {
      nodeVec[i].resize(nComp);
      for (int n=0; n<nComp; ++n)
        nodeVec[i][n] = *nodeData++;
    }
    nodeFab.clear();

    Vector<int> connData(nElts*MYLEN,0);
    ifs.read((char*)connData.dataPtr(),sizeof(int)*connData.size());
    if (ParallelDescriptor::IOProcessor())
      cerr << "..." << nElts << " elements read from data file" << endl;

    // Rearrange elt data into handy format (this sure takes a while...)
    eltVec.resize(nElts);
    int cnt = 0;
    for (int i=0; i<nElts; ++i)
    {
      eltVec[i].resize(MYLEN);
      for (int n=0; n<MYLEN; ++n)
        eltVec[i][n]=connData[cnt++];
    }
    connData.clear();

    if (ParallelDescriptor::IOProcessor())
      cerr << "...finished reading data" << endl;
#endif

    vector<int> binComps;
    int nc;
    if ((nc = pp.countval("binComps")))
    {
      binComps.resize(nc);
      pp.getarr("binComps",binComps,0,nc);
      for (int i=0; i<nc; ++i)
        if (binComps[i]>=nComp)
          Abort("At least one element in binComps out of range");
    }
    else
      Abort("Need to specify binComps array");

    int nmin = pp.countval("binMin");
    vector<Real> binMin(nc);
    if (nmin==nc)
    {
      binMin.resize(nc);
      pp.getarr("binMin",binMin,0,nc);
    }
    else
      Abort("Number of binMin components must match number of binComps components");

    int nmax = pp.countval("binMax");
    vector<Real> binMax(nc);
    if (nmax==nc)
    {
      binMax.resize(nc);
      pp.getarr("binMax",binMax,0,nc);
    }
    else
      Abort("Number of binMax components must match number of binComps components");

    int nsize = pp.countval("nBins");
    vector<int> nBins(nc);
    if (nsize==nc)
    {
      nBins.resize(nc);
      pp.getarr("nBins",nBins,0,nc);
    }
    else
      Abort("Number of nBins components must match number of binComps components");

    // Conditional binning?
    //   condApply == true?  Do it
    //   condSgn=(-,0,+) --> op=(<,==,>)
    //   condition:
    //    area added if ALL vertices, v,  satisfy v[condComp] op condVal
    bool condApply = false; pp.query("condApply",condApply);
    int condComp = 0;
    Real condVal = 0.;
    int condSgn = 0;
    if (condApply)
    {
      // Force the user to supply conditional parameters
      pp.get("condComp",condComp);
      pp.get("condVal",condVal);
      pp.get("condSgn",condSgn);
    }

    vector<Real> dBin(nc);
    for (int i=0; i<nc; ++i)
      dBin[i] = (binMax[i]-binMin[i])/nBins[i];


    bool dumpBins = false; pp.query("dumpBins",dumpBins);
    vector<vector<Real> > binLO(nc);
    for (int j=0; j<nc; ++j)
    {
      BL_ASSERT(nBins[j]>0);
      binLO[j].resize(nBins[j]);
      for (int i=0; i<nBins[j]; ++i)
        binLO[j][i] = binMin[j]+i*dBin[j];

      if (dumpBins && ParallelDescriptor::IOProcessor())
      {
        cout << "bin: " << binComps[j] << " bounds: " << endl;
        for (int i=0; i<nBins[j]; ++i)
        {
          Real lo = binLO[j][i];
          Real hi = (i==nBins[j]-1 ?  binMax[j] : binLO[j][i+1]);
          cout << "         bin: [" << lo << "," << hi << "]" << endl;
        }
        cout << endl;
      }
    }

    Real areaEps = 1.e-20; pp.query("areaEps",areaEps);



    //
    // Let each CPU process a non-intersecting group of the triangles.
    //
    const int MyProc = ParallelDescriptor::MyProc();
    const int NProcs = ParallelDescriptor::NProcs();
    const int IOProc = ParallelDescriptor::IOProcessorNumber();
    const int M      = nElts / NProcs;
    const int lo     =  M * MyProc;
    const int hi     = (MyProc == NProcs-1) ? nElts-1 : (MyProc+1)*M-1;

    Real area = 0;
    map<vector<int>,Real > bins;

    int idx = 0;
    for (int i=0; i<nElts; ++i, ++idx)
    {
      if (idx >= lo && idx <= hi)
      {
        const vector<int>&  E = eltVec[i];
        const vector<Real>& A = nodeVec[E[0]-1];
        const vector<Real>& B = nodeVec[E[1]-1];

        vector<int> Abin = getBin(A,binComps,binLO,binMax);
        vector<int> Bbin = getBin(B,binComps,binLO,binMax);

        if (MYLEN == 2)
        {
          // 2D contour: elements are segments and the measure is arc length.
          area += segmentLength(A,B);

          processSegment(A,Abin,B,Bbin,bins,binLO,binMax,areaEps,binComps,
                         condComp,condVal,condSgn,condApply);
        }
        else
        {
          const vector<Real>& C = nodeVec[E[2]-1];
          vector<int> Cbin = getBin(C,binComps,binLO,binMax);

          area += triangleArea(A,B,C);

          processTriangle(A,Abin,B,Bbin,C,Cbin,bins,binLO,binMax,areaEps,binComps,
                          condComp,condVal,condSgn,condApply);
        }
      }
    }
    //
    // Communicate results back to IOProc.
    //
    //std::cerr << ParallelDescriptor::MyProc() << ": finished processing triangles." << std::endl;
    ParallelDescriptor::ReduceRealSum(area, IOProc);

    vector<int>  binIdx(bins.size()*nc);
    vector<Real> binDat(bins.size());
    int icnt = 0;
    for (std::map<vector<int>,Real >::const_iterator it=bins.begin(); it!=bins.end(); ++it, ++icnt)
    {
      const vector<int>& k=it->first;
      for (int j=0; j<nc; ++j)
        binIdx[icnt*nc+j] = k[j];
      binDat[icnt] = it->second;
    }
    //std::cerr << ParallelDescriptor::MyProc() << ": building maps, size = " << bins.size() << std::endl;
    //
    // Now get data to IOProc.
    //
    std::vector<int> pSizes = ParallelDescriptor::Gather(int(bins.size()),IOProc);

    for (int i=0; i<ParallelDescriptor::NProcs(); ++i)
    {
      if (ParallelDescriptor::IOProcessor())
      {
        if (i!=0)
        {
          std::vector<int> tmpI(pSizes[i]*nc); ParallelDescriptor::Recv(tmpI,i,101);
          std::vector<Real> tmpR(pSizes[i]);   ParallelDescriptor::Recv(tmpR,i,102);
          //
          // Add to map on IOProc.
          //
          vector<int> itmp(nc);
          for (int j=0; j<pSizes[i]; ++j)
          {
            for (int k=0; k<nc; ++k)
              itmp[k] = tmpI[j*nc+k];
            bins[itmp] += tmpR[j];
          }
        }

        //std::cerr << "Bin data received from proc: " << i << std::endl;
      }
      else if (ParallelDescriptor::MyProc()==i)
      {
        ParallelDescriptor::Send(binIdx,IOProc,101);
        ParallelDescriptor::Send(binDat,IOProc,102);
      }

      ParallelDescriptor::Barrier();
    }

    if (ParallelDescriptor::IOProcessor())
    {
      cerr << "number of nonempty bins: " << bins.size() << endl;

      // Sum bins
      Real binSum = 0;
      for (std::map<vector<int>,Real >::const_iterator it=bins.begin(); it!=bins.end(); ++it)
        binSum += it->second;


      //
      // Dump binned data.
      //
      if (dumpFab && nc<=2)
      {
        string outFabFile = fabFileBase + ".fab";

        Box box;
        if (nc==1)
        {
          box = Box(IntVect::TheZeroVector(),IntVect(AMREX_D_DECL(nBins[0]-1,0,0)));
        }
        else
        {
          box = Box(IntVect::TheZeroVector(),IntVect(AMREX_D_DECL(nBins[0]-1,nBins[1]-1,0)));
        }

        FArrayBox outFab(box,1);
        outFab.setVal(0.0);
        for (std::map<vector<int>,Real >::const_iterator it=bins.begin(); it!=bins.end(); ++it)
        {
          IntVect iv = IntVect::TheZeroVector();
          const vector<int>& myBins = it->first;
          for (int j=0; j<myBins.size(); ++j)
          {
            iv[j] = myBins[j];
          }
          outFab(iv,0) = it->second;
        }


        bool normalize = false; pp.query("normalize",normalize);
        if (normalize && dumpFab && nc<=2) {
          outFab.mult(1./binSum);
        }

        std::ofstream ofs;
        ofs.open(outFabFile.c_str(),std::ios::out|std::ios::trunc|std::ios::binary);
        outFab.writeOn(ofs);
        ofs.close();
      }
      else
      {
        for (std::map<vector<int>,Real >::const_iterator it=bins.begin(); it!=bins.end(); ++it)
        {
          const vector<int>& idxs = it->first;
          for (int j=0; j<idxs.size(); ++j)
          {
#if 0
            Real lo = binLO[j][idxs[j]];
            Real hi = (idxs[j]==nBins[j]-1 ?  binMax[j] : binLO[j][idxs[j]+1]);
            cout << "[" << lo << "," << hi << "]: ";
#else
            Real lo = binLO[j][idxs[j]];
            Real hi = (idxs[j]==nBins[j]-1 ?  binMax[j] : binLO[j][idxs[j]+1]);
            cout << 0.5*(lo+hi) << " ";
#endif
          }
          cout << it->second << endl;
        }
      }

      // "area" for a 3D surface, "length" for a 2D contour. The surrounding
      // format is deliberately unchanged, since it is what downstream scripts
      // parse to check that the bins account for the whole surface.
#if AMREX_SPACEDIM==2
      const char* measureLabel = "length";
#else
      const char* measureLabel = "area";
#endif
      cerr << "Total " << measureLabel << " of this surface: " << area
           << " (sum of bins: " << binSum << ")" << endl;
      if (condApply)
        cerr << "   " << measureLabel << " outside condition: " << areaOutsideCondition
             << " (total: " << areaOutsideCondition + binSum << ")" << endl;
    }

    if (ParallelDescriptor::NProcs()>1) {
      if (ParallelDescriptor::IOProcessor())
        std::cerr << "Load balance: " << std::endl;

      for (int i=0; i<ParallelDescriptor::NProcs(); ++i)
      {
        if (i==ParallelDescriptor::MyProc())
          std::cerr << "    " << i << ": " << NmyElements() << std::endl;

        ParallelDescriptor::Barrier();
      }

      ParallelDescriptor::ReduceLongSum(NmyTriangles);
      ParallelDescriptor::ReduceLongSum(NmySegments);

      if (ParallelDescriptor::IOProcessor())
        std::cerr << "  Total: " << NmyElements() << std::endl;
    }
  }
  Finalize();
  return 0;
}
