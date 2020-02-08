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

static long NmyTriangles = 0;

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

    if (area < areaEps)
        return;

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

int
main (int   argc,
      char* argv[])
{
  Initialize(argc,argv);
  {
    ParmParse pp;

    std::string infile; pp.get("infile",infile);
    bool dumpFab = false; pp.query("dumpFab",dumpFab);
    std::string fabFileBase="bin"; pp.query("fabFileBase",fabFileBase);

    std::ifstream ifs;
    ifs.open(infile.c_str());

    const std::string title = parseTitle(ifs);
    const std::vector<std::string> names = parseVarNames(ifs);
    const int nComp = names.size();

    size_t nElts;
    int MYLEN;
    ifs >> nElts;
    ifs >> MYLEN;
    if (ParallelDescriptor::IOProcessor())
      cerr << "...finished reading data header" << endl;

    vector<int> binComps;
    int nc;
    if (nc = pp.countval("binComps"))
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

    Real areaEps = 1.e-12; pp.query("areaEps",areaEps);


    // Read MEF data
    FArrayBox nodeFab;
    nodeFab.readFrom(ifs);
    Real* nodeData = nodeFab.dataPtr();
    int nPts = nodeFab.box().numPts();
    if (ParallelDescriptor::IOProcessor())
      cerr << "..." << nPts << " nodes read from data file (nComp=" << nComp << ")" << endl;

    // Rotate data to be accessible pointwise for triangle stuff below
    vector<vector<Real> > nodeVec(nPts);
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
    vector<vector<int> > eltVec(nElts);
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
        const vector<Real>& C = nodeVec[E[2]-1];
            
        vector<int> Abin = getBin(A,binComps,binLO,binMax);
        vector<int> Bbin = getBin(B,binComps,binLO,binMax);
        vector<int> Cbin = getBin(C,binComps,binLO,binMax);
            
        area += triangleArea(A,B,C);
            
        processTriangle(A,Abin,B,Bbin,C,Cbin,bins,binLO,binMax,areaEps,binComps,
                        condComp,condVal,condSgn,condApply);
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
          box = Box(IntVect::TheZeroVector(),IntVect(D_DECL(nBins[0]-1,0,0)));
        }
        else
        {
          box = Box(IntVect::TheZeroVector(),IntVect(D_DECL(nBins[0]-1,nBins[1]-1,0)));
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
        
      cerr << "Total area of this surface: " << area << " (sum of bins: " << binSum << ")" << endl;
      if (condApply)
        cerr << "   area outside condition: " << areaOutsideCondition
             << " (total: " << areaOutsideCondition + binSum << ")" << endl;
    }

    if (ParallelDescriptor::NProcs()>1) {
      if (ParallelDescriptor::IOProcessor())
        std::cerr << "Load balance: " << std::endl;
        
      for (int i=0; i<ParallelDescriptor::NProcs(); ++i)
      {
        if (i==ParallelDescriptor::MyProc())
          std::cerr << "    " << i << ": " << NmyTriangles << std::endl;
            
        ParallelDescriptor::Barrier();
      }
        
      ParallelDescriptor::ReduceLongSum(NmyTriangles); 
        
      if (ParallelDescriptor::IOProcessor())
        std::cerr << "  Total: " << NmyTriangles << std::endl;
    }
  }
  Finalize();
  return 0;
}
