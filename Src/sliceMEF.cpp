#include "winstd.H"

#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <map>
#include "ParmParse.H"
#include "MultiFab.H"
#include "VisMF.H"
#include "ParallelDescriptor.H"
#include "DataServices.H"
#include "Geometry.H"
#include "Utility.H"
using std::vector;
using std::map;
using std::set;
using std::list;
using std::pair;
using std::string;
using std::cerr;
using std::endl;
using std::cout;
using std::ofstream;

static Real epsilon_DEF = 1.e-8;


void
read_iso(const std::string& infile,
         FArrayBox&         nodes,
         Array<int>&        faceData,
         int&               nElts,
         vector<string>&    names,
         string&            label);

void
write_iso(const std::string&    outfile,
          const FArrayBox&      nodes,
          const Array<int>&     faceData,
          int                   nElts,
          const vector<string>& names,
          const string&         label);

vector<std::string>
Tokenize (const std::string& instr, const std::string& separators)
{
    vector<char*> ptr;
    //
    // Make copy of line that we can modify.
    //
    char* line = new char[instr.size()+1];

    (void) strcpy(line, instr.c_str());

    char* token = 0;

    if (!((token = strtok(line, separators.c_str())) == 0))
    {
        ptr.push_back(token);
        while (!((token = strtok(0, separators.c_str())) == 0))
            ptr.push_back(token);
    }

    vector<std::string> tokens(ptr.size());

    for (int i = 1; i < ptr.size(); i++)
    {
        char* p = ptr[i];

        while (strchr(separators.c_str(), *(p-1)) != 0)
            *--p = 0;
    }

    for (int i = 0; i < ptr.size(); i++)
        tokens[i] = ptr[i];

    delete line;

    return tokens;
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

static string rootName(const string inStr)
{
#ifdef WIN32
    const string dirSep("\\");
#else
    const string dirSep("/");
#endif
    vector<string> res = Tokenize(inStr,dirSep);
    res = Tokenize(res[res.size()-1],string("."));
    string result = res[0];
    for (int i=1; i<res.size()-1; ++i)
        result = result + string(".") + res[i];
    return result;
}

typedef Array<Real> Point;
typedef list<Point> PointList;
typedef vector<Point> PointVec;
typedef pair<int,int> Edge;
typedef map<Edge, pair<Point,int> > PMap;
typedef PMap::iterator PMapIt;
static PMap vertCache;
vector<PMapIt> vertVec;




struct Segment
{
    Segment() : p(2), mLength(-1) {}
    Real Length ();
    const PMapIt& operator[] (int n) const { return p[n]; }
    PMapIt& operator[] (int n) { return p[n]; }
    int ID_l () const {return (*p[0]).second.second;}
    int ID_r () const {return (*p[1]).second.second;}
    void flip ();
    static int xComp,yComp;
    Array<PMapIt> p;
    
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

PMapIt VertexInterp(Real isoVal,int isoComp,const vector<Point>& pts,int p1,int p2);

Array<Segment> Segmentise(const vector<Point>& pts,
                          const vector<int>&   elt,
                          Real                 isoVal,
                          int                  isoComp);


SegList::iterator FindMySeg(SegList& segs, int idx);

int Nreuse=0;

int
main (int   argc,
      char* argv[])
{
    BoxLib::Initialize(argc,argv);

    ParmParse pp;

    int nElts;
    string infile; pp.get("infile",infile);

    FArrayBox nodes;
    Array<int> faceData;
    vector<string> names;
    string label;
    read_iso(infile,nodes,faceData,nElts,names,label);
    int nodesPerElt = faceData.size() / nElts;
    BL_ASSERT(nodesPerElt*nElts == faceData.size());

    nElts = faceData.size() / nodesPerElt;
    int nComp = nodes.nComp();
    BL_ASSERT(nElts*nodesPerElt == faceData.size());

    int dir=0; pp.query("dir",dir);
    Array<Real> loc(1,0);
    if (int nloc = pp.countval("locs"))
    {
        loc.resize(nloc);
        pp.getarr("locs",loc,0,nloc);
    }

    // Put data into more convenient structure (so I can steal code from elsewhere for this)
    Array<const Real*> dat(nComp);
    for (int i=0; i<dat.size(); ++i)
        dat[i] = nodes.dataPtr(i);

    Real* nodeData = nodes.dataPtr();
    const int nNodes = nodes.box().length(0);
    vector<Point> GridPts(nNodes);
    for (int i=0; i<nNodes; ++i) {
        GridPts[i].resize(nComp);
        for (int j=0; j<nComp; ++j) {
            GridPts[i][j] = dat[j][i];
        }
    }

    for (int k=0; k<loc.size(); ++k)
    {        
        SegList segments;
        vector<int> elt(3);
        for (int i=0; i<nElts; ++i)
        {
            for (int j=0; j<3; ++j)
            {
                elt[j] = faceData[i*nodesPerElt + j ] - 1;
            }

            Array<Segment> eltSegs = Segmentise(GridPts,elt,loc[k],dir);

            for (int j=0; j<eltSegs.size(); ++j)
            {
                segments.push_back(eltSegs[j]);
            }
        }
        cerr << names[dir] << "=" << loc[k] << " slice computed." << endl;

        // Number the isosurface vertices
        int cnt=0;
        for (PMapIt it=vertCache.begin(); it!=vertCache.end(); ++it)
            (*it).second.second = cnt++;

        // Build a reverse vector into the vertCache
        vertVec.resize(vertCache.size());
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
                SegList::iterator segIt = FindMySeg(segList,idx);
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

        string outfileRoot = rootName(infile);

        char buf[72];
        sprintf(buf, "%g", std::abs(loc[k])); 
        string locStr = (loc[k]<0 ? "n" : (loc[k] > 0 ? "p" : "" ) ) + string(buf);
        outfileRoot += "_"+names[dir]+"_"+locStr;
        bool write_tec = true; pp.query("write_tec",write_tec);
        if (write_tec)
        {
            ofstream os;
            string outfile=outfileRoot + ".dat";
            cerr << "Writing contour (" << vertCache.size() << " nodes, "
                 << segments.size() << " segments, " << cLines.size() << " lines" << ") to " << outfile << endl;
        
            os.open(outfile.c_str(),std::ios::out);
            os << "VARIABLES = ";
            for (int i=0; i<names.size(); ++i)
                os << "\"" << names[i] << "\" ";
            os << endl;
        
            int Ccnt = 0;
            for (std::list<SegList>::iterator it = cLines.begin(); it!=cLines.end(); it++)
            {
                char buf[72];
                sprintf(buf, "%g", loc[k]); 
                string zoneName = rootName(infile)+"_"+names[dir]+"_"+string(buf);
                sprintf(buf, "%d", Ccnt++);
                zoneName += "_"+string(buf);
            
                os << "ZONE T=\"" << zoneName << "\", I=" << (*it).size() + 1 << endl;
            
                int lcnt=0;
                for (SegList::const_iterator its=(*it).begin(); its!=(*it).end(); ++its)
                {
                    if (its==(*it).begin())
                    {
                        const Point& p = (*(*its)[0]).second.first;
                        for (int i=0; i<p.size(); ++i)
                            os << p[i] << " ";
                        os << endl;
                    }
                    const Point& p = (*(*its)[1]).second.first;  
                    for (int i=0; i<p.size(); ++i)
                        os << p[i] << " ";
                    os << endl;
                    lcnt++;
                }
            }
        }

        bool write_mef = true; pp.query("write_mef",write_mef);
        if (write_mef)
        {
            int Nslice_nodes = vertVec.size();
            const Box slice_box(IntVect::TheZeroVector(),Nslice_nodes*BoxLib::BASISV(0));
            FArrayBox slice_nodes(slice_box,nComp);
            Array<Real*> slice_dat(nComp);
            for (int i=0; i<slice_dat.size(); ++i)
                slice_dat[i] = slice_nodes.dataPtr(i);        

            for (int i=0; i<Nslice_nodes; ++i)
            {
                const Point& p = vertVec[i]->second.first;
                for (int j=0; j<nComp; ++j)
                    slice_dat[j][i] = p[j];
            }

            int Nsegs = segments.size();
            Array<int> slice_segs(2*Nsegs);
            cnt = 0;
            for (SegList::const_iterator it=segments.begin(); it!=segments.end(); ++it)
            {
                slice_segs[cnt++] = it->ID_l() + 1;
                slice_segs[cnt++] = it->ID_r() + 1;
            }
            string outfile=outfileRoot + ".mef";
            write_iso(outfile,slice_nodes,slice_segs,Nsegs,names,label);
        }

        vertCache.clear();
    }

    BoxLib::Finalize();
    return 0;
}

void
write_iso(const std::string&    outfile,
          const FArrayBox&      nodes,
          const Array<int>&     faceData,
          int                   nElts,
          const vector<string>& names,
          const string&         label)
{
    // Rotate data to vary quickest on component
    int nCompSurf = nodes.nComp();
    int nNodes = nodes.box().numPts();

    FArrayBox tnodes(nodes.box(),nCompSurf);
    const Real** np = new const Real*[nCompSurf];
    for (int j=0; j<nCompSurf; ++j)
        np[j] = nodes.dataPtr(j);
    Real* ndat = tnodes.dataPtr();
    for (int i=0; i<nNodes; ++i)
    {
        for (int j=0; j<nCompSurf; ++j)
        {
            ndat[j] = np[j][i];
        }
        ndat += nCompSurf;
    }
    delete [] np;

    std::ofstream ofs;
    ofs.open(outfile.c_str(),std::ios::out|std::ios::trunc|std::ios::binary);
    ofs << label << endl;
    for (int i=0; i<nCompSurf; ++i)
    {
        ofs << names[i];
        if (i < nCompSurf-1)
            ofs << " ";
        else
            ofs << std::endl;
    }
    int nodesPerElt = faceData.size() / nElts;
    ofs << nElts << " " << nodesPerElt << endl;
    tnodes.writeOn(ofs);
    ofs.write((char*)faceData.dataPtr(),sizeof(int)*faceData.size());
    ofs.close();
}

void
read_iso(const std::string& infile,
         FArrayBox&         nodes,
         Array<int>&        faceData,
         int&               nElts,
         vector<string>&    names,
         string&            label)
{
    std::ifstream ifs;
    ifs.open(infile.c_str(),std::ios::in|std::ios::binary);
    label = parseTitle(ifs);
    names = parseVarNames(ifs);
    const int nCompSurf = names.size();

    int nodesPerElt;
    ifs >> nElts;
    ifs >> nodesPerElt;

    FArrayBox tnodes;
    tnodes.readFrom(ifs);
    const int nNodes = tnodes.box().numPts();

    // "rotate" the data so that the components are 'in the right spot for fab data'
    nodes.resize(tnodes.box(),nCompSurf);
    Real** np = new Real*[nCompSurf];
    for (int j=0; j<nCompSurf; ++j)
        np[j] = nodes.dataPtr(j);

    Real* ndat = tnodes.dataPtr();
    for (int i=0; i<nNodes; ++i)
    {
        for (int j=0; j<nCompSurf; ++j)
        {
            np[j][i] = ndat[j];
        }
        ndat += nCompSurf;
    }
    delete [] np;
    tnodes.clear();

    faceData.resize(nElts*nodesPerElt,0);
    ifs.read((char*)faceData.dataPtr(),sizeof(int)*faceData.size());
}

SegList::iterator FindMySeg(SegList& segs, int idx)
{
    const PMapIt& vertToFind = vertVec[idx];
    for (SegList::iterator it=segs.begin(); it!=segs.end(); ++it)
    {
        if ( ((*it)[0] == vertToFind) || ((*it)[1] == vertToFind) )
            return it;
    }
    return segs.end();
}


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
    const Point& p0 = (*p[0]).second.first;
    const Point& p1 = (*p[1]).second.first;
    mLength = std::sqrt(((p1[xComp] - p0[xComp])*(p1[xComp] - p0[xComp]))
                        +((p1[yComp] - p0[yComp])*(p1[yComp] - p0[yComp])));
}    

void
Segment::flip()
{
    PMapIt ptmp = p[0];
    p[0] = p[1];
    p[1] = ptmp;
}    

/*
   Linearly interpolate the position where an isosurface cuts
   an edge between two vertices, each with their own scalar value
*/
Point VI_doIt(Real isoVal,int isoComp,const vector<Point>& pts,int p1,int p2)
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
    
    Point res(pt1.size());
    const Real mu = (isoVal - valp1) / (valp2 - valp1);
    for (int j=0; j<res.size(); ++j)
        res[j] = pt1[j] + mu * (pt2[j] - pt1[j]);
    
    return res;
}

PMapIt VertexInterp(Real isoVal,int isoComp,const vector<Point>& pts,int p1,int p2)
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
            Nreuse++;
            return rev;
        }
    }
    Nreuse++;
    return fwd;
}

/*
   Given a cell and an isoVal, calculate the line segments
   required to represent the contour through the cell.
   Return an array of line segments
*/
Array<Segment> Segmentise(const vector<Point>&  pts,
                          const vector<int>&    elt,
                          Real                  isoVal,
                          int                   isoComp)
{
   BL_ASSERT(elt.size()==3);
   Array<PMapIt> vertlist(2);
   Array<Segment> segments;

   const int p0 = elt[0];
   const int p1 = elt[1];
   const int p2 = elt[2];

   int segCase = 0;
   if (pts[p0][isoComp] < isoVal) segCase |= 1;
   if (pts[p1][isoComp] < isoVal) segCase |= 2;
   if (pts[p2][isoComp] < isoVal) segCase |= 4;

   if (segCase==0 || segCase==7)
       return segments;

   switch (segCase) {
   case 1:
   case 6:
       vertlist[0] = VertexInterp(isoVal,isoComp,pts,p0,p1);
       vertlist[1] = VertexInterp(isoVal,isoComp,pts,p0,p2);
       break;
   case 2:
   case 5:
       vertlist[0] = VertexInterp(isoVal,isoComp,pts,p1,p0);
       vertlist[1] = VertexInterp(isoVal,isoComp,pts,p1,p2);
       break;
   case 3:
   case 4:
       vertlist[0] = VertexInterp(isoVal,isoComp,pts,p2,p0);
       vertlist[1] = VertexInterp(isoVal,isoComp,pts,p2,p1);
       break;
   }

   segments.resize(1);
   segments[0][0] = vertlist[0];
   segments[0][1] = vertlist[1];

   return segments;
}
