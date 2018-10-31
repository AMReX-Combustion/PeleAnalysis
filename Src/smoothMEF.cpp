#include "winstd.H"

#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <set>
#include "ParmParse.H"
#include "MultiFab.H"
#include "VisMF.H"
#include "ParallelDescriptor.H"
#include "DataServices.H"
#include "Geometry.H"
#include "Utility.H"
using std::vector;
using std::set;
using std::string;
using std::cerr;
using std::endl;
using std::cout;
using std::ofstream;

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

static 
Array<Array<int> > 
buildNodeNeighbors(const Array<int>& faceData,
                   int               nNodes)
{
    Array<Array<int> > nodeNeighbors(nNodes);
    int nElts = faceData.size()/3;
    for (int i=0; i<nElts; ++i)
    {
        int offset=i*3;
        for (int j=0; j<3; ++j)
        {
            Array<int>& neighbors = nodeNeighbors[ faceData[offset+j]-1 ];
            neighbors.resize(neighbors.size()+1);
            neighbors[neighbors.size()-1] = i; // b/c of how this is done, we know the list will be unique
        }
    }

    Array<Array<int> > cellNeighbors(nElts);
    for (int i=0; i<nElts; ++i)
    {
        int offset=i*3;
        set<int> n;
        for (int j=0; j<3; ++j)
        {
            const Array<int>& neighbors = nodeNeighbors[ faceData[offset+j]-1 ];
            for (int k=0; k<neighbors.size(); ++k)
            {
                int nc = neighbors[k];
                if (nc != i)
                    n.insert(nc);
            }
        }
        cellNeighbors[i].resize(n.size());
        int cnt=0;
        for (std::set<int>::const_iterator it=n.begin(); it!=n.end(); ++it)
            cellNeighbors[i][cnt++] = *it;
    }
    return cellNeighbors;
}

void
smoothVals(Array<Real>&              newVals,
           const Array<Real>&        vals,
           const Array<Real>&        area,
           const Array<Array<int> >& nodeNeighbors)
{
    int nElts = vals.size();
    if (nElts!=newVals.size() || nElts!=area.size())
        BoxLib::Abort("Data wrong size");

    for (int i=0; i<nElts; ++i)
    {
        const Array<int> neighbors = nodeNeighbors[i];

        Real accumArea = area[i];
        for (int j=0; j<neighbors.size(); ++j)
            accumArea += area[ neighbors[j] ];

        Real accumWt = vals[i];
        for (int j=0; j<neighbors.size(); ++j)
            accumWt += vals[ neighbors[j] ];

        newVals[i] = accumWt / accumArea;
    }
}

void
triangle_area(const Real*       x,
              const Real*       y,
              const Real*       z,
              const Array<int>& faceData,
              int               nElts,
              int               nNodes,
              Real*             area)
{
    Real R1[3], R2[3], R3[3];
    const Real* loc[3];
    loc[0] = x;
    loc[1] = y;
    loc[2] = z;

    for (int i=0; i<nElts; ++i)
    {
        int offset=i*3;
        int pa = faceData[offset+0] - 1;
        int pb = faceData[offset+1] - 1;
        int pc = faceData[offset+2] - 1;
        for (int j=0; j<3; ++j)
        {
            R1[j] = loc[j][pb] - loc[j][pa];
            R2[j] = loc[j][pc] - loc[j][pa];
        }
        R3[0] = R1[1]*R2[2] - R2[1]*R1[2];
        R3[1] = R1[2]*R2[0] - R2[2]*R1[0];
        R3[2] = R1[0]*R2[1] - R2[0]*R1[1];

        Real res=0;
        for (int j=0; j<3; ++j)
            res += R3[j]*R3[j];
        area[i] = 0.5*std::sqrt(res);
    }
}

int
main (int   argc,
      char* argv[])
{
    BoxLib::Initialize(argc,argv);

    ParmParse pp;

    int nElts;
    string infile; pp.get("infile",infile);
    string outfile; pp.get("outfile",outfile);

    FArrayBox nodes;
    Array<int> faceData;
    vector<string> names;
    string label;
    read_iso(infile,nodes,faceData,nElts,names,label);
    int nodesPerElt = faceData.size() / nElts;
    BL_ASSERT(nodesPerElt*nElts == faceData.size());

    nElts = faceData.size() / nodesPerElt;
    int nNodes = nodes.box().numPts();
    int nCompMEF = nodes.nComp();
    BL_ASSERT(nElts*nodesPerElt == faceData.size());

    int areaComp=-1; pp.query("areaComp",areaComp); // Treated as an area, rather than computing here on the fly
    int comp=-1; pp.get("comp",comp); // Component to smooth

    Array<Real> vals(nElts); // Begin with the assumption that we associate values on this surface with the elements (vs. nodes)
    Array<Real> area(nElts);

    // Initialize data by averaging node-based data to elements
    Real* dataN = nodes.dataPtr(comp);
    Array<Real> tmp(nElts);
    Real* areaN;
    if (areaComp>=0 && areaComp<nCompMEF)
    {
        areaN = nodes.dataPtr(areaComp);
    }
    else
    {
        int idX=0;
        int idY=1;
        int idZ=2;

        areaN = tmp.dataPtr();
        triangle_area(nodes.dataPtr(idX),nodes.dataPtr(idY),nodes.dataPtr(idZ),faceData,nElts,nNodes,areaN);
    }

    for (int i=0; i<nElts; ++i)
    {
        int offset=i*nodesPerElt;
        vals[i] = 0;
        area[i] = 0;
        for (int j=0; j<nodesPerElt; ++j)
        {
            area[i] += areaN[ faceData[offset+j]-1 ];
            vals[i]  += dataN[ faceData[offset+j]-1 ] * area[i];
        }
        vals[i] *= 1./nodesPerElt;
        area[i] *= 1./nodesPerElt;
    }

    Array<Array<int> > nodeNeighbors = buildNodeNeighbors(faceData,nodes.box().numPts());

    Array<Real> newVals(nElts);

    int nSmooth=1; pp.query("nSmooth",nSmooth);

    for (int j=0; j<nSmooth; ++j)
    {
        smoothVals(newVals,vals,area,nodeNeighbors);
        for (int i=0; i<nElts; ++i)
            vals[i] = newVals[i];
    }

    for (int i=0; i<nElts; ++i)
        dataN[i] = vals[i]/area[i];

    write_iso(outfile,nodes,faceData,nElts,names,label);

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

