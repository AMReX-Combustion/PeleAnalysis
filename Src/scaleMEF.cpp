#include <string>
#include <iostream>
#include <set>
#include <map>
#include <vector>

#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Utility.H>
#include <AMReX_VisMF.H>

using namespace amrex;

using std::vector;
using std::map;
using std::string;
using std::cerr;
using std::endl;
using std::cout;
using std::ofstream;

void
read_iso(const std::string& infile,
         FArrayBox&         nodes,
         Vector<int>&       faceData,
         int&               nElts,
         vector<string>&    names,
         string&            label);

void
write_iso(const std::string&    outfile,
          const FArrayBox&      nodes,
          const Vector<int>&    faceData,
          int                   nElts,
          const vector<string>& names,
          const string&         label);

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

int
main (int   argc,
      char* argv[])
{
    amrex::Initialize(argc,argv);

    ParmParse pp;

    int nElts;
    string infile; pp.get("infile",infile);
    string outfile; pp.get("outfile",outfile);

    FArrayBox nodes;
    Vector<int> faceData;
    vector<string> names;
    string label;
    read_iso(infile,nodes,faceData,nElts,names,label);
    int nodesPerElt = faceData.size() / nElts;
    BL_ASSERT(nodesPerElt*nElts == faceData.size());

    nElts = faceData.size() / nodesPerElt;
    int nCompMEF = nodes.nComp();
    BL_ASSERT(nElts*nodesPerElt == faceData.size());

    Vector<int> comps;
    if (pp.countval("comps"))
    {
        int nComp = pp.countval("comps");
        comps.resize(nComp);
        pp.getarr("comps",comps,0,nComp);
    }
    else
    {
        int sComp = 0;
        pp.query("sComp",sComp);
        int nc = 1;
        pp.query("nComp",nc);
        BL_ASSERT(sComp+nc <= nCompMEF);
        comps.resize(nc);
        for (int i=0; i<nc; ++i)
            comps[i] = sComp + i;
    }
    int nComp = comps.size();

    Vector<Real> vals(nComp);
    BL_ASSERT(pp.countval("vals")==nComp);
    pp.getarr("vals",vals,0,nComp);

    int nNodes = nodes.box().numPts();
    for (int j=0; j<nComp; ++j)
    {
        BL_ASSERT(comps[j]<nCompMEF);
        Real* dat=nodes.dataPtr(comps[j]);
        for (int i=0; i<nNodes; ++i)
            dat[i] = dat[i]*vals[j];
    }

    if (pp.countval("newNames"))
    {
        int numNew = pp.countval("newNames");
        int numNewComps = pp.countval("newComps");
        BL_ASSERT(numNew==numNewComps);

        Vector<std::string> newNames(numNew);
        Vector<int> newComps(numNewComps);
        pp.getarr("newNames",newNames,0,newNames.size());
        pp.getarr("newComps",newComps,0,newComps.size());

        for (int i=0; i<numNew; ++i) {
            BL_ASSERT(newComps[i] <= names.size());
            names[newComps[i]] = newNames[i];
        }       
    }

    write_iso(outfile,nodes,faceData,nElts,names,label);

    amrex::Finalize();
    return 0;
}

void
write_iso(const std::string&    outfile,
          const FArrayBox&      nodes,
          const Vector<int>&    faceData,
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
         Vector<int>&       faceData,
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

