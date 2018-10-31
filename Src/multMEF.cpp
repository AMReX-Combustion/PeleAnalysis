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
    int nCompMEF = nodes.nComp();
    BL_ASSERT(nElts*nodesPerElt == faceData.size());

    Array<int> comps;
    int nComp=0;
    if (nComp = pp.countval("comps"))
    {
        comps.resize(nComp);
        pp.getarr("comps",comps,0,nComp);
    }
    else
    {
        int sComp = 0;
        pp.query("sComp",sComp);
        nComp = 1;
        pp.query("nComp",nComp);
        BL_ASSERT(sComp+nComp <= nCompMEF);
        comps.resize(nComp);
        for (int i=0; i<nComp; ++i)
            comps[i] = sComp + i;
    }

    FArrayBox nodesOut(nodes.box(),1);
    vector<string> namesOut(1);
    namesOut[0] = "product"; pp.query("nameOut",namesOut[0]);

    nodesOut.setVal(1.0); // Initialize to 1
    Real* datOut = nodesOut.dataPtr(0);
    int nNodes = nodes.box().numPts();
    for (int j=0; j<nComp; ++j)
    {
        BL_ASSERT(comps[j]<nCompMEF);
        Real* dat=nodes.dataPtr(comps[j]);
        for (int i=0; i<nNodes; ++i)
            datOut[i] *= dat[i];
    }

    write_iso(outfile,nodesOut,faceData,nElts,namesOut,label);

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

