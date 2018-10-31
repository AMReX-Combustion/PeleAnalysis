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

    string infileL; pp.get("infileL",infileL);
    string infileR; pp.get("infileR",infileR);
    string outfile; pp.get("outfile",outfile);

    int nEltsL,nEltsR;
    FArrayBox nodesL,nodesR;
    Array<int> faceDataL,faceDataR;
    vector<string> namesL,namesR;
    string labelL,labelR;
    read_iso(infileL,nodesL,faceDataL,nEltsL,namesL,labelL);
    read_iso(infileR,nodesR,faceDataR,nEltsR,namesR,labelR);
    int nodesPerElt = faceDataL.size() / nEltsL;
    BL_ASSERT(nodesPerElt*nElts == faceData.size());

    nEltsL = faceDataL.size() / nodesPerElt;
    nEltsR = faceDataR.size() / nodesPerElt;

    if (nEltsL != nEltsR)
    {
        BoxLib::Abort("MEF files not compible");
    }

    int nCompMEFL = nodesL.nComp();
    BL_ASSERT(nEltsL*nodesPerElt == faceDataL.size());
    int nCompMEFR = nodesR.nComp();
    BL_ASSERT(nEltsR*nodesPerElt == faceDataR.size());

    Array<int> compsL;
    int nCompL = pp.countval("compsL");
    if (nCompL>0)
    {
        compsL.resize(nCompL);
        pp.getarr("compsL",compsL,0,nCompL);
    }
    else
    {
        int sCompL = 0;
        pp.query("sCompL",sCompL);
        nCompL = namesL.size();
        pp.query("nCompL",nCompL);
        BL_ASSERT(sCompL+nCompL <= nCompMEFL);
        compsL.resize(nCompL);
        for (int i=0; i<nCompL; ++i)
            compsL[i] = sCompL + i;
    }

    Array<int> compsR;
    int nCompR=pp.countval("compsR");
    if (nCompR > 0)
    {
        compsR.resize(nCompR);
        pp.getarr("compsR",compsR,0,nCompR);
    }
    else
    {
        int sCompR = 0;
        pp.query("sCompR",sCompR);
        nCompR = namesR.size();
        pp.query("nCompR",nCompR);
        BL_ASSERT(sCompR+nCompR <= nCompMEFR);
        compsR.resize(nCompR);
        for (int i=0; i<nCompR; ++i)
            compsR[i] = sCompR + i;
    }

    int nCompOut = compsL.size()+compsR.size();
    FArrayBox nodesOut(nodesL.box(),nCompOut);
    vector<string> namesOut(nCompOut);
    for (int i=0; i<compsL.size(); ++i)
    {
        std::cout << "(L): " << namesL[compsL[i]] << std::endl;
        namesOut[i] = namesL[compsL[i]];
        nodesOut.copy(nodesL,compsL[i],i,1);
    }
    for (int i=0; i<compsR.size(); ++i)
    {
        std::cout << "(R): " << namesR[compsR[i]] << std::endl;
        namesOut[compsL.size()+i] = namesR[compsR[i]];
        nodesOut.copy(nodesR,compsR[i],compsL.size()+i,1);
    }

    write_iso(outfile,nodesOut,faceDataL,nEltsL,namesOut,labelL);

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

