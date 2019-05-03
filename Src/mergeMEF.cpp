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
void
merge_iso(FArrayBox&        nodesM,
          Array<int>&       faceDataM,
          int&              nEltsM,
          const FArrayBox&  nodes,
          const Array<int>& faceData,
          int               nElts,
          Real              eps,
          int               verbose,
          bool              remDupNodes);

int
main (int   argc,
      char* argv[])
{
    BoxLib::Initialize(argc,argv);
    {
    ParmParse pp;

    int verbose=0; pp.query("verbose",verbose);
    Real eps=1.e-8; pp.query("eps",eps);
    bool remDupNodes = false; pp.query("remDupNodes",remDupNodes);
    
    Array<string> infiles;
    int nf=pp.countval("infiles");
    infiles.resize(nf);
    pp.getarr("infiles",infiles,0,nf);
    string outfile; pp.get("outfile",outfile);

    FArrayBox nodesM,nodes;
    Array<int> faceDataM, faceData;
    vector<string> namesM,names;
    string labelM, label;
    int nEltsM, nElts;
    for (int i=0; i<infiles.size(); ++i)
    {
        if (verbose)
            cout << "Reading " << infiles[i] << endl;

        read_iso(infiles[i],nodes,faceData,nElts,names,label);
        if (i==0)
        {
            int nComp = nodes.nComp();
            nodesM.resize(nodes.box(),nComp);
            nodesM.copy(nodes,0,0,nComp);
            namesM.resize(nComp);
            for (int j=0; j<nComp; ++j)
            {
                namesM[j] = names[j];
            }
            nEltsM = nElts;
            labelM = label;
            faceDataM.resize(faceData.size());
            for (int j=0; j<faceData.size(); ++j)
            {
                faceDataM[j] = faceData[j];
            }
        }
        else
        {
            // For now, simply assert that the files have the same components
            bool sameComps = true;
            if (namesM.size() == names.size())
            {
                for (int j=0; j<names.size() && sameComps; ++j)
                    sameComps &= names[j]==namesM[j];
            }
            if (!sameComps)
                BoxLib::Abort("Incompatible files");

            if (verbose)
                cout << "  merging surfaces" << endl;

            merge_iso(nodesM,faceDataM,nEltsM,nodes,faceData,nElts,eps,verbose,remDupNodes);
        }
    }

    cout << "Writing " << outfile << endl;

    write_iso(outfile,nodesM,faceDataM,nEltsM,namesM,labelM);
    }
    BoxLib::Finalize();
    return 0;
}

static
void
merge_iso(FArrayBox&        nodesM,
          Array<int>&       faceDataM,
          int&              nEltsM,
          const FArrayBox&  nodes,
          const Array<int>& faceData,
          int               nElts,
          Real              eps,
          int               verbose,
          bool              remDupNodes)
{
    int nodesPerEltM = (int)(faceDataM.size() / nEltsM);
    if (nodesPerEltM * nEltsM != faceDataM.size())
        BoxLib::Abort("Non-constant nodes/elt not supported");

    int nodesPerElt = (int)(faceData.size() / nElts);
    if (nodesPerElt * nElts != faceData.size())
        BoxLib::Abort("Non-constant nodes/elt not supported");

    if (nodesPerEltM != nodesPerElt)
        BoxLib::Abort("Incompatible surfaces");

    const Box& box = nodes.box();
    int nNodes;
    if (box.numPtsOK() && box.numPts()==(int)(box.numPts()))
        nNodes = box.numPts();
    else
        BoxLib::Abort("New surface too big");

    const Box& boxM = nodesM.box();
    long nNodesM;
    if (boxM.numPtsOK() && boxM.numPts()==(int)(boxM.numPts()))
        nNodesM = boxM.numPts();
    else
        BoxLib::Abort("Old surface too big");

    int nComp = nodes.nComp();
    vector<Real*> dptrM(nComp);
    vector<const Real*> dptr(nComp);
    for (int j=0; j<nComp; ++j)
    {
        dptrM[j] = nodesM.dataPtr(j);
        dptr[j] = nodes.dataPtr(j);
    }
    Array<int> newNodes(nNodes);
    Real epsSq = eps*eps;
    int cnt = nNodesM;
    for (int i=0; i<nNodes; ++i)
    {
        bool found = false;
        if (remDupNodes)
        {
            for (int j=0; j<nNodesM && !found; ++j)
            {
                Real dist = 0;
                for (int k=0; k<BL_SPACEDIM; ++k)
                {
                    Real dx = dptr[k][i] - dptrM[k][j];
                    dist += dx*dx;
                }
                if (dist <= epsSq)
                {
                    found = true;
                    newNodes[i] = j;
                }
            }
        }
        if (!found)
            newNodes[i] = cnt++;
    }

    if (cnt>nNodesM)
    {
        FArrayBox old(nodesM.box(),nComp);
        old.copy(nodesM);
        nodesM.resize(Box(boxM).growHi(0,cnt-nNodesM),nComp);
        nodesM.copy(old);
        old.clear();
        for (int k=0; k<nComp; ++k)
        {
            Real* d = nodesM.dataPtr(k);
            for (int i=0; i<newNodes.size(); ++i)
            {
                if (newNodes[i]>=nNodesM)
                {
                    d[newNodes[i]] = dptr[k][i];
                }
            }
        }

        int offset = faceDataM.size();
        Array<int> oldFD(offset);
        for (int i=0; i<offset; ++i)
            oldFD[i] = faceDataM[i];

        faceDataM.resize(offset+faceData.size());
        for (int i=0; i<offset; ++i)
            faceDataM[i] = oldFD[i];
        for (int i=0; i<faceData.size(); ++i)
        {
            faceDataM[offset+i] = newNodes[faceData[i] - 1] + 1;
        }
        oldFD.clear();
        nEltsM = (int) (faceDataM.size() / nodesPerElt);
    }
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

