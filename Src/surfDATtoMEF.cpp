#include <string>
#include <iostream>

#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Utility.H>
using std::string;
using std::cout;
using std::endl;
using std::cerr;
using std::map;
using std::vector;
using namespace amrex;

#define SIZET int

static 
bool isNumberLine(const std::string& line)
{
    int i=0;
    while (i<line.size() && isspace(line[i]))
        i++;
    if (i==line.size()-1)
        return false;
    else if (line[i]=='-' || line[i]=='.')
        i++;

    if (i==line.size()-1)
        return false;

    return isdigit(line[i]) != 0;
}

Vector<string>
GetVarNames(std::istream& is, string& currentLine)
{
    std::getline(is,currentLine);
    vector<string> tokens;
    string buf;
    bool foundZone = false;
    while (!foundZone)
    {
        std::getline(is,buf);
        tokens = Tokenize(buf,"= ,");
        if (tokens[0]=="ZONE")
            foundZone = true;
        else
            currentLine += string(" ") + buf;
    }
    tokens = Tokenize(currentLine," =");

    Vector<std::string> names;
    for (int i=0; i<tokens.size(); ++i)
    {
        if (tokens[i] == "VARIABLES")
        {
            int nComp = tokens.size()-i-1;
            names.resize(nComp);
            for (int j=0; j<nComp; ++j)
            {
                names[j] = tokens[i+1+j];
            }
        }
    }
    currentLine = buf;
    return names;
}

map<std::string,std::string>
GetZoneParams(std::istream& is, int nComp, string& currentLine)
{
    map<std::string,std::string> zoneParams; // The result

    string txtLine = currentLine;
    bool foundData = false;
    while (!foundData)
    {
        std::getline(is,currentLine);

        if (!is)
            return zoneParams;

        if (isNumberLine(currentLine))
            foundData = true;
        else
            txtLine = txtLine + string(" ") + currentLine;
    }
    
    vector<string> tokens = Tokenize(txtLine,"\"");

    if (tokens.size()>1)
    {
        for (int i=0; i<tokens.size(); i=i+2)
        {
            if (i+2<=tokens.size())
            {
                vector<std::string> subTokens=Tokenize(tokens[i],", ");

                vector<std::string> tmp = Tokenize(subTokens[subTokens.size()-1],"=");
                
                string key = tmp[tmp.size()-1]; // key associated with this string parameter
                
                zoneParams[key] = "\"" + tokens[i+1] + "\"";
            }
        
            vector<std::string> subTokens=Tokenize(tokens[i],", =");
            int jst = (i==0 ? 1 : 0);
            for (int j=jst; j<subTokens.size()-1; j=j+2)
            {
                string key = subTokens[j];

                if (key=="DT")
                {
                    string tmpBuf;
                    for (int jj=1; jj<nComp+1; ++jj)
                        tmpBuf += subTokens[j+jj] + " ";
                    int nSkip = nComp;
                    if (subTokens[j+nComp+1] == ")")
                    {
                        nSkip++;
                        tmpBuf += ")";
                    }
                    zoneParams[key] = tmpBuf;
                    j = j+nComp;
                }
                else
                    zoneParams[key] = subTokens[j+1];
            }
        }
    }
    else
    {
        tokens = Tokenize(txtLine,", =");
        for (int i=1; i<tokens.size()-1; i=i+2)  // Skip the ZONE word
        {
            string key = tokens[i];
            if (key=="DT")
            {
                std::string dataTypeTokens=Tokenize(tokens[i+1],"()")[1];
                i = i+nComp;
            }
            else
                zoneParams[tokens[i]] = tokens[i+1];
        }
    }
    return zoneParams;
}

static
Real triangleArea(const Vector<Real>& p0,
                  const Vector<Real>& p1,
                  const Vector<Real>& p2)
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

class FABdata
{
public:
    FABdata(SIZET i, SIZET n)
        {fab.resize(Box(IntVect::TheZeroVector(),
                        IntVect(AMREX_D_DECL(i-1,0,0))),n);}
    Real& operator[](SIZET i) {return fab.dataPtr()[i];}
    FArrayBox fab;
    SIZET boxSize;
};

void
writeBin(const Vector<string>&         names,
         const Vector<Vector<Real> >&  nodeVec,
         const Vector<Vector<SIZET> >& eltVec,
         const std::string&            outfile,
         const std::string&            label);

int
main (int   argc,
      char* argv[])
{
    Initialize(argc,argv);

    ParmParse pp;

    std::string infile; pp.get("infile",infile);

    // Build outfile name (may modify below if more than one dataset)
    vector<string> infileTokens = Tokenize(infile,".");
    string outfile;
    outfile = infileTokens[0];
    for (int i=1; i<infileTokens.size()-1; ++i)
        outfile += string(".") + infileTokens[i];
    outfile += ".mef";
    pp.query("outfile",outfile);

    std::ifstream ifs;
    ifs.open(infile.c_str());
    string buf;

    //Read first zone info
    Vector<string> names = GetVarNames(ifs,buf);
    int nComp = names.size();

    map<std::string,std::string> zoneParams = GetZoneParams(ifs,nComp,buf);

    Real areaEps = 1.e-12; pp.query("areaEps",areaEps);

    int zoneID = 0;
    while (zoneParams.size() > 0)
    {
        // Find N and E in zone parameters
        SIZET N,E;
        string eltType("undefined");
        string label("LabelHere");
        
        for (std::map<std::string,std::string>::const_iterator it=zoneParams.begin(); it!=zoneParams.end(); ++it)
        {
            if (it->first == "N")
                N = atoi(it->second.c_str());
            
            if (it->first == "E")
                E = atoi(it->second.c_str());
            
            if (it->first == "ET")
                eltType = it->second.c_str();

            if (it->first == "T")
                label = it->second.c_str();
        }

        if (ParallelDescriptor::IOProcessor())
            cerr << "From " << infile << ": " << N << " nodes and " << E << " elements with "
                 << nComp << " unknowns per node.  Dataset labelled: " << label << endl;
        
        Vector<Vector<Real> > nodeVec(N);
        vector<string> tokens;
        for (SIZET j=0; j<N; ++j)
        {
            Vector<Real>& v = nodeVec[j];
            v.resize(nComp);
            tokens = Tokenize(buf,", ");
            BL_ASSERT(tokens.size()==nComp);
            for (int k=0; k<nComp; ++k)
            {
                v[k] = atof(tokens[k].c_str());
            }
            std::getline(ifs,buf);
        }
        
        int nodesPerElt = -1;
        if (eltType == "TRIANGLE")
            nodesPerElt = 3;
        
        Vector<Vector<SIZET> > eltVec(E);
        for (SIZET j=0; j<E; ++j)
        {
            Vector<SIZET>& e = eltVec[j];
            tokens = Tokenize(buf,", ");
            BL_ASSERT(tokens.size()==nodesPerElt);
            e.resize(nodesPerElt);
            for (int k=0; k<nodesPerElt; ++k)
                e[k] = atoi(tokens[k].c_str());
            std::getline(ifs,buf);
        }

        if (ParallelDescriptor::IOProcessor())
            cerr << "...finished reading data" << endl;

        Real area = 0;
        map<Vector<int>,Real > bins;
        for (int i=0; i<E; ++i)
        {
            const Vector<SIZET>& elt = eltVec[i];
            const Vector<Real>& A = nodeVec[elt[0]-1];
            const Vector<Real>& B = nodeVec[elt[1]-1];
            const Vector<Real>& C = nodeVec[elt[2]-1];
                        
            area += triangleArea(A,B,C);
        }

        std::cout << "zoneID, area = " << zoneID << ", " << area << std::endl;

        string binout = outfile;
        if (zoneID>0)
        {
            char buf1[72];
            sprintf(buf1, "%d", zoneID); 
            binout = infileTokens[0];
            for (int i=1; i<infileTokens.size()-1; ++i)
                outfile += string(".") + infileTokens[i];
            binout += string("_") + string(buf1) + ".mef";
        }

        writeBin(names,nodeVec,eltVec,binout,label);

        zoneID++;
        zoneParams = GetZoneParams(ifs,nComp,buf);

    }

    Finalize();
    return 0;
}

void
writeBin(const Vector<string>&         names,
         const Vector<Vector<Real> >&  nodeVec,
         const Vector<Vector<SIZET> >& eltVec,
         const std::string&            outfile,
         const std::string&            label)
{
    std::string vars;
    for (int j=0; j<names.size(); ++j)
        vars += " " + names[j];

    FABdata tmpData(nodeVec.size(),nodeVec[0].size());
    SIZET cnt = 0;
    for (SIZET i=0; i<nodeVec.size(); ++i)
    {
        const Vector<Real>& vec = nodeVec[i];
        for (int n=0; n<vec.size(); ++n)
            tmpData[cnt++] = vec[n];
    }


    SIZET nElts = eltVec.size();
    cnt = 0;
    int nodesPerElt = eltVec[0].size();
    Vector<SIZET> connData(nodesPerElt*nElts);
    for (SIZET i=0; i<nElts; ++i)
    {
        const Vector<SIZET>& elt = eltVec[i];
        BL_ASSERT(elt.size()==nodesPerElt);
        for (int j=0; j<elt.size(); ++j)
            connData[cnt++] = elt[j];
    }

    std::ofstream ofs;
    std::ostream* os =
        (outfile=="-" ? (std::ostream*)(&std::cout) : (std::ostream*)(&ofs) );
    if (outfile!="-")
        ofs.open(outfile.c_str(),std::ios::out|std::ios::trunc|std::ios::binary);
    (*os) << label << endl;
    (*os) << vars << endl;
    (*os) << nElts << " " << nodesPerElt << endl;
    tmpData.fab.writeOn(*os);
    (*os).write((char*)connData.dataPtr(),sizeof(SIZET)*connData.size());
    if (outfile!="-")
        ofs.close();
}

