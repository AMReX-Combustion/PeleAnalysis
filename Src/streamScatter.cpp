#include "AMReX_ParmParse.H"
#include "AMReX_MultiFab.H"
#include "AMReX_VisMF.H"
#include "AMReX_ParallelDescriptor.H"
#include "AMReX_Utility.H"
using std::set;
using std::list;
using std::vector;
using std::map;
using std::string;
using std::cerr;
using std::endl;
using std::cout;
using std::ofstream;

using namespace amrex;

struct MLloc
{
    MLloc() :
        amr_lev(-1), box_idx(-1), pt_idx(-1) {}
    MLloc(int lev,int box,int pt) :
        amr_lev(lev), box_idx(box), pt_idx(pt) {}
    int amr_lev, box_idx, pt_idx;
};

int
ml_Nlevels(const std::string& infile);

void
read_ml_streamline_data(const std::string& infile,
                        Vector<MultiFab*>& data,
                        vector<string>&    names,
                        Vector<int>&       faceData,
                        int&               nElts,
                        Vector<Vector<Vector<int> > >& inside_nodes);
std::vector<MLloc>
build_nodeMap(const Vector<Vector<Vector<int> > >& inside_nodes);

int get_nPts(const Vector<MultiFab*>& data);
int get_jlo(const Vector<MultiFab*>& data);

Real
wedge_volume_int(const MLloc& p1,
                 const MLloc& p2,
                 const MLloc& p3,
                 int          ptOnStr,
                 int          comp,
                 const Vector<MultiFab*>& data,
                 int idX, int idY, int idZ);

int
main (int   argc,
      char* argv[])
{
    Initialize(argc,argv);

    ParmParse pp;

    std::string infile; pp.get("infile",infile);

    // Build outfile name (may modify below if more than one dataset)
    vector<string> infileTokens = Tokenize(infile,".");
    string outfileBase = infileTokens[0];
    for (int i=1; i<infileTokens.size()-1; ++i)
        outfileBase += string(".") + infileTokens[i];
    outfileBase += "/JPDF_"; 
    pp.query("outfileBase",outfileBase);

    int verbose=0; pp.query("verbose",verbose);

    int Nlev = ml_Nlevels(infile);
    Vector<MultiFab*> streamlines(Nlev);
    Vector<std::string> names;
    Vector<int> faceData;
    int nElts;
    Vector<Vector<Vector<int> > > inside_nodes;
    if (verbose)
        std::cerr << "Reading stream file ...\n";
    read_ml_streamline_data(infile,streamlines,names,faceData,nElts,inside_nodes);
    if (verbose)
        std::cerr << "...finished reading stream file \n";

    Vector<int> comps;
    int nc = pp.countval("vars");
    Vector<string> vars(nc);
    comps.resize(nc,-1);
    pp.getarr("vars",vars,0,nc);
    for (int i=0; i<nc; ++i) {
        for (int j=0; comps[i]<0 && j<names.size(); ++j) {
            if (names[j]==vars[i]) {
                comps[i] = j;
            }
        }
        if (comps[i] < 0) {
            Abort(string(string("Variable not found: ")+vars[i]).c_str());
        }
    }

    // Build a structure that for each node points to where in the streamlines to get the data
    std::vector<MLloc> nodeMap = build_nodeMap(inside_nodes);

    int condComp = -1; pp.query("condComp",condComp);
    string condVar = ""; pp.query("condVar",condVar);
    if (condVar != "" ) {
        for (int j=0; j<names.size(); ++j) {
            if (names[j]==condVar) {
                condComp=j;
            }
        }
        if (condComp < 0) {
            Abort(string(string("Conditioning variable not found: ")+condVar).c_str());
        }
    }
        
    Real condValMoreThan = 0; pp.query("condValMoreThan",condValMoreThan);
    Real condValLessThan = 0; pp.query("condValLessThan",condValLessThan);

    for (int i=0; i<nodeMap.size(); ++i)
    {
        const MLloc& p = nodeMap[i];
        const int L = p.amr_lev;
        const int B = p.box_idx;
        const int Pt = p.pt_idx;
        
        const FArrayBox& fab = (*streamlines[L])[B];
        const Box& strBox = fab.box();
        
        int loPtOnStr = strBox.smallEnd()[1];
        int hiPtOnStr = strBox.bigEnd()[1];
        
        int ptOnStr = 0.5*(loPtOnStr + hiPtOnStr);
        Real hiVal = fab(IntVect(Pt,ptOnStr,0),condComp);
        for (int j=loPtOnStr; j<=hiPtOnStr; ++j)
        {
            IntVect iv = IntVect(Pt,j,0);
            const Real val = fab(iv,condComp);
            if (val > hiVal) 
            {
                hiVal = val;
                ptOnStr = j;
            }
        }
        
        if ( (hiVal >= condValMoreThan) && (hiVal < condValLessThan) )
        {
            // write scatter data
            IntVect iv(Pt,ptOnStr,0);
            for (int j=0; j<nc; ++j)
            {
                std::cout << fab(iv,comps[j]) << " ";
            }
            std::cout << std::endl;
        }
    }
    
    Finalize();
    return 0;
}

int
get_jlo(const Vector<MultiFab*>& data)
{
  int jlo = (*data[0])[0].box().smallEnd(1);
    for (int i=0; i<data.size(); ++i)
        for (int j=0; j<data[i]->size(); ++j)
            jlo = std::min(jlo,(*data[i])[j].box().smallEnd(1));
    return jlo;
}

int
get_nPts(const Vector<MultiFab*>& data)
{
    int nPts = 0;
    for (int i=0; i<data.size(); ++i)
        for (int j=0; j<data[i]->size(); ++j)
            nPts = std::max(nPts,(*data[i])[j].box().length(1));
    return nPts;
}

Real
Vdot(const vector<Real>& A,
     const vector<Real>& B)
{
    Real res=0;
    for (int i=0; i<3; ++i)
        res += A[i]*B[i];
    return res;
}

vector<Real>
Vcross(const vector<Real>& A,
       const vector<Real>& B)
{
    vector<Real> res(3);
    res[0] = A[1]*B[2] - B[1]*A[2];
    res[1] = A[2]*B[0] - B[2]*A[0];
    res[2] = A[0]*B[1] - B[0]*A[1];
    return res;
}

vector<Real>
Vminus(const vector<Real>& A,
       const vector<Real>& B)
{
    vector<Real> res(3);
    for (int i=0; i<3; ++i)
        res[i] = A[i] - B[i];
    return res;
}

std::vector<MLloc>
build_nodeMap(const Vector<Vector<Vector<int> > >& inside_nodes)
{
    // Count number of nodes total
    int num_nodes = 0;
    int Nlev = inside_nodes.size();
    for (int i=0; i<Nlev; ++i)
        for (int j=0; j<inside_nodes[i].size(); ++j)
            num_nodes += inside_nodes[i][j].size();

    std::vector<MLloc> nodeMap(num_nodes);
    for (int i=0; i<Nlev; ++i)
    {
        for (int j=0; j<inside_nodes[i].size(); ++j)
        {
            for (int k=0; k<inside_nodes[i][j].size(); ++k)
            {
                BL_ASSERT(inside_nodes[i][j][k]<=num_nodes);
                BL_ASSERT(inside_nodes[i][j][k]>0);
                nodeMap[inside_nodes[i][j][k] - 1] = MLloc(i,j,k); // Remember that inside_nodes is 1-based
            }
        }
    }

    return nodeMap;    
}

int
ml_Nlevels(const std::string& infile)
{
    // Read header for path traces
    const string Header = infile + std::string("/Header");
    std::ifstream ifs;
    ifs.open(Header.c_str());
    std::string typeLabel; ifs >> typeLabel;
    int NlevPath; ifs >> NlevPath;
    ifs.close();
    return NlevPath;
}

void
read_ml_streamline_data(const std::string& infile,
                        Vector<MultiFab*>& data,
                        vector<string>&    names,
                        Vector<int>&       faceData,
                        int&               nElts,
                        Vector<Vector<Vector<int> > >& inside_nodes)
{
    // Read header for path traces
    const string Header = infile + string("/Header");
    std::ifstream ifh;
    ifh.open(Header.c_str());
    std::string typeLabel; ifh >> typeLabel;
    int NlevPath; ifh >> NlevPath;
    int nCompPath; ifh >> nCompPath;
    names.resize(nCompPath);
    for (int i=0; i<nCompPath; ++i)
        ifh >> names[i];
    ifh.close();

    // Read elements
    const string iElements = infile + string("/Elements");
    std::ifstream ife;
    ife.open(iElements.c_str());
    ife >> nElts;
    int nodesPerElt; ife >> nodesPerElt;
    faceData.resize(nElts*nodesPerElt);
    for (int i=0; i<nElts*nodesPerElt; ++i)
        ife >> faceData[i];

    for (int lev=0; lev<NlevPath; ++lev)
    {
        char buf[64];
        sprintf(buf, "/Level_%d", lev);
        std::string ThisStreamFile = infile + std::string(buf);
        ThisStreamFile += "/Str";
        data[lev] = new MultiFab();
        VisMF::Read(*data[lev],ThisStreamFile);
    }
    //
    // Force other processors to wait until parallel data has been read
    //
    ParallelDescriptor::Barrier();

    // Read element distribution (after we know Nlev from path traces)
    inside_nodes.resize(NlevPath);
    for (int lev=0; lev<NlevPath; ++lev)
    {
        inside_nodes[lev].resize(data[lev]->size());

        // Get number of nonempty boxes
        int num_non_zero;
        ife >> num_non_zero;

        for (int j=0; j<num_non_zero; ++j)
        {
            // Get index of non-empty box, and number of entries
            int box_id, num_ids;
            ife >> box_id >> num_ids;

            inside_nodes[lev][box_id].resize(num_ids);

            for (int k=0; k<num_ids; ++k)
                ife >> inside_nodes[lev][box_id][k];
        }
    }
    ife.close();
}

Real
tetVol(const Real* A,
       const Real* B,
       const Real* C,
       const Real* D)
{
    Real R1[3], R2[3], R3[3], R4[3];
    for (int i=0; i<3; ++i)
    {
        R1[i] = D[i] - A[i];
        R2[i] = B[i] - A[i];
        R3[i] = C[i] - A[i];
    }
    for (int i=0; i<3; ++i)
    {
        R4[0] = R2[1]*R3[2] - R3[1]*R2[2];
        R4[1] = R2[2]*R3[0] - R3[2]*R2[0];
        R4[2] = R2[0]*R3[1] - R3[0]*R2[1];
    }
    Real res=0;
    for (int i=0; i<3; ++i)
        res += R1[i]*R4[i];
    return std::abs(res);
}

Real
wedge_volume_int(const MLloc& p1,
                 const MLloc& p2,
                 const MLloc& p3,
                 int          ptOnStr,
                 int          comp,
                 const Vector<MultiFab*>& data,
                 int idX, int idY, int idZ)
{
    const int L1 = p1.amr_lev;
    const int L2 = p2.amr_lev;
    const int L3 = p3.amr_lev;
    
    const int B1 = p1.box_idx;
    const int B2 = p2.box_idx;
    const int B3 = p3.box_idx;
    
    const int P1 = p1.pt_idx;
    const int P2 = p2.pt_idx;
    const int P3 = p3.pt_idx;

    Real A[3], B[3], C[3], D[3], E[3], F[3];

    const IntVect iv1(P1,ptOnStr,0);
    const IntVect iv2(P2,ptOnStr,0);
    const IntVect iv3(P3,ptOnStr,0);
    const IntVect iv1p(P1,ptOnStr+1,0);
    const IntVect iv2p(P2,ptOnStr+1,0);
    const IntVect iv3p(P3,ptOnStr+1,0);

    A[0] = (*data[L1])[B1](iv1,idX);
    A[1] = (*data[L1])[B1](iv1,idY);
    A[2] = (*data[L1])[B1](iv1,idZ);

    B[0] = (*data[L2])[B2](iv2,idX);
    B[1] = (*data[L2])[B2](iv2,idY);
    B[2] = (*data[L2])[B2](iv2,idZ);

    C[0] = (*data[L3])[B3](iv3,idX);
    C[1] = (*data[L3])[B3](iv3,idY);
    C[2] = (*data[L3])[B3](iv3,idZ);

    D[0] = (*data[L1])[B1](iv1p,idX);
    D[1] = (*data[L1])[B1](iv1p,idY);
    D[2] = (*data[L1])[B1](iv1p,idZ);

    E[0] = (*data[L2])[B2](iv2p,idX);
    E[1] = (*data[L2])[B2](iv2p,idY);
    E[2] = (*data[L2])[B2](iv2p,idZ);

    F[0] = (*data[L3])[B3](iv3p,idX);
    F[1] = (*data[L3])[B3](iv3p,idY);
    F[2] = (*data[L3])[B3](iv3p,idZ);

    // The following are actually 6 times the tet volume
    const Real vol_EABC = tetVol(A,B,C,E);
    const Real vol_ADEF = tetVol(A,D,E,F);
    const Real vol_ACEF = tetVol(C,E,F,A);

    Real result;
    if (comp<0)
    {
        result = (vol_EABC + vol_ADEF + vol_ACEF)/6.;
    }
    else
    {
        // The following are actually 6 times the tet volume
        const Real vol_DABC = tetVol(A,B,C,D);
        const Real vol_FABC = tetVol(A,B,C,F);
        const Real vol_BDEF = tetVol(B,D,E,F);
        const Real vol_CDEF = tetVol(C,D,E,F);
        const Real vol_ACED = tetVol(C,E,D,A);
        const Real vol_BCDF = tetVol(B,C,D,F);
        const Real vol_BCDE = tetVol(B,C,D,E);
        const Real vol_ABDF = tetVol(B,D,F,A);
        const Real vol_ABEF = tetVol(B,E,F,A);

        const Real vA = (*data[L1])[B1](IntVect(P1,ptOnStr,  0),comp);
        const Real vB = (*data[L2])[B2](IntVect(P2,ptOnStr,  0),comp);
        const Real vC = (*data[L3])[B3](IntVect(P3,ptOnStr,  0),comp);                    
        const Real vD = (*data[L1])[B1](IntVect(P1,ptOnStr+1,0),comp);
        const Real vE = (*data[L2])[B2](IntVect(P2,ptOnStr+1,0),comp);
        const Real vF = (*data[L3])[B3](IntVect(P3,ptOnStr+1,0),comp);

        // These are actually 24 times the integral
        const Real int_1 = ( (vD+vA+vB+vC)*vol_DABC + 
                             (vB+vD+vE+vF)*vol_BDEF + 
                             (vB+vC+vD+vF)*vol_BCDF );

        const Real int_2 = ( (vD+vA+vB+vC)*vol_DABC + 
                             (vC+vD+vE+vF)*vol_CDEF + 
                             (vB+vC+vD+vE)*vol_BCDE );

        const Real int_3 = ( (vE+vA+vB+vC)*vol_EABC + 
                             (vA+vD+vE+vF)*vol_ADEF + 
                             (vA+vC+vE+vF)*vol_ACEF );

        const Real int_4 = ( (vE+vA+vB+vC)*vol_EABC + 
                             (vC+vD+vE+vF)*vol_CDEF + 
                             (vA+vC+vE+vD)*vol_ACED );

        const Real int_5 = ( (vF+vA+vB+vC)*vol_FABC + 
                             (vA+vD+vE+vF)*vol_ADEF + 
                             (vA+vB+vE+vF)*vol_ABEF );

        const Real int_6 = ( (vF+vA+vB+vC)*vol_FABC + 
                             (vB+vD+vE+vF)*vol_BDEF + 
                             (vA+vB+vD+vF)*vol_ABDF );

        // integrals need to be scaled by 1/24, but also must scale average by 6
        result = (int_1 + int_2 + int_3 + int_4 + int_5 + int_6)/144.;
    }

    return result;
}
