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
                        Vector<MultiFab*>&  data,
                        const Vector<int>& comps,
                        Vector<int>&        faceData,
                        int&               nElts,
                        Vector<Vector<Vector<int> > >& inside_nodes,
                        const Vector<int>&  eltIDs);
void
write_ml_streamline_data(const std::string&       outfile,
                         const Vector<MultiFab*>& data,
                         const Vector<string>&    names,
                         const Vector<int>&       faceData,
                         int                      nElts,
                         const Vector<Vector<Vector<int> > >& inside_nodes);

std::vector<MLloc>
build_nodeMap(const Vector<Vector<Vector<int> > >& inside_nodes);

int get_nPts(const Vector<MultiFab*>& data);
int get_jlo(const Vector<MultiFab*>& data);

Real Vdot(const vector<Real>& A,
          const vector<Real>& B);

vector<Real> Vcross(const vector<Real>& A,
                    const vector<Real>& B);

vector<Real> Vminus(const vector<Real>& A,
                    const vector<Real>& B);

Vector<string>
read_streamline_names(const std::string& infile)
{
    // Read header for path traces
    const string Header = infile + string("/Header");
    std::ifstream ifh;
    ifh.open(Header.c_str());
    std::string typeLabel; ifh >> typeLabel;
    int Nlev; ifh >> Nlev;
    int nComp; ifh >> nComp;
    Vector<string> names(nComp);
    for (int i=0; i<nComp; ++i)
        ifh >> names[i];
    ifh.close();
    return names;
}

int
main (int   argc,
      char* argv[])
{
    Initialize(argc,argv);

    ParmParse pp;

    std::string infile; pp.get("infile",infile);

    // Build outfile name (may modify below if more than one dataset)
    vector<string> infileTokens = Tokenize(infile,".");
    string outfile = infileTokens[0];
    for (int i=1; i<infileTokens.size()-1; ++i)
        outfile += string(".") + infileTokens[i];
    outfile += "_new"; pp.query("outfile",outfile);

    int verbose=0; pp.query("verbose",verbose);

    Vector<MultiFab*> streamlines;
    Vector<int> faceData;
    int nElts;
    Vector<Vector<Vector<int> > > inside_nodes;
    Vector<int> strComps;

    Vector<int> eltIDs;
    if (int ne = pp.countval("eltIDs"))
    {
        eltIDs.resize(ne);
        pp.getarr("eltIDs",eltIDs,0,ne);
    }
    else
    {
        int sElt = 0;
        pp.query("sElt",sElt);
        int nElt = 1;
        pp.query("nElt",nElt);
        eltIDs.resize(nElt);
        for (int i=0; i<nElt; ++i)
            eltIDs[i] = sElt + i;
    }


    Vector<string> readNames;
    int nc= pp.countval("comps");
    if (nc > 0)
    {
        readNames.resize(nc);
        pp.getarr("comps",readNames,0,nc);
    }
    Vector<int> readComps(nc);

    Vector<string> strNames = read_streamline_names(infile);
    if (nc==0) {
        nc = strNames.size();
        readComps.resize(nc);
        readNames.resize(nc);
        for (int i=0; i<nc; ++i) {
            readComps[i] = i;
            readNames[i] = strNames[i];
        }
    } else {
        int cnt=0;
        for (int i=0; i<strNames.size(); ++i)
        {
            for (int j=0; j<readNames.size(); ++j)
                if (strNames[i] == readNames[j])
                    readComps[cnt++] = i;
        }
    }

    if (verbose)
    {
        std::cerr << "Reading components: ";
        for (int j=0; j<readComps.size(); ++j)
            std::cerr << readComps[j] << " ";
        std::cerr << " from stream file ...\n";
    }
    read_ml_streamline_data(infile,streamlines,readComps,faceData,nElts,inside_nodes,eltIDs);
    if (verbose)
        std::cerr << "...finished reading stream file \n";

    write_ml_streamline_data(outfile,streamlines,readNames,faceData,nElts,inside_nodes);

    Finalize();
    return 0;
}

int
get_jlo(const Vector<MultiFab*>& data)
{
    bool found_one_yet = false;
    int jlo=0;
    for (int i=0; i<data.size(); ++i)
    {
        for (int j=0; j<data[i]->size(); ++j)
        {
            if (data[i]->defined(j))
            {
                if (!found_one_yet)
                    jlo = (*data[i])[j].box().smallEnd(1);
                else
                    jlo = std::min(jlo,(*data[i])[j].box().smallEnd(1));
            }
        }
    }
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
    for (int i=0; i<Nlev; ++i) {
        for (int j=0; j<inside_nodes[i].size(); ++j) {
            for (int k=0; k<inside_nodes[i][j].size(); ++k) {
                if (inside_nodes[i][j][k]<=0 || inside_nodes[i][j][k]>num_nodes)
                    cout << "inside_nodes: " << inside_nodes[i][j][k] << '\n';
                nodeMap[inside_nodes[i][j][k]-1] = MLloc(i,j,k);
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

static
void
ReadMF (MultiFab&          mf,
        const Vector<int>& comps,
        const Vector<int>& boxIDs,
        const std::string& infile)
{
    VisMF vismf(infile);

    BL_ASSERT(comps.size()<=vismf.nComp());

    const BoxArray& ba = vismf.boxArray();

    BoxArray new_ba(boxIDs.size());
    for (int i=0; i<boxIDs.size(); ++i) {
        new_ba.set(i,ba[boxIDs[i]]);
    }
    DistributionMapping new_dm(new_ba);
    mf.define(new_ba,new_dm,comps.size(),vismf.nGrow());

    if (boxIDs.size()!=0)
    {
        for (int j=0; j<comps.size(); ++j)
        {
            BL_ASSERT(comps[j]>=0 && comps[j]<vismf.nComp());
        }
        
        for (MFIter mfi(mf); mfi.isValid(); ++mfi)
        {
            cerr << "Reading mf.  idx: " << boxIDs[mfi.index()] << endl;
            Box gbox = Box(mfi.validbox()).grow(vismf.nGrow());
            for (int j=0; j<comps.size(); ++j)
            {
                const FArrayBox& fab = vismf.GetFab(boxIDs[mfi.index()], comps[j]);
                mf[mfi].copy(fab,0,j,1);
                vismf.clear(boxIDs[mfi.index()], comps[j]);
            }
        }
    }
}

void
read_ml_streamline_data(const std::string& infile,
                        Vector<MultiFab*>& data,
                        const Vector<int>& comps,
                        Vector<int>&       faceData,
                        int&               nElts,
                        Vector<Vector<Vector<int> > >& inside_nodes,
                        const Vector<int>& eltIDs)
{
    // Read header for path traces
    const string Header = infile + string("/Header");
    std::ifstream ifh;
    ifh.open(Header.c_str());
    std::string typeLabel; ifh >> typeLabel;
    int NlevPath; ifh >> NlevPath;
    int nCompPath; ifh >> nCompPath;
    vector<string> namesPath(nCompPath);
    for (int i=0; i<nCompPath; ++i)
        ifh >> namesPath[i];
    ifh.close();

    // Read elements
    const string iElements = infile + string("/Elements");
    std::ifstream ife;
    ife.open(iElements.c_str());
    ife >> nElts;
    int nodesPerElt; ife >> nodesPerElt;
    Vector<int> fileFaceData(nElts*nodesPerElt);
    for (int i=0; i<nElts*nodesPerElt; ++i)
        ife >> fileFaceData[i];

    nElts = eltIDs.size();
    faceData.resize(nElts*nodesPerElt);
    for (int j=0; j<nElts; ++j) {
        for (int i=0; i<nodesPerElt; ++i) {
            faceData[i] = fileFaceData[eltIDs[j]*nodesPerElt+i];
        }
    }

    Vector<Vector<Vector<int> > > file_inside_nodes(NlevPath);
    for (int lev=0; lev<NlevPath; ++lev)
    {
        // Get the name of the mf file
        char buf[64];
        sprintf(buf, "/Level_%d", lev);
        std::string ThisStreamFile = infile + std::string(buf);
        ThisStreamFile += "/Str";

        VisMF vismf(ThisStreamFile);
        file_inside_nodes[lev].resize(vismf.boxArray().size());

        // Get number of nonempty boxes
        int num_non_zero;
        ife >> num_non_zero;

        for (int j=0; j<num_non_zero; ++j)
        {
            // Get index of non-empty box, and number of entries
            int box_id, num_ids;
            ife >> box_id >> num_ids;

            file_inside_nodes[lev][box_id].resize(num_ids);

            for (int k=0; k<num_ids; ++k)
                ife >> file_inside_nodes[lev][box_id][k];
        }
    }
    ife.close();

    // Find just the boxes we need to read
    std::vector<MLloc> nodeMap = build_nodeMap(file_inside_nodes);
    Vector<Vector<int> > boxIDs(NlevPath);
    Vector<set<int> > boxes_needed(NlevPath);

    for (int k=0; k<eltIDs.size(); ++k) {
        int eltID = eltIDs[k];
        int offset = eltID*nodesPerElt;
        vector<MLloc> p(nodesPerElt);
        for (int i=0; i<nodesPerElt; ++i) {
            const MLloc& pi = nodeMap[ faceData[offset+i] - 1 ];
            boxes_needed[pi.amr_lev].insert(pi.box_idx);
        }
    }
    for (int i=0; i<NlevPath; ++i) {
        boxIDs[i].resize(boxes_needed[i].size());
        int cnt = 0;
        for (set<int>::const_iterator sit=boxes_needed[i].begin(); sit!=boxes_needed[i].end(); ++sit) {
            boxIDs[i][cnt++] = *sit;
        }
    }

    data.resize(NlevPath);
    for (int lev=0; lev<NlevPath; ++lev)
    {
        char buf[64];
        sprintf(buf, "/Level_%d", lev);
        std::string ThisStreamFile = infile + std::string(buf);
        ThisStreamFile += "/Str";
        data[lev] = new MultiFab();
        cerr << "Possibly reading mf level " << lev << endl;
        ReadMF(*data[lev],comps,boxIDs[lev],ThisStreamFile);
    }

    // Build reduced inside_nodes, use reduced numbering

    int tot = 0;
    for (int lev=0; lev<NlevPath; ++lev) {
        for (int i=0; i<boxIDs[lev].size(); ++i) {
            tot += file_inside_nodes[lev][boxIDs[lev][i]].size();
        }
    }
    map<int,int> newNum;
    int cnt = 1; // 1-based
    for (int lev=0; lev<NlevPath; ++lev) {
        for (int i=0; i<boxIDs[lev].size(); ++i) {
            for (int j=0; j<file_inside_nodes[lev][boxIDs[lev][i]].size(); ++j) {
                newNum[file_inside_nodes[lev][boxIDs[lev][i]][j]] = cnt++;
            }
        }
    }
    
    inside_nodes.resize(NlevPath);
    for (int lev=0; lev<NlevPath; ++lev) {
        inside_nodes[lev].resize(boxIDs[lev].size());
        for (int i=0; i<boxIDs[lev].size(); ++i) {
            int size = file_inside_nodes[lev][boxIDs[lev][i]].size();
            inside_nodes[lev][i].resize(size);
            for (int j=0; j<size; ++j) {
                inside_nodes[lev][i][j] = newNum[ file_inside_nodes[lev][boxIDs[lev][i]][j] ];
            }
        }
    }
    // Reset elements
    for (int i=0; i<faceData.size(); ++i) {
        faceData[i] = newNum[faceData[i]];
    }

    //
    // Force other processors to wait until parallel data has been read
    //
    ParallelDescriptor::Barrier();
}

void
write_ml_streamline_data(const std::string&       outfile,
                         const Vector<MultiFab*>& data,
                         const Vector<string>&    names,
                         const Vector<int>&       faceData,
                         int                      nElts,
                         const Vector<Vector<Vector<int> > >& inside_nodes)
{
    if (ParallelDescriptor::IOProcessor())
    {
        if (!UtilCreateDirectory(outfile, 0755))
            CreateDirectoryFailed(outfile);

        // Write Header
        const string oHeader = outfile + string("/Header");
        std::ofstream ofh;
        ofh.open(oHeader.c_str());
        ofh << "Oddball-multilevel-connected-data-format" << '\n';
        ofh << data.size() << '\n';
        ofh << names.size() << '\n';
        for (int i=0; i<names.size(); ++i)
            ofh << names[i] << '\n';
        ofh.close();

        // Write elements
        int nodesPerElt = faceData.size() / nElts;
        const string oElements = outfile + string("/Elements");
        std::ofstream ofe;
        ofe.open(oElements.c_str());
        ofe << nElts << '\n';
        ofe << nodesPerElt << '\n';
        for (int i=0; i<faceData.size(); ++i)
            ofe << faceData[i] << " ";
        ofe << '\n';
        // Write element distribution
        for (int i=0; i<inside_nodes.size(); ++i)
        {
            const Vector<Vector<int> >& inside_nodes_i = inside_nodes[i];
 
            // Count the non-zero-length entries
            int num_non_zero = 0;
            for (int j=0; j<inside_nodes_i.size(); ++j)
                if (inside_nodes_i[j].size() > 0)
                    num_non_zero++;
 
            ofe << num_non_zero << '\n';
            for (int j=0; j<inside_nodes_i.size(); ++j)
            {
                const Vector<int>& inside_nodes_i_j = inside_nodes[i][j];
                if (inside_nodes_i[j].size() > 0)
                {
                    ofe << j << " " << inside_nodes_i_j.size();
                    for (int k=0; k<inside_nodes_i_j.size(); ++k)
                        ofe << " " << inside_nodes_i_j[k];
                    ofe << '\n';
                }
            }
        }
        ofe.close();
    }
    //
    // Force other processors to wait till directory is built and Header is written.
    //
    ParallelDescriptor::Barrier();

    // Write data
    for (int lev=0; lev<data.size(); ++lev)
    {
        char buf[64];
        sprintf(buf, "/Level_%d", lev);
        std::string DataFile = outfile + std::string(buf);
        
        if (ParallelDescriptor::IOProcessor())
            if (!UtilCreateDirectory(DataFile, 0755))
                CreateDirectoryFailed(DataFile);
        //
        // Force other processors to wait till directory is built.
        //
        ParallelDescriptor::Barrier();
        
        std::string FabDataFile = DataFile + "/Str";

        VisMF::Write(*data[lev],FabDataFile);
    }
}

