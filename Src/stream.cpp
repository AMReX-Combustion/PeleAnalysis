#include <string>
#include <iostream>
#include <set>

#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include <AMReX_DataServices.H>
#include <AMReX_BCRec.H>
#include <AMReX_Interpolater.H>
#include <AMReX_Extrapolater.H>
#include <AMReX_BLFort.H>

using namespace amrex;

extern "C" {
  void vtrace(const Real* T,    const int* T_lo,    const int* T_hi,    const int* nT,
              const Real* loc,  const int* loc_lo,  const int* loc_hi,  const int* nl,
              const int*  ids,  const int* n_ids,
              Real*       g,    const int* g_lo,    const int* g_hi,    const int* computeVec,
              Real*       strm, const int* strm_lo, const int* strm_hi, const int* ncs,
              const Real* dx, const Real* plo, const Real* hRK, const int* errFlag);

}
static const Real eps = 1.e-4; // distance inside domain surface points are pushed if they fall outside

struct MLloc
{
    MLloc() :
        amr_lev(-1), box_idx(-1), pt_idx(-1) {}
    MLloc(int lev,int box,int pt) :
        amr_lev(lev), box_idx(box), pt_idx(pt) {}
    int amr_lev, box_idx, pt_idx;
};

static
void 
print_usage (int,
             char* argv[])
{
    std::cerr << "usage:\n";
    std::cerr << argv[0] << " infile plotfile=<string> [options] \n\tOptions:\n";
    std::cerr << " isoFile=<string>  OR  seedLoc=<real real [real]\n";
    std::cerr << " streamFile=<string>  OR  outFile=<string>\n"; 
    std::cerr << " is_per=<int int int> (DEF=1 1 1)\n";
    std::cerr << " finestLevel=<int> (DEF=finest level in plotfile)\n";
    std::cerr << " progressName=<string> (DEF=temp)\n";
    std::cerr << " traceAlongV=<bool> (DEF=0)\n";
    std::cerr << " buildAltSurf=<bool> (DEF=0)\n";
    std::cerr << "     (if true, requires altVal=<real>, also takes dt=<real> (DEF=0) and altIsoFile=<string>)\n";
    std::cerr << " nRKsteps=<int> (DEF=51)\n";
    std::cerr << " hRK=<real> (DEF=.1 (*dx_finest in plotfile)\n";
    std::cerr << " nGrow=<int> (DEF=4)\n";
    std::cerr << " bounds=<float * 4> (DEF=NULL)\n";
    exit(1);
}

void
FillCFgrowCells(AmrData& amrd, const BoxArray& fine_ba,Vector<Geometry*>& geoms, int lev,
                MultiFab& crseSrc, int sComp, int nComp, int nGrow, MultiFab& growScr)
{
    int refRatio = amrd.RefRatio()[lev-1];
    //const BoxArray cfBoxesCRSE = CFBoxes(BoxArray(fine_ba).coarsen(refRatio),nGrow,*geoms[lev-1]);
    const BoxArray cfBoxesCRSE = BoxArray(GetBndryCells(BoxArray(fine_ba).coarsen(refRatio),nGrow));
    const DistributionMapping dmScrCRSE(cfBoxesCRSE);
    MultiFab growScrCRSE(cfBoxesCRSE,dmScrCRSE,nComp,0);

    crseSrc.FillBoundary(sComp,nComp,geoms[lev-1]->periodicity());
    
    BoxArray gcba = BoxArray(crseSrc.boxArray()).grow(nGrow);
    MultiFab scrCRSE(gcba,crseSrc.DistributionMap(),nComp,0); // Build growLess mf to get f-f and p-f grow data through parallel copy
    for (MFIter mfi(scrCRSE); mfi.isValid(); ++mfi)
        scrCRSE[mfi].copy(crseSrc[mfi],sComp,0,nComp);

    // Now we have good grow data in scrCRSE, copy to fill interp source data
    growScrCRSE.copy(scrCRSE,0,0,nComp);
    
    // Do interp from coarse data to fine
    const BoxArray cfBoxes = BoxArray(cfBoxesCRSE).refine(refRatio); // Get matching fine boxes for interp
    MultiFab growScrFC(cfBoxes,dmScrCRSE,nComp,0); // Container on matching fine boxes
    Vector<BCRec> bc(1); // unused for pc_interp...
    int idummy1=0, idummy2=0;
    for (MFIter mfi(growScrFC);mfi.isValid();++mfi) {
        pc_interp.interp(growScrCRSE[mfi],0,growScrFC[mfi],0,1,growScrFC[mfi].box(),
                         refRatio*IntVect::TheUnitVector(),*geoms[lev-1],*geoms[lev],bc,idummy1,idummy2,RunOn::Cpu);
    }    
    // Finally, build correct c-f boxes and copy-on-intersect from interp data
    const BoxArray cfBoxesFINE = BoxArray(GetBndryCells(fine_ba,nGrow));
    DistributionMapping dmBoxesFINE(cfBoxesFINE);
    growScr.define(cfBoxesFINE,dmBoxesFINE,nComp,0);
    growScr.copy(growScrFC); // Get c-f data into just-right size container
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

std::string
getFileRoot(const std::string& infile)
{
    vector<std::string> tokens = Tokenize(infile,std::string("/"));
    return tokens[tokens.size()-1];
}

void
push_nodes_inside(FArrayBox&          nodes,
                  const Vector<Real>& plo,
                  const Vector<Real>& phi,
                  const Real          epsPush)
{
    const int nNodes = nodes.box().numPts();
    for (int d=0; d<BL_SPACEDIM; ++d)
    {
        Real* x = nodes.dataPtr(d);

        for (int i=0; i<nNodes; ++i)
        {
            x[i] = std::max(plo[d]+epsPush, std::min(phi[d]-epsPush, x[i]));
        }
    }        
}                  

Vector<int>
setInsideNodes(const FArrayBox&                         nodes,
               const Box&                               box,
               const std::vector< std::pair<int,Box> >& fboxes,
               const Vector<Real>&                       dx,
               const Vector<Real>&                       plo)
{
    const Box& nbox = nodes.box();

    vector<Real> lo(BL_SPACEDIM);
    vector<Real> hi(BL_SPACEDIM);
    for (int i=0; i<BL_SPACEDIM; ++i)
    {
        lo[i] = plo[i] + (box.loVect()[i]     )*dx[i];
        hi[i] = plo[i] + (box.hiVect()[i] + 1.)*dx[i];
    }

    // Get bounds of coarsened fine boxes as well
    const int Nfboxes = fboxes.size();
    vector<vector<Real> > flo(Nfboxes);
    vector<vector<Real> > fhi(Nfboxes);
    for (int n=0; n<Nfboxes; ++n)
    {
        flo[n].resize(BL_SPACEDIM);
        fhi[n].resize(BL_SPACEDIM);

        for (int i=0; i<BL_SPACEDIM; ++i)
        {
            flo[n][i] = plo[i] + (fboxes[n].second.loVect()[i]     )*dx[i];
            fhi[n][i] = plo[i] + (fboxes[n].second.hiVect()[i] + 1.)*dx[i];
        }
    }

    Vector<const Real*> loc(BL_SPACEDIM);
    for (int i=0; i<BL_SPACEDIM; ++i)
        loc[i] = nodes.dataPtr(i);

    std::list<int> id_list; // id of nodes inside, to be passed back in an array of ints made afterward
    int cnt = 0;
    for (IntVect iv=nbox.smallEnd(); iv<=nbox.bigEnd(); nbox.next(iv), cnt++)
    {
        bool isIn = true;

        for (int i=0; i<BL_SPACEDIM; ++i)
        {
            isIn &= (loc[i][cnt]>=lo[i] && loc[i][cnt]<hi[i]);
        }

        if (isIn)
        {
            // Now see if this spot is convered by the fine grid
            for (int n=0; n<Nfboxes && isIn; ++n)
            {
                bool inThisBox = true;
                for (int i=0; i<BL_SPACEDIM; ++i)
                {
                    inThisBox &= (loc[i][cnt]>=flo[n][i] && loc[i][cnt]<fhi[n][i]);
                }

                isIn = !(inThisBox);  // Point in finer box, dont add to this guys list
            }

            if (isIn)
                id_list.push_back(cnt);
        }
    }
    // Make a vector of ids for easier access later, make this structure "1"-based
    Vector<int> inside_nodes(id_list.size());
    cnt = 0;
    for (std::list<int>::const_iterator it=id_list.begin(); it!=id_list.end(); ++it)
    {
        inside_nodes[cnt] = *it + 1;
        cnt++;
    }
    return inside_nodes;
}

void
trim_surface(const vector<Real>& bbll,
             const vector<Real>& bbur,
             FArrayBox&          nodes,
             Vector<int>&         faceData,
             int                 nodesPerElt)
{
    const Box& nbox = nodes.box();
    const Real* xp = nodes.dataPtr(0);
    const Real* yp = nodes.dataPtr(1);
    const Real* zp = nodes.dataPtr(2);

    int cnt = 0;
    int cnt_new = 0;
    vector<int> nodeIdx;
    for (IntVect iv=nbox.smallEnd(); iv<=nbox.bigEnd(); nbox.next(iv), cnt++)
    {
        const Real x=xp[cnt];
        const Real y=yp[cnt];
        const Real z=zp[cnt];
        bool remove_this_node = x<bbll[0] || x>bbur[0] || y<bbll[1] || y>bbur[1] || z<bbll[2] || z>bbur[2];
        
        nodeIdx.push_back( remove_this_node ? -1 : cnt_new++ );
    }
    int nNodesNEW = cnt_new;

    // Build new nodes structure with bad points removed
    Box newBox(IntVect::TheZeroVector(),(nNodesNEW-1)*amrex::BASISV(0));
    int nNodesOLD = nbox.numPts();
    int nComp = nodes.nComp();
    FArrayBox newNodes(newBox,nComp);
    Vector<Real*> odp(nComp);
    Vector<Real*> ndp(nComp);
    for (int j=0; j<nComp; ++j)
    {
        odp[j] = nodes.dataPtr(j);
        ndp[j] = newNodes.dataPtr(j);
    }
    cnt_new = 0;
    for (int i=0; i<nNodesOLD; ++i)
    {
        int newNode = nodeIdx[i];
        if (newNode>=0)
            for (int j=0; j<nComp; ++j)
                ndp[j][newNode] = odp[j][i];
    }
    nodes.resize(newBox,nComp);
    nodes.copy(newNodes);
    newNodes.clear(); // Make some space...not really necessary, but what the heck

    const int nElts = faceData.size() / nodesPerElt;
    BL_ASSERT(nElts * nodesPerElt == faceData.size()); // Idiot check

    // Remove elements that refer to removed nodes
    vector<int> newFaceData;
    cnt_new = 0;
    for (int i=0; i<nElts; ++i)
    {
        int offset = nodesPerElt * i;
        bool eltGood = true;
        for (int j=0; j<nodesPerElt; ++j)
            eltGood &= (nodeIdx[ faceData[offset+j] - 1 ] >= 0); // Remember that faceData is 1-based

        if (eltGood)
        {
            for (int j=0; j<nodesPerElt; ++j)
                newFaceData.push_back(nodeIdx[ faceData[offset+j] - 1 ] + 1); // Make sure new faceData is 1-based
            cnt_new++;
        }
    }

    faceData.resize(cnt_new*nodesPerElt);
    for (int i=0; i<faceData.size(); ++i)
        faceData[i] = newFaceData[i];
}

void
build_surface_at_isoVal(const Vector<MultiFab*>&   paths,
                        const vector<string>&     pathNames,
                        int                       xComp,
                        int                       isoComp,
                        Real                      isoVal,
                        const Vector<Vector<Vector<int> > >& inside_nodes,
                        const Vector<int>&         comps,
                        FArrayBox&                altSurfNodes,
                        vector<string>&           altSurfNames,
                        string&                   distanceVarName);

void
add_thermal_thickness_to_surf(const Vector<MultiFab*>&   paths,
                              const vector<string>&     pathNames,
                              int                       xComp,
                              const Vector<Vector<Vector<int> > >& inside_nodes,
                              FArrayBox&                altSurfNodes,
                              vector<string>&           altSurfNames,
                              const string&             thickName,
                              int                       thickComp,
                              Real                      loVal,
                              Real                      hiVal);

void
add_cold_strain_to_surf(const Vector<MultiFab*>&   paths,
                        const vector<string>&     pathNames,
                        int                       xComp,
                        const Vector<Vector<Vector<int> > >& inside_nodes,
                        FArrayBox&                altSurfNodes,
                        vector<string>&           altSurfNames,
                        const string&             coldStrainName,
                        int                       strainComp,
                        int                       TComp,
                        Real                      TVal);

void
add_angle_to_surf(const Vector<MultiFab*>&   paths,
                  const vector<string>&     pathNames,
                  int                       xComp,
                  const Vector<Vector<Vector<int> > >& inside_nodes,
                  FArrayBox&                altSurfNodes,
                  vector<string>&           altSurfNames,
                  const string&             angleName);

void
write_ml_streamline_data(const std::string&       outfile,
                         const Vector<MultiFab*>& data,
                         int                      sComp,
                         const vector<string>&    names,
                         const Vector<int>&       faceData,
                         int                      nElts,
                         const Vector<Vector<Vector<int> > >& inside_nodes,
                         const AmrData&           amrdata);

void
dump_ml_streamline_data(const std::string&       outfile,
                        const Vector<MultiFab*>& data,
                        int                      sComp,
                        const vector<string>&    names,
                        const Vector<int>&       faceData,
                        int                      nElts,
                        const Vector<Vector<Vector<int> > >& inside_nodes);

void
write_iso(const std::string&    outfile,
          const FArrayBox&      nodes,
          const Vector<int>&     faceData,
          int                   nElts,
          const vector<string>& names,
          const string&         label);

std::string box_string(const Box& box) {
    std::string res;
    const IntVect& se = box.smallEnd();
    const IntVect& be = box.bigEnd();
    Vector<std::string> ses(BL_SPACEDIM), bes(BL_SPACEDIM);
    for (int d=0; d<BL_SPACEDIM; ++d) {
        std::string sroot = ( se[d] < 0 ? "n" : "");
        std::string broot = ( be[d] < 0 ? "n" : "");
        ses[d] = amrex::Concatenate(sroot,std::abs(se[d]),0);
        bes[d] = amrex::Concatenate(broot,std::abs(be[d]),0);
    }
    res = ses[0];
    for (int d=1; d<BL_SPACEDIM; ++d) res += "_" + ses[d];
    res += "__" + bes[0];
    for (int d=1; d<BL_SPACEDIM; ++d) res += "_" + bes[d];
    return res;
}

int
main (int   argc,
      char* argv[])
{
    amrex::Initialize(argc,argv);
    {
    Real strt_io, io_time = 0;

    if (argc < 2)
        print_usage(argc,argv);

    ParmParse pp;

    if (pp.contains("help"))
        print_usage(argc,argv);

    int verbose=0; pp.query("verbose",verbose);
    if (verbose>1)
        AmrData::SetVerbose(true);

    std::string plotfile; pp.get("plotfile",plotfile);

    strt_io = ParallelDescriptor::second();
    DataServices::SetBatchMode();
    Amrvis::FileType fileType(Amrvis::NEWPLT);

    DataServices dataServices(plotfile, fileType);
    if( ! dataServices.AmrDataOk()) {
        DataServices::Dispatch(DataServices::ExitRequest, NULL);
        // ^^^ this calls ParallelDescriptor::EndParallel() and exit()
    }
    AmrData& amrData = dataServices.AmrDataRef();
    io_time += ParallelDescriptor::second() - strt_io;

    int finestLevel = amrData.FinestLevel();
    pp.query("finestLevel",finestLevel);
    int Nlev = finestLevel + 1;

    std::string progressName = "temp"; pp.query("progressName",progressName);

    int idC = -1;
    const Vector<std::string>& plotVarNames = amrData.PlotVarNames();
    for (int i=0; i<plotVarNames.size(); ++i)
    {
        if (plotVarNames[i] == progressName) idC = i;
    }
    if (ParallelDescriptor::IOProcessor() && idC<0)
        std::cerr << "Cannot find required data in pltfile" << std::endl;

    // Read in seed pt locations
    FArrayBox nodes;
    Vector<int> faceData;
    int nElts;
    int nodesPerElt;
    int nSeedNodes, nCompSeedNodes;
    std::vector<std::string> surfNames;
    if (pp.countval("isoFile")>0)
    {
      // Read in isosurface
      std::string isoFile; pp.get("isoFile",isoFile);
      if (ParallelDescriptor::IOProcessor())
        std::cerr << "Reading isoFile... " << isoFile << std::endl;
    
      strt_io = ParallelDescriptor::second();

      std::ifstream ifs;
      ifs.open(isoFile.c_str(),std::ios::in|std::ios::binary);
      const std::string title = parseTitle(ifs);
      surfNames = parseVarNames(ifs);
      nCompSeedNodes = surfNames.size();

      ifs >> nElts;
      ifs >> nodesPerElt;

      FArrayBox tnodes;
      tnodes.readFrom(ifs);
      nSeedNodes = tnodes.box().numPts();

      // transpose the data so that the components are 'in the right spot for fab data'
      nodes.resize(tnodes.box(),nCompSeedNodes);
      Real** np = new Real*[nCompSeedNodes];
      for (int j=0; j<nCompSeedNodes; ++j)
        np[j] = nodes.dataPtr(j);

      Real* ndat = tnodes.dataPtr();
      for (int i=0; i<nSeedNodes; ++i)
      {
        for (int j=0; j<nCompSeedNodes; ++j)
        {
          np[j][i] = ndat[j];
        }
        ndat += nCompSeedNodes;
      }
      delete [] np;
      tnodes.clear();

      faceData.resize(nElts*nodesPerElt,0);
      ifs.read((char*)faceData.dataPtr(),sizeof(int)*faceData.size());
      ifs.close();
    }
    else
    {
      nSeedNodes = 1;
      nodes.resize(Box(IntVect::TheZeroVector(),IntVect::TheZeroVector()),BL_SPACEDIM);
      nodesPerElt = 1;
      nElts = 1;
      faceData.resize(nElts*nodesPerElt,1);
      surfNames = {D_DECL("X", "Y", "Z")};
      Vector<Real> loc(BL_SPACEDIM); 
      pp.getarr("seedLoc",loc,0,BL_SPACEDIM);
      for (int i=0; i<BL_SPACEDIM; ++i) nodes(IntVect::TheZeroVector(),i) = loc[i];
    }
    io_time += ParallelDescriptor::second() - strt_io;
    if (ParallelDescriptor::IOProcessor() && verbose>0)
        std::cerr << "...finished reading isoFile " << std::endl;

    // Enforce that all nodes be inside bounds of plotfile domain
    const Vector<Real>& plo = amrData.ProbLo();
    const Vector<Real>& phi = amrData.ProbHi();
    Real epsilon = eps * amrData.DxLevel()[finestLevel][0];
    push_nodes_inside(nodes,plo,phi,epsilon);

    if (int nx=pp.countval("bounds"))
    {
        Vector<Real> barr;
        pp.getarr("bounds",barr,0,nx);
        int d=BL_SPACEDIM;
        BL_ASSERT(barr.size()==2*d);
        vector<Real> bbll(BL_SPACEDIM);
        vector<Real> bbur(BL_SPACEDIM);
        for (int i=0; i<d; ++i)
        {
            bbll[i] = barr[i];
            bbur[i] = barr[d+i];
        }
        if (ParallelDescriptor::IOProcessor() && verbose>0)
            std::cerr << "...limiting stream seed points to physical range = ("
                 << bbll[0] << "," << bbll[1] << "," << bbll[2] << ") to ("
                 << bbur[0] << "," << bbur[1] << "," << bbur[2] << ")" << std::endl;

        trim_surface(bbll,bbur,nodes,faceData,nodesPerElt);
        nSeedNodes = nodes.box().numPts();
        nElts = faceData.size() / nodesPerElt;
        BL_ASSERT(nElts*nodesPerElt == faceData.size());

    }


    // Process pltfile
    bool traceAlongV=false; pp.query("traceAlongV",traceAlongV);
    bool needV=traceAlongV;
    Real altVal=-1;
    Real dt=0.; // dt is time increment to advect alt surf, if built

    int xCompStr = 0;
    int isoCompSt = 0;
    int isoCompStr = xCompStr + BL_SPACEDIM;
    int vCompSt = -1;
    int vCompStr = -1;
    int nCompStr = isoCompStr + 1;
    int nCompSt = isoCompSt + 1;

    bool buildAltSurf=false; pp.query("buildAltSurf",buildAltSurf);
    string thickCompName="null";
    string strainCompName="null";
    string TCompName="null";
    Real thickLo, thickHi, TVal;
    bool addAngle=false;
    if (buildAltSurf)
    {
        pp.get("altVal",altVal); // This is the value of progress searched for to find the "alt" velocity
        pp.query("dt",dt);
        needV=true; // Might need velocities at altVal (if going to advect the surface)

        pp.query("thickCompName",thickCompName);
        if (thickCompName!="null")
        {
            pp.get("thickLo",thickLo);
            pp.get("thickHi",thickHi);
        }
        
        pp.query("strainCompName",strainCompName);
        if (strainCompName!="null")
        {
            pp.get("TCompName",TCompName);
            pp.get("TVal",TVal);
        }

        pp.query("addAngle",addAngle);
    }
    if (needV)
    {
        vCompSt  = isoCompSt + 1;
        vCompStr = isoCompStr + 1;
        nCompStr = vCompStr + BL_SPACEDIM;
        nCompSt = vCompSt + BL_SPACEDIM;
    }

    Vector<std::string> inVarNames(nCompSt);
    inVarNames[isoCompSt] = progressName;
    if (needV)
    {
        inVarNames[vCompSt  ] = "x_velocity";
        inVarNames[vCompSt+1] = "y_velocity";
#if BL_SPACEDIM==3
        inVarNames[vCompSt+2] = "z_velocity";
#endif
    }
    vector<string> strNames(nCompStr);
    for (int i=0; i<BL_SPACEDIM; ++i)
        strNames[i] = surfNames[i];
    for (int i=0; i<nCompSt; ++i)
        strNames[BL_SPACEDIM+i] = inVarNames[i];

    Vector<int> destFillComps(nCompSt);
    for (int i=0; i<nCompSt; ++i)
        destFillComps[i] = i;

    // Grab auxiliary variables from plotfile
    Vector<int> auxComps;
    if (int nc = pp.countval("aux_comps"))
    {
        auxComps.resize(nc);
        pp.getarr("aux_comps",auxComps,0,nc);
    }
    else
    {
        int sComp = 0;
        pp.query("aux_sComp",sComp);
        int nComp = 0;
        pp.query("aux_nComp",nComp);
        BL_ASSERT(sComp+nComp <= amrData.NComp());
        auxComps.resize(nComp);
        for (int i=0; i<nComp; ++i)
            auxComps[i] = sComp + i;
    }
    if (auxComps.size()>0)
    {
        inVarNames.resize(inVarNames.size()+auxComps.size());
        destFillComps.resize(destFillComps.size()+auxComps.size());
        strNames.resize(strNames.size()+auxComps.size());
        for (int i=0; i<auxComps.size(); ++i)
        {
            BL_ASSERT(auxComps[i]<amrData.NComp());
            inVarNames[nCompSt+i] = plotVarNames[auxComps[i]];
            strNames[nCompStr+i] = plotVarNames[auxComps[i]];
            destFillComps[nCompSt+i] = nCompSt+i;
        }
        nCompSt += auxComps.size();
        nCompStr += auxComps.size();
    }

    // Find various names in strNames
    int thickComp =-1;
    int strainComp =-1;
    int TComp=-1;
    for (int i=0; i<strNames.size(); ++i)
    {
        if (strNames[i]==thickCompName)
            thickComp=i;

        if (strNames[i]==strainCompName)
            strainComp=i;

        if (strNames[i]==TCompName)
            TComp=i;
    }

    if (ParallelDescriptor::IOProcessor())
    {
        if (thickComp>=0)
            std::cout << "Thermal thickness component id: " << thickComp << std::endl;
        if (strainComp>=0)
        {
            std::cout << "Strain component id: " << strainComp << std::endl;
            std::cout << "   T component for strain eval at id: " << TComp << std::endl;
        }
    }

    // Stream trace parameters
    int nRKsteps = 51; pp.query("nRKsteps",nRKsteps);
    int nRKh = (nRKsteps-1)/2;
    Real hRK = .1; pp.query("hRK",hRK); // Fraction of finest grid spacing to step in RK stream integration
    //hRK = hRK * amrData.DxLevel()[amrData.FinestLevel()][0];
    int nGrow = (int)(hRK*nRKh) + 2; pp.query("nGrow",nGrow);
    if (ParallelDescriptor::IOProcessor()) {
        std::cout << "nGrow = " << nGrow << std::endl;
    }
    RealBox rb(&(amrData.ProbLo()[0]),
               &(amrData.ProbHi()[0]));
    int coord = 0;
    Vector<int> is_per(BL_SPACEDIM,1);
    pp.queryarr("is_per",is_per,0,BL_SPACEDIM);
    Print() << "Periodicity assumed for this case: ";
    for (int i=0; i<BL_SPACEDIM; ++i) {
        Print() << is_per[i] << " ";
    }
    Print() << std::endl;

    hRK = hRK * amrData.ProbSize()[0] / amrData.ProbDomain()[finestLevel].length(0);

    // Data structure holding node id of each path line (globally visible)
    Vector<Vector<Vector<int> > > inside_nodes(Nlev);

    vector<BoxArray> bal_str(Nlev); // The "BoxArray" to hold the streamlines
    Vector<Real> delta(BL_SPACEDIM);
    for (int lev=0; lev<Nlev; ++lev)
    {
        const BoxArray ba = amrData.boxArray(lev);
        //const Vector<Real>& delta = amrData.DxLevel()[lev];
        // We do not trust the dx in the plotfile, make our own
        for (int i=0; i<BL_SPACEDIM; ++i)
            delta[i] = amrData.ProbSize()[i] / amrData.ProbDomain()[lev].length(i);

        BoxArray baf_c;
        if (lev < finestLevel)
            baf_c = BoxArray(amrData.boxArray(lev+1)).coarsen(amrData.RefRatio()[lev]);

        BoxList bl(IndexType(D_DECL(IndexType::CELL,
                                    IndexType::CELL,
                                    IndexType::CELL)));

        inside_nodes[lev].resize(ba.size());
        for (int i=0; i<ba.size(); ++i)
        {
            const Box& box = ba[i];

            // Get list of points inside so that we can allocate space to hold resulting traces
            std::vector< std::pair<int,Box> > isects;
            if (lev < finestLevel)
                isects = baf_c.intersections(box);
            
            inside_nodes[lev][i] = setInsideNodes(nodes,box,isects,delta,plo);
            const int num_inside = inside_nodes[lev][i].size();

            if (num_inside == 0)
            {
                bl.push_back(Box(IntVect::TheZeroVector(),IntVect::TheZeroVector()));
            }
            else
            {
                bl.push_back(Box(IntVect(D_DECL(0,-nRKh,0)),IntVect(D_DECL(num_inside-1,-nRKh+nRKsteps-1,0))));
            }
        }

        // Because this box array is created to be the same length as the level box array
        //  is is guaranteed to have the same distribution mapping (so the next part will work)
        bal_str[lev] = BoxArray(bl);
    }

    Vector<MultiFab*> state(Nlev);
    Vector<Geometry*> geoms(Nlev);
    Vector<MultiFab*> streamlines(Nlev);

    for (int lev=0; lev<Nlev; ++lev)
    {
        const BoxArray ba = amrData.boxArray(lev);
        DistributionMapping dm(ba);
        state[lev] = new MultiFab(ba,dm,nCompSt,nGrow);
        //const Vector<Real>& delta = amrData.DxLevel()[lev];
        const Box& probDomain = amrData.ProbDomain()[lev];
        for (int i=0; i<BL_SPACEDIM; ++i)
            delta[i] = amrData.ProbSize()[i] / probDomain.length(i);

        if (ParallelDescriptor::IOProcessor())
            std::cerr << "Reading data for level " << lev << std::endl;

        if (ParallelDescriptor::IOProcessor() && lev==0)
        {
            std::cerr << "InVarNames:" << '\n';
            for (int i=0; i<inVarNames.size(); ++i)
                std::cerr << " " << inVarNames[i];
            std::cerr << '\n';
        }

        strt_io = ParallelDescriptor::second();
        amrData.FillVar(*state[lev],lev,inVarNames,destFillComps);
        
        for (int i=0; i<inVarNames.size(); ++i)
            amrData.FlushGrids(amrData.StateNumber(inVarNames[i]));
        io_time += ParallelDescriptor::second() - strt_io;

        Print() << "Data has been read for level " << lev << " ("
                << destFillComps.size() << " components)" << std::endl;

        //const bool do_corners = true;
        geoms[lev] = new Geometry(amrData.ProbDomain()[lev],&rb,coord,&(is_per[0]));

        MultiFab bigMF(BoxArray(ba).grow(nGrow),dm,nCompSt,0); // handy structure to simplify parallel copies

        // Fill progress variable grow cells using interp over c-f boundaries
        {
            Extrapolater::FirstOrderExtrap(*state[lev],*geoms[lev],isoCompSt,nCompSt);
            state[lev]->FillBoundary(isoCompSt,nCompSt,geoms[lev]->periodicity());
            
            if (lev>0)
            {
                MultiFab growScr; // Container in c-f grow cells

                FillCFgrowCells(amrData,ba,geoms,lev,*state[lev-1],isoCompSt,nCompSt,nGrow,growScr);

                // Get data into grow-free mf to allow parallel copy for grabbing c-f data
                for (MFIter mfi(*state[lev]); mfi.isValid(); ++mfi)
                    bigMF[mfi].copy((*state[lev])[mfi],isoCompSt,0,nCompSt); // get valid data, and extrap grow data
                
                bigMF.copy(growScr,0,0,nCompSt); // Overwrite c-f with preferred c-f interp data
                
                for (MFIter mfi(bigMF); mfi.isValid(); ++mfi)
                    (*state[lev])[mfi].copy(bigMF[mfi],0,isoCompSt,nCompSt); // Put result back into isoCompSt
            }
        }

        // Fix up fine-fine and periodic for idSmProg
        state[lev]->FillBoundary(isoCompSt,nCompSt,geoms[lev]->periodicity());

        if (ParallelDescriptor::IOProcessor())
            std::cerr << "Progress variable grow cells filled on level " << lev << std::endl;

        streamlines[lev] = new MultiFab(bal_str[lev],dm,nCompStr,0);

        FArrayBox tmp;
        for (MFIter mfi(*state[lev]); mfi.isValid(); ++mfi)
        {
            const FArrayBox& fab = (*state[lev])[mfi];
            FArrayBox& strm = (*streamlines[lev])[mfi];
            const Vector<int>& ids = inside_nodes[lev][mfi.index()];
            int n_ids = ids.size();

            if (n_ids==0)
            {
                strm.setVal(0.0);
            }
            else
            {
                const Box gbox = Box(fab.box()).grow(-1);                
                FArrayBox* vec;
                int computeVec, vecComp;
                if (traceAlongV)
                {
                    vec = (FArrayBox*)(&fab); // Cast away constness....to satisfy arg, promise not to change...
                    computeVec = 0;
                    vecComp = vCompSt;
                }
                else
                {
                    vec = &tmp;
                    tmp.resize(gbox,BL_SPACEDIM);
                    computeVec = 1;
                    vecComp = 0;
                }

                int errFlag = 0;

                vtrace(BL_TO_FORTRAN_N_ANYD(fab,isoCompSt),  &nCompSt,
                       BL_TO_FORTRAN_N_ANYD(nodes,      0),  &nCompSeedNodes,
                       ids.dataPtr(), &n_ids,
                       BL_TO_FORTRAN_N_ANYD(*vec, vecComp),  &computeVec,
                       BL_TO_FORTRAN_N_ANYD(strm,       0),  &nCompStr,
                       delta.dataPtr(), plo.dataPtr(), &hRK, &errFlag);
                         
                if (errFlag > 0) {
                    amrex::Abort("Problem with interpolation");
                }
            }
        }
        if (ParallelDescriptor::IOProcessor())
            std::cerr << "Streamlines computed on level " << lev << std::endl;

        if (ParallelDescriptor::IOProcessor())
            std::cerr << "Derive finished for level " << lev << std::endl;
    }

    if (pp.countval("streamFile") > 0)
    {
      ParallelDescriptor::Barrier();
      if (ParallelDescriptor::IOProcessor())
        std::cerr << "Writing the streamline data " << std::endl;

      // Get output filename
      std::string streamFile; pp.get("streamFile",streamFile);
      if (!streamFile.empty() && streamFile[streamFile.length()-1] != '/')
        streamFile += '/';
    
      int sComp = 0; // write stream lines ...
      strt_io = ParallelDescriptor::second();
      write_ml_streamline_data(streamFile,streamlines,sComp,strNames,faceData,nElts,inside_nodes,amrData);
      io_time += ParallelDescriptor::second() - strt_io;

      if (ParallelDescriptor::IOProcessor())
        std::cerr << "...done writing the streamline data " << std::endl;
    }
    else
    {
      std::string outFile;
      if (pp.countval("outFile")==0) Abort("Must specify streamFile or outFile");
      pp.get("outFile",outFile);
      dump_ml_streamline_data(outFile,streamlines,0,strNames,faceData,nElts,inside_nodes);
    }
    
    if (buildAltSurf)
    {
        FArrayBox altSurfNodes;
        vector<string> altSurfNames;
        Vector<int> comps;

        bool advectColdIso=false; pp.query("advectColdIso",advectColdIso);
        if (!advectColdIso)
        {
            // In this case, altVal is intended to be the "flame" and isoVal at the cold place
            //  We are not going to advect the flame, and we are going to need the distance
            //  separating the cold and hot isotherms.  The routine to build the alt surface
            //  will add a variable on the surface that is the distance from isoVal to altVal
            //  along the respective streamline.  We keep the isoComp, just in case the
            //  streamlines don't quite make it to altVal
            comps.resize(1);
            comps[0] = isoCompStr;
        }
        else
        {
            // In this case, altVal is intended to be the "cold" value, and we are making a surface
            //  in the cold region.  We will advect this surface, and so need to add the velocity
            //  to the quantities needed at the surface.  We also keep the isoComp, just in case
            //  the streamlines don't quite make it to altVal
            comps.resize(BL_SPACEDIM+1);
            for (int d=0; d<BL_SPACEDIM; ++d)
                comps[d] = vCompStr + d;
            comps[BL_SPACEDIM] = isoCompStr;
        }
        
        // Extract a new isosurface at altVal (hold in the form of a FAB, whose component names are "altSurfNames")
        // NOTE: only IOProc actually holds the surface data
        if (ParallelDescriptor::IOProcessor() && verbose>0)
            std::cout << "building new surface where " << strNames[isoCompStr] << " = "  << altVal << std::endl;

        string distanceVarName("distance_iso_to_alt");
        build_surface_at_isoVal(streamlines,strNames,xCompStr,isoCompStr,altVal,inside_nodes,
                                comps,altSurfNodes,altSurfNames,distanceVarName);

        // Add thermal thickness to alt surface
        if (thickComp>=0)
	{
            string thickName("thermalThickness");
            if (!advectColdIso)
                thickName = "thermalThickness_notAdv";
            add_thermal_thickness_to_surf(streamlines,strNames,xCompStr,inside_nodes,
                                          altSurfNodes,altSurfNames,thickName,thickComp,thickLo,thickHi);
        }

        // Add cold strain to alt surface
        if (TComp>=0)
	{
            string coldStrainName("coldStrain");
            add_cold_strain_to_surf(streamlines,strNames,xCompStr,inside_nodes,
                                    altSurfNodes,altSurfNames,coldStrainName,strainComp,TComp,TVal);
        }

        if (addAngle)
        {
            string angleName("angleWRTvert");
            add_angle_to_surf(streamlines,strNames,xCompStr,inside_nodes,
                              altSurfNodes,altSurfNames,angleName);
        }

        // Advect and write the surface...only I/O proc can do this
        if (ParallelDescriptor::IOProcessor())
        {
            // Advect the surface
            if (!advectColdIso)
            {
                // Grab the distance function from the original alt surface
                //  (which has the same connectivity as the new one) and the
                //  present one, and overwrite the latter with their difference.
                int dCompIso = -1;
                for (int i=0; i<surfNames.size(); ++i)
                    if (surfNames[i]==distanceVarName)
                        dCompIso = i;
                int dCompAlt = -1;
                for (int i=0; i<altSurfNames.size(); ++i)
                    if (altSurfNames[i]==distanceVarName)
                        dCompAlt = i;

                if ((dCompAlt>=0) && (dCompIso>=0))
                    altSurfNodes.plus(nodes,dCompIso,dCompAlt,1);

                // Rename this variable
                altSurfNames[dCompAlt] = "delta";
            }
            else
            {
                // Find the location, velocity
                int xCompSurf=-1, vCompSurf=xCompSurf;
                int nCompAltSurf = altSurfNodes.nComp();
                for (int i=0; i<nCompAltSurf; ++i)
                {
                    if (altSurfNames[i]==strNames[xCompStr])
                        xCompSurf=i;
                    if (altSurfNames[i]==strNames[vCompStr])
                        vCompSurf=i;
                }
                
                int nAltNodes = altSurfNodes.box().numPts();
                Real** csd = new Real*[nCompAltSurf];
                for (int d=0; d<nCompAltSurf; ++d)
                    csd[d] = altSurfNodes.dataPtr(d);
                
                // Move the points
                for (int i=0; i<nAltNodes; ++i)
                    for (int d=0; d<BL_SPACEDIM; ++d)
                        csd[xCompSurf+d][i] += csd[vCompSurf+d][i] * dt;

            }

            // Write the alt surface file
            strt_io = ParallelDescriptor::second();
            string altIsoFile = "surf";
            string label;
            if (advectColdIso)
            {
                altIsoFile += "_alt.mef";
                label="advected alt surface";
            }
            else
            {
                altIsoFile += "_new_flame.mef";
                label="new flame surface from advected alt";                
            }
            pp.query("altIsoFile",altIsoFile);
            write_iso(altIsoFile,altSurfNodes,faceData,nElts,altSurfNames,label);
            io_time += ParallelDescriptor::second() - strt_io;
        }

        if (ParallelDescriptor::IOProcessor())
            std::cerr << "Finished writing streamline data " << std::endl;
    }

    long nInsideNodes = 0;
    for (int i=0; i<inside_nodes.size(); ++i) {
        for (int j=0; j<inside_nodes[i].size(); ++j) {
            nInsideNodes += inside_nodes[i][j].size();
        }
    }
#if 0
    long min_nonfab_bytes = 
        faceData.size() * sizeof(int) 
        + nInsideNodes * sizeof(int);
    long max_nonfab_bytes = min_nonfab_bytes;
    const int IOProc = ParallelDescriptor::IOProcessorNumber();
    ParallelDescriptor::ReduceLongMin(min_nonfab_bytes,IOProc);
    ParallelDescriptor::ReduceLongMax(max_nonfab_bytes,IOProc);
    if (ParallelDescriptor::IOProcessor()) {
        std::cout << "\nNon-FAB byte spread across MPI nodes: ["
                  << min_nonfab_bytes
                  << " ... "
                  << max_nonfab_bytes
                  << "]\n";
    }

    long min_fab_bytes = amrex::total_bytes_allocated_in_fabs_hwm;
    long max_fab_bytes = amrex::total_bytes_allocated_in_fabs_hwm;
    
    ParallelDescriptor::ReduceLongMin(min_fab_bytes,IOProc);
    ParallelDescriptor::ReduceLongMax(max_fab_bytes,IOProc);
    if (ParallelDescriptor::IOProcessor()) {
        std::cout << "\nFAB byte spread across MPI nodes: ["
                  << min_fab_bytes
                  << " ... "
                  << max_fab_bytes
                  << "]\n";
    }

    const Real end_time = ParallelDescriptor::second();
    const Real run_time = end_time - strt_time - io_time;

    Real run_time_max, run_time_min; run_time_max = run_time_min = run_time;
    Real io_time_max, io_time_min; io_time_max = io_time_min = io_time;

    ParallelDescriptor::ReduceRealMax(run_time_max,IOProc);
    ParallelDescriptor::ReduceRealMax(io_time_max,IOProc);
    ParallelDescriptor::ReduceRealMin(run_time_min,IOProc);
    ParallelDescriptor::ReduceRealMin(io_time_min,IOProc);
    if (ParallelDescriptor::IOProcessor()) {
        std::cout << "Max Compute time: " << run_time_max << '\n';
        std::cout << "Min Compute time: " << run_time_min << '\n';
        std::cout << "Max I/O time: " << io_time_max << '\n';
        std::cout << "Min I/O time: " << io_time_min << '\n';
    }
#endif
    }
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
    ofs << label << std::endl;
    for (int i=0; i<nCompSurf; ++i)
    {
        ofs << names[i];
        if (i < nCompSurf-1)
            ofs << " ";
        else
            ofs << std::endl;
    }
    int nodesPerElt = faceData.size() / nElts;
    ofs << nElts << " " << nodesPerElt << std::endl;
    tnodes.writeOn(ofs);
    ofs.write((char*)faceData.dataPtr(),sizeof(int)*faceData.size());
    ofs.close();
}

void
add_angle_to_surf(const Vector<MultiFab*>&  paths,
                  const vector<string>&     pathNames,
                  int                       xComp,
                  const Vector<Vector<Vector<int> > >& inside_nodes,
                  FArrayBox&                altSurfNodes,
                  vector<string>&           altSurfNames,
                  const string&             angleName)
{
    int Nlev = paths.size();
    Box ZBOX(IntVect::TheZeroVector(),IntVect::TheZeroVector());

    int nCompAdd = 1;
    int nCompSurfOld = altSurfNodes.nComp();
    int nCompSurf = nCompSurfOld + nCompAdd;
    altSurfNames.resize(nCompSurf);
    altSurfNames[nCompSurfOld + 0] = angleName;

    int my_nNodes = 0;
    vector<Real> my_altSurfData, my_altSurfDatatmp;
    for (int lev = 0; lev<Nlev; ++lev)
    {
        for (MFIter mfi(*paths[lev]); mfi.isValid(); ++mfi)
        {
            const Box& box = mfi.validbox();
            const FArrayBox& pth = (*paths[lev])[mfi];
            
            if (box == ZBOX)
            {
                // Nothing to do
            }
            else
            {
                const int Npaths = box.length(0);
                my_nNodes += Npaths;

                // Use two points around the middle index of the line
                const int hpath = 0.5*(box.smallEnd(1) + box.bigEnd(1));
                const int loPath = hpath-1;
                const int hiPath = hpath+1;

                // For each path, find angle wrt vertical
                for (int i=0; i<Npaths; ++i)
                {
                    Real mag = 0;
                    Vector<Real> dx(BL_SPACEDIM);

                    for (int d=0; d<BL_SPACEDIM; ++d)
                    {
                        dx[d] = pth(IntVect(D_DECL(i,loPath,0)),xComp+d)
                            -   pth(IntVect(D_DECL(i,hiPath,0)),xComp+d);
                        mag += dx[d]*dx[d];
                    }
                    mag = std::sqrt(mag);
                    Real angle = std::acos(dx[2] / mag);

                    my_altSurfData.push_back(angle);
                }
            }
        }
    }

    if (ParallelDescriptor::IOProcessor())
        std::cerr << "Angle computed, communicating info to ioproc " << std::endl;

    // Communicate all data back to IOProc
    const int IOProc = ParallelDescriptor::IOProcessorNumber();        
    
    Vector<int> nNodes(ParallelDescriptor::NProcs(),0);
    Vector<int> offset(ParallelDescriptor::NProcs(),0);

#if BL_USE_MPI
    //
    // Tell root CPU how many tags each CPU will be sending.
    //
    MPI_Gather(&my_nNodes,
               1,
               ParallelDescriptor::Mpi_typemap<int>::type(),
               nNodes.dataPtr(),
               1,
               ParallelDescriptor::Mpi_typemap<int>::type(),
               IOProc,
               ParallelDescriptor::Communicator());
    
    if (ParallelDescriptor::IOProcessor())
    {
        for (int i = 0; i < nNodes.size(); i++)
            //
            // Convert from count of node data to add to count of reals to expect.
            //
            nNodes[i] *= nCompAdd;
        
        int N=0;
        for (int i = 0; i < nNodes.size(); i++)
            N += nNodes[i];

        my_altSurfDatatmp.resize(N);

        for (int i = 1; i < offset.size(); i++)
            offset[i] = offset[i-1] + nNodes[i-1];
    }
    //
    // Gather all the tags to IOProc into TheCollateSpace.
    //
    Real* bap = my_nNodes==0 ? 0 : &my_altSurfData[0];
    MPI_Gatherv(bap,
                my_nNodes*nCompAdd,
                ParallelDescriptor::Mpi_typemap<Real>::type(),
                (my_altSurfDatatmp.size()>0 ? &my_altSurfDatatmp[0] : 0),
                nNodes.dataPtr(),
                offset.dataPtr(),
                ParallelDescriptor::Mpi_typemap<Real>::type(),
                IOProc,
                ParallelDescriptor::Communicator());

    my_altSurfDatatmp.swap(my_altSurfData);

    if (ParallelDescriptor::IOProcessor())
        std::cerr << "Additional surface data sent" << std::endl;

#endif

    if (ParallelDescriptor::IOProcessor())
    {
        // Resize node fab to hold more data
        Box nBox = altSurfNodes.box();
        FArrayBox tmpFab(nBox,nCompSurfOld);
        tmpFab.copy(altSurfNodes,0,0,nCompSurfOld);
        altSurfNodes.resize(nBox,nCompSurf);
        altSurfNodes.copy(tmpFab,0,0,nCompSurfOld);

        Real** csd = new Real*[nCompAdd];
        for (int d=0; d<nCompAdd; ++d)
            csd[d] = altSurfNodes.dataPtr(nCompSurfOld + d);
        
        // De-scramble data so it looks like a isosurface node fab
        for (int lev=0; lev<Nlev; ++lev)
        {
            for (int i=0; i<paths[lev]->size(); ++i)
            {
                int proc = paths[lev]->DistributionMap()[i];
                int nPaths = inside_nodes[lev][i].size();
                for (int j=0; j<nPaths; ++j)
                {
                    int idx = inside_nodes[lev][i][j] - 1; // inside_nodes is 1-based

                    for (int d=0; d<nCompAdd; ++d)
                        csd[d][idx] = my_altSurfData[offset[proc] + d];

                    offset[proc]+= nCompAdd; // point to next box-ints on this proc
                }
            }
        }
        
        std::cerr << "New surface data added on ioproc " << std::endl;
    }
}

void
add_cold_strain_to_surf(const Vector<MultiFab*>&  paths,
                        const vector<string>&     pathNames,
                        int                       xComp,
                        const Vector<Vector<Vector<int> > >& inside_nodes,
                        FArrayBox&                altSurfNodes,
                        vector<string>&           altSurfNames,
                        const string&             coldStrainName,
                        int                       strainComp,
                        int                       TComp,
                        Real                      TVal)
{
    int Nlev = paths.size();
    Box ZBOX(IntVect::TheZeroVector(),IntVect::TheZeroVector());

    int nCompAdd = 1;
    int nCompSurfOld = altSurfNodes.nComp();
    int nCompSurf = nCompSurfOld + nCompAdd;
    altSurfNames.resize(nCompSurf);
    altSurfNames[nCompSurfOld + 0] = coldStrainName;

    int my_nNodes = 0;
    vector<Real> my_altSurfData, my_altSurfDatatmp;
    for (int lev = 0; lev<Nlev; ++lev)
    {
        for (MFIter mfi(*paths[lev]); mfi.isValid(); ++mfi)
        {
            const Box& box = mfi.validbox();
            const FArrayBox& pth = (*paths[lev])[mfi];
            
            if (box == ZBOX)
            {
                // Nothing to do
            }
            else
            {
                const int Npaths = box.length(0);
                my_nNodes += Npaths;

                const int loPath = box.smallEnd(1);
                const int hiPath = box.bigEnd(1);

                // For each path, search for where TComp equals TVal
                for (int i=0; i<Npaths; ++i)
                {
                    IntVect sIdx(D_DECL(i,0,0));
                    bool foundIt = false;
                    int lIdx=loPath, rIdx=lIdx+1;
                    Real lVal = pth(IntVect(D_DECL(i,lIdx,0)),TComp), rVal=lVal;
                    Real frac = 0.0;
                    
                    if (pth(IntVect(D_DECL(i,loPath,0)),TComp) > pth(IntVect(D_DECL(i,loPath,0)),TComp))
                        amrex::Abort("Path oriented backwards");
                    
                    if (TVal > pth(IntVect(D_DECL(i,hiPath,0)),TComp))
                    {
                        lIdx=hiPath-1;
                        rIdx=hiPath;
                        lVal = pth(IntVect(D_DECL(i,lIdx,0)),TComp); rVal=lVal;
                        frac = 1.0;
                    }
                    else if (TVal > lVal)
                    {
                        for (int j=loPath; (j<hiPath) && (!foundIt); ++j)
                        {
                            rIdx = lIdx + 1;  rVal = pth(IntVect(D_DECL(i,rIdx,0)),TComp);
                            
                            if ( (TVal >= lVal) && (TVal < rVal) )
                            {
                                foundIt = true;
                            }
                            else
                            {                                // Increment to next data pt
                                lIdx = rIdx; lVal = rVal;
                            }
                        }
                        frac = ( rVal==lVal ? 0.0 : (TVal-lVal)/(rVal-lVal) );
                    }

                    // Interpolate strainComp to this location
                    const Real& xL = pth(IntVect(D_DECL(i,lIdx,0)),strainComp);
                    const Real& xR = pth(IntVect(D_DECL(i,rIdx,0)),strainComp);
                    my_altSurfData.push_back(xL + (xR-xL) * frac);
                }
            }
        }
    }

    if (ParallelDescriptor::IOProcessor())
        std::cerr << "Cold Strain computed, communicating info to ioproc " << std::endl;

    // Communicate all data back to IOProc
    const int IOProc = ParallelDescriptor::IOProcessorNumber();        
    
    Vector<int> nNodes(ParallelDescriptor::NProcs(),0);
    Vector<int> offset(ParallelDescriptor::NProcs(),0);

#if BL_USE_MPI
    //
    // Tell root CPU how many tags each CPU will be sending.
    //
    MPI_Gather(&my_nNodes,
               1,
               ParallelDescriptor::Mpi_typemap<int>::type(),
               nNodes.dataPtr(),
               1,
               ParallelDescriptor::Mpi_typemap<int>::type(),
               IOProc,
               ParallelDescriptor::Communicator());
    
    if (ParallelDescriptor::IOProcessor())
    {
        for (int i = 0; i < nNodes.size(); i++)
            //
            // Convert from count of node data to add to count of reals to expect.
            //
            nNodes[i] *= nCompAdd;
        
        int N=0;
        for (int i = 0; i < nNodes.size(); i++)
            N += nNodes[i];

        my_altSurfDatatmp.resize(N);

        for (int i = 1; i < offset.size(); i++)
            offset[i] = offset[i-1] + nNodes[i-1];
    }
    //
    // Gather all the tags to IOProc into TheCollateSpace.
    //
    Real* bap = my_nNodes==0 ? 0 : &my_altSurfData[0];
    MPI_Gatherv(bap,
                my_nNodes*nCompAdd,
                ParallelDescriptor::Mpi_typemap<Real>::type(),
                (my_altSurfDatatmp.size()>0 ? &my_altSurfDatatmp[0] : 0),
                nNodes.dataPtr(),
                offset.dataPtr(),
                ParallelDescriptor::Mpi_typemap<Real>::type(),
                IOProc,
                ParallelDescriptor::Communicator());

    my_altSurfDatatmp.swap(my_altSurfData);

    if (ParallelDescriptor::IOProcessor())
        std::cerr << "Additional surface data sent" << std::endl;

#endif

    if (ParallelDescriptor::IOProcessor())
    {
        // Resize node fab to hold more data
        Box nBox = altSurfNodes.box();
        FArrayBox tmpFab(nBox,nCompSurfOld);
        tmpFab.copy(altSurfNodes,0,0,nCompSurfOld);
        altSurfNodes.resize(nBox,nCompSurf);
        altSurfNodes.copy(tmpFab,0,0,nCompSurfOld);

        Real** csd = new Real*[nCompAdd];
        for (int d=0; d<nCompAdd; ++d)
            csd[d] = altSurfNodes.dataPtr(nCompSurfOld + d);
        
        // De-scramble data so it looks like a isosurface node fab
        for (int lev=0; lev<Nlev; ++lev)
        {
            for (int i=0; i<paths[lev]->size(); ++i)
            {
                int proc = paths[lev]->DistributionMap()[i];
                int nPaths = inside_nodes[lev][i].size();
                for (int j=0; j<nPaths; ++j)
                {
                    int idx = inside_nodes[lev][i][j] - 1; // inside_nodes is 1-based

                    for (int d=0; d<nCompAdd; ++d)
                        csd[d][idx] = my_altSurfData[offset[proc] + d];

                    offset[proc]+= nCompAdd; // point to next box-ints on this proc
                }
            }
        }
        
    }
    if (ParallelDescriptor::IOProcessor())
        std::cerr << "New surface data added on ioproc " << std::endl;
}

void
add_thermal_thickness_to_surf(const Vector<MultiFab*>&  paths,
                              const vector<string>&     pathNames,
                              int                       xComp,
                              const Vector<Vector<Vector<int> > >& inside_nodes,
                              FArrayBox&                altSurfNodes,
                              vector<string>&           altSurfNames,
                              const string&             thickName,
                              int                       thickComp,
                              Real                      loVal,
                              Real                      hiVal)
{
    int Nlev = paths.size();
    Box ZBOX(IntVect::TheZeroVector(),IntVect::TheZeroVector());

    int nCompAdd = 1;
    int nCompSurfOld = altSurfNodes.nComp();
    int nCompSurf = nCompSurfOld + nCompAdd;
    altSurfNames.resize(nCompSurf);
    altSurfNames[nCompSurfOld + 0] = thickName;

    int my_nNodes = 0;
    vector<Real> my_altSurfData, my_altSurfDatatmp;
    for (int lev = 0; lev<Nlev; ++lev)
    {
        for (MFIter mfi(*paths[lev]); mfi.isValid(); ++mfi)
        {
            const Box& box = mfi.validbox();
            const FArrayBox& pth = (*paths[lev])[mfi];
            
            if (box == ZBOX)
            {
                // Nothing to do
            }
            else
            {
                const int Npaths = box.length(0);
                my_nNodes += Npaths;

                const int loPath = box.smallEnd(1);
                const int hiPath = box.bigEnd(1);

                // For each path, search for loVal, hiVal
                for (int i=0; i<Npaths; ++i)
                {
                    Real loLoc, hiLoc;
                    {
                        IntVect sIdx(D_DECL(i,0,0));
                        bool foundIt = false;
                        int lIdx=loPath, rIdx=lIdx+1;
                        Real lVal = pth(IntVect(D_DECL(i,lIdx,0)),thickComp), rVal=lVal;
                        Real frac = 0.0;
                        
                        if (pth(IntVect(D_DECL(i,loPath,0)),thickComp) > pth(IntVect(D_DECL(i,loPath,0)),thickComp))
                            amrex::Abort("Path oriented backwards");
                        
                        if (loVal > pth(IntVect(D_DECL(i,hiPath,0)),thickComp))
                        {
                            lIdx=hiPath-1;
                            rIdx=hiPath;
                            lVal = pth(IntVect(D_DECL(i,lIdx,0)),thickComp); rVal=lVal;
                            frac = 1.0;
                        }
                        else if (loVal > lVal)
                        {
                            for (int j=loPath; (j<hiPath) && (!foundIt); ++j)
                            {
                                rIdx = lIdx + 1;  rVal = pth(IntVect(D_DECL(i,rIdx,0)),thickComp);
                                
                                if ( (loVal >= lVal) && (loVal < rVal) )
                                {
                                    foundIt = true;
                                }
                                else
                                {
                                    // Increment to next data pt
                                    lIdx = rIdx; lVal = rVal;
                                }
                            }
                            frac = ( rVal==lVal ? 0.0 : (loVal-lVal)/(rVal-lVal) );
                        }
                        
                        // set distance to base pt of stream line
                        Real distance = 0;
                        if (lIdx<0)
                        {
                            for (int j=lIdx; j<0; ++j)
                            {
                                Real dd = 0;
                                for (int d=0; d<BL_SPACEDIM; ++d)
                                {
                                    const Real& xL = pth(IntVect(D_DECL(i,j  ,0)),xComp+d);
                                    const Real& xR = pth(IntVect(D_DECL(i,j+1,0)),xComp+d);
                                    dd += (xR-xL)*(xR-xL);
                                }
                                distance -= (j==lIdx ? (1.-frac) : 1.0) * std::sqrt(dd);
                            }
                        }
                        else
                        {
                            for (int j=0; j<rIdx; ++j)
                            {
                                Real dd = 0;
                                for (int d=0; d<BL_SPACEDIM; ++d)
                                {
                                    const Real& xL = pth(IntVect(D_DECL(i,j  ,0)),xComp+d);
                                    const Real& xR = pth(IntVect(D_DECL(i,j+1,0)),xComp+d);
                                    dd += (xR-xL)*(xR-xL);
                                }
                                distance += (j==rIdx-1 ? frac : 1.0) * std::sqrt(dd);
                            }
                        }
                        loLoc = distance;
                    }
                    {
                        IntVect sIdx(D_DECL(i,0,0));
                        bool foundIt = false;
                        int lIdx=loPath, rIdx=lIdx+1;
                        Real lVal = pth(IntVect(D_DECL(i,lIdx,0)),thickComp), rVal=lVal;
                        Real frac = 0.0;
                        
                        if (pth(IntVect(D_DECL(i,loPath,0)),thickComp) > pth(IntVect(D_DECL(i,loPath,0)),thickComp))
                            amrex::Abort("Path oriented backwards");
                        
                        if (hiVal > pth(IntVect(D_DECL(i,hiPath,0)),thickComp))
                        {
                            lIdx=hiPath-1;
                            rIdx=hiPath;
                            lVal = pth(IntVect(D_DECL(i,lIdx,0)),thickComp); rVal=lVal;
                            frac = 1.0;
                        }
                        else if (hiVal > lVal)
                        {
                            for (int j=loPath; (j<hiPath) && (!foundIt); ++j)
                            {
                                rIdx = lIdx + 1;  rVal = pth(IntVect(D_DECL(i,rIdx,0)),thickComp);
                                
                                if ( (hiVal >= lVal) && (hiVal < rVal) )
                                {
                                    foundIt = true;
                                }
                                else
                                {
                                    // Increment to next data pt
                                    lIdx = rIdx; lVal = rVal;
                                }
                            }
                            frac = ( rVal==lVal ? 0.0 : (hiVal-lVal)/(rVal-lVal) );
                        }
                        
                        // set distance to base pt of stream line
                        Real distance = 0;
                        if (lIdx<0)
                        {
                            for (int j=lIdx; j<0; ++j)
                            {
                                Real dd = 0;
                                for (int d=0; d<BL_SPACEDIM; ++d)
                                {
                                    const Real& xL = pth(IntVect(D_DECL(i,j  ,0)),xComp+d);
                                    const Real& xR = pth(IntVect(D_DECL(i,j+1,0)),xComp+d);
                                    dd += (xR-xL)*(xR-xL);
                                }
                                distance -= (j==lIdx ? (1.-frac) : 1.0) * std::sqrt(dd);
                            }
                        }
                        else
                        {
                            for (int j=0; j<rIdx; ++j)
                            {
                                Real dd = 0;
                                for (int d=0; d<BL_SPACEDIM; ++d)
                                {
                                    const Real& xL = pth(IntVect(D_DECL(i,j  ,0)),xComp+d);
                                    const Real& xR = pth(IntVect(D_DECL(i,j+1,0)),xComp+d);
                                    dd += (xR-xL)*(xR-xL);
                                }
                                distance += (j==rIdx-1 ? frac : 1.0) * std::sqrt(dd);
                            }
                        }

                        hiLoc = distance;
                    }
                    my_altSurfData.push_back(hiLoc - loLoc);
                }
            }
        }
    }

    if (ParallelDescriptor::IOProcessor())
        std::cerr << "Thermal thickness computed, communicating info to ioproc " << std::endl;

    // Communicate all data back to IOProc
    const int IOProc = ParallelDescriptor::IOProcessorNumber();        
    
    Vector<int> nNodes(ParallelDescriptor::NProcs(),0);
    Vector<int> offset(ParallelDescriptor::NProcs(),0);

#if BL_USE_MPI
    //
    // Tell root CPU how many tags each CPU will be sending.
    //
    MPI_Gather(&my_nNodes,
               1,
               ParallelDescriptor::Mpi_typemap<int>::type(),
               nNodes.dataPtr(),
               1,
               ParallelDescriptor::Mpi_typemap<int>::type(),
               IOProc,
               ParallelDescriptor::Communicator());
    
    if (ParallelDescriptor::IOProcessor())
    {
        for (int i = 0; i < nNodes.size(); i++)
            //
            // Convert from count of node data to add to count of reals to expect.
            //
            nNodes[i] *= nCompAdd;
        
        int N=0;
        for (int i = 0; i < nNodes.size(); i++)
            N += nNodes[i];

        my_altSurfDatatmp.resize(N);

        for (int i = 1; i < offset.size(); i++)
            offset[i] = offset[i-1] + nNodes[i-1];
    }
    //
    // Gather all the tags to IOProc into TheCollateSpace.
    //
    Real* bap = my_nNodes==0 ? 0 : &my_altSurfData[0];
    MPI_Gatherv(bap,
                my_nNodes*nCompAdd,
                ParallelDescriptor::Mpi_typemap<Real>::type(),
                (my_altSurfDatatmp.size()>0 ? &my_altSurfDatatmp[0] : 0),
                nNodes.dataPtr(),
                offset.dataPtr(),
                ParallelDescriptor::Mpi_typemap<Real>::type(),
                IOProc,
                ParallelDescriptor::Communicator());

    my_altSurfDatatmp.swap(my_altSurfData);

    if (ParallelDescriptor::IOProcessor())
        std::cerr << "Additional surface data sent" << std::endl;

#endif

    if (ParallelDescriptor::IOProcessor())
    {
        // Resize node fab to hold more data
        Box nBox = altSurfNodes.box();
        FArrayBox tmpFab(nBox,nCompSurfOld);
        tmpFab.copy(altSurfNodes,0,0,nCompSurfOld);
        altSurfNodes.resize(nBox,nCompSurf);
        altSurfNodes.copy(tmpFab,0,0,nCompSurfOld);

        Real** csd = new Real*[nCompAdd];
        for (int d=0; d<nCompAdd; ++d)
            csd[d] = altSurfNodes.dataPtr(nCompSurfOld + d);
        
        // De-scramble data so it looks like a isosurface node fab
        for (int lev=0; lev<Nlev; ++lev)
        {
            for (int i=0; i<paths[lev]->size(); ++i)
            {
                int proc = paths[lev]->DistributionMap()[i];
                int nPaths = inside_nodes[lev][i].size();
                for (int j=0; j<nPaths; ++j)
                {
                    int idx = inside_nodes[lev][i][j] - 1; // inside_nodes is 1-based

                    for (int d=0; d<nCompAdd; ++d)
                        csd[d][idx] = my_altSurfData[offset[proc] + d];

                    offset[proc]+= nCompAdd; // point to next box-ints on this proc
                }
            }
        }
        
    }
    if (ParallelDescriptor::IOProcessor())
        std::cerr << "New surface data added on ioproc " << std::endl;
}

void
build_surface_at_isoVal(const Vector<MultiFab*>&  paths,
                        const vector<string>&     pathNames,
                        int                       xComp,
                        int                       isoComp,
                        Real                      isoVal,
                        const Vector<Vector<Vector<int> > >& inside_nodes,
                        const Vector<int>&         comps,
                        FArrayBox&                altSurfNodes,
                        vector<string>&           altSurfNames,
                        string&                   distanceVarName)
{
    int Nlev = paths.size();
    Box ZBOX(IntVect::TheZeroVector(),IntVect::TheZeroVector());

    int nComp = comps.size();

    // build a structure to hold the locations of isoVal and other scalars, in the format of streamlines
    int nCompSurf = BL_SPACEDIM+nComp+1;
    altSurfNames.resize(nCompSurf);
    for (int d=0; d<BL_SPACEDIM; ++d)
        altSurfNames[d] = pathNames[xComp+d];
    for (int d=0; d<comps.size(); ++d)
        altSurfNames[BL_SPACEDIM+d] = pathNames[comps[d]];
    altSurfNames[BL_SPACEDIM+nComp] = distanceVarName;

    int my_nNodes = 0;
    vector<Real> my_altSurfData;
    vector<Real> my_altSurfDatatmp;
    for (int lev = 0; lev<Nlev; ++lev)
    {
        if (ParallelDescriptor::IOProcessor())
            std::cerr << "building new surface at level: " << lev << std::endl;

        for (MFIter mfi(*paths[lev]); mfi.isValid(); ++mfi)
        {
            const Box& box = mfi.validbox();
            const FArrayBox& pth = (*paths[lev])[mfi];
            
            if (box == ZBOX)
            {
                // Nothing to do
            }
            else
            {
                const int Npaths = box.length(0);
                my_nNodes += Npaths;

                const int loPath = box.smallEnd(1);
                const int hiPath = box.bigEnd(1);

                // For each path, search for isoVal
                for (int i=0; i<Npaths; ++i)
                {
                    IntVect sIdx(D_DECL(i,0,0));
                    bool foundIt = false;
                    int lIdx=loPath, rIdx=lIdx+1;
                    Real lVal = pth(IntVect(D_DECL(i,lIdx,0)),isoComp), rVal=lVal;
                    Real frac = 0.0;

                    if (pth(IntVect(D_DECL(i,loPath,0)),isoComp) > pth(IntVect(D_DECL(i,loPath,0)),isoComp))
                        amrex::Abort("Path oriented backwards");

                    if (isoVal > pth(IntVect(D_DECL(i,hiPath,0)),isoComp))
                    {
                        lIdx=hiPath-1;
                        rIdx=hiPath;
                        lVal = pth(IntVect(D_DECL(i,lIdx,0)),isoComp); rVal=lVal;
                        frac = 1.0;
                    }
                    else if (isoVal > lVal)
                    {
                        for (int j=loPath; (j<hiPath) && (!foundIt); ++j)
                        {
                            rIdx = lIdx + 1;  rVal = pth(IntVect(D_DECL(i,rIdx,0)),isoComp);

                            if ( (isoVal >= lVal) && (isoVal < rVal) )
                            {
                                foundIt = true;
                            }
                            else
                            {
                                // Increment to next data pt
                                lIdx = rIdx; lVal = rVal;
                            }
                        }
                        frac = ( rVal==lVal ? 0.0 : (isoVal-lVal)/(rVal-lVal) );
                    }

                    // Compute position at isoVal
                    for (int d=0; d<BL_SPACEDIM; ++d)
                    {
                        const Real& xL = pth(IntVect(D_DECL(i,lIdx,0)),xComp+d);
                        const Real& xR = pth(IntVect(D_DECL(i,rIdx,0)),xComp+d);
                        my_altSurfData.push_back(xL + (xR-xL) * frac);
                    }
                    // Compute other values at isoVal as well
                    for (int d=0; d<nComp; ++d)
                    {
                        const Real& xL = pth(IntVect(D_DECL(i,lIdx,0)),comps[d]);
                        const Real& xR = pth(IntVect(D_DECL(i,rIdx,0)),comps[d]);
                        my_altSurfData.push_back(xL + (xR-xL) * frac);
                    }
                    // set distance to base pt of stream line (isoComp == isoVal)
                    Real distance = 0;
                    if (lIdx<0)
                    {
                        for (int j=lIdx; j<0; ++j)
                        {
                            Real dd = 0;
                            for (int d=0; d<BL_SPACEDIM; ++d)
                            {
                                const Real& xL = pth(IntVect(D_DECL(i,j  ,0)),xComp+d);
                                const Real& xR = pth(IntVect(D_DECL(i,j+1,0)),xComp+d);
                                dd += (xR-xL)*(xR-xL);
                            }
                            distance -= (j==lIdx ? (1.-frac) : 1.0) * std::sqrt(dd);
                        }
                    }
                    else
                    {
                        for (int j=0; j<rIdx; ++j)
                        {
                            Real dd = 0;
                            for (int d=0; d<BL_SPACEDIM; ++d)
                            {
                                const Real& xL = pth(IntVect(D_DECL(i,j  ,0)),xComp+d);
                                const Real& xR = pth(IntVect(D_DECL(i,j+1,0)),xComp+d);
                                dd += (xR-xL)*(xR-xL);
                            }
                            distance += (j==rIdx-1 ? frac : 1.0) * std::sqrt(dd);
                        }
                    }
                    my_altSurfData.push_back(distance);
                }
            }
        }
    }

    if (ParallelDescriptor::IOProcessor())
        std::cerr << "New surface built, communicating info to ioproc " << std::endl;

    // Communicate all data back to IOProc
    const int IOProc = ParallelDescriptor::IOProcessorNumber();        
    
    Vector<int> nNodes(ParallelDescriptor::NProcs(),0);
    Vector<int> offset(ParallelDescriptor::NProcs(),0);

#if BL_USE_MPI
    //
    // Tell root CPU how many tags each CPU will be sending.
    //
    MPI_Gather(&my_nNodes,
               1,
               ParallelDescriptor::Mpi_typemap<int>::type(),
               nNodes.dataPtr(),
               1,
               ParallelDescriptor::Mpi_typemap<int>::type(),
               IOProc,
               ParallelDescriptor::Communicator());
    
    if (ParallelDescriptor::IOProcessor())
    {
        for (int i = 0; i < nNodes.size(); i++)
            //
            // Convert from count of nodes to count of reals to expect.
            //
            nNodes[i] *= nCompSurf;
        
        int N=0;
        for (int i = 0; i < nNodes.size(); i++)
            N += nNodes[i];

        my_altSurfDatatmp.resize(N);

        for (int i = 1; i < offset.size(); i++)
            offset[i] = offset[i-1] + nNodes[i-1];
    }
    //
    // Gather all the tags to IOProc into TheCollateSpace.
    //
    Real* bap = my_nNodes==0 ? 0 : &my_altSurfData[0];
    MPI_Gatherv(bap,
                my_nNodes*nCompSurf,
                ParallelDescriptor::Mpi_typemap<Real>::type(),
                (my_altSurfDatatmp.size()>0 ? &my_altSurfDatatmp[0] : 0),
                nNodes.dataPtr(),
                offset.dataPtr(),
                ParallelDescriptor::Mpi_typemap<Real>::type(),
                IOProc,
                ParallelDescriptor::Communicator());

    my_altSurfDatatmp.swap(my_altSurfData);

#endif

    if (ParallelDescriptor::IOProcessor())
        std::cerr << "Data sent, building new iso" << std::endl;

    if (ParallelDescriptor::IOProcessor())
    {
        // Count total number of nodes
        int num_nodes = 0;
        for (int i=0; i<Nlev; ++i)
            for (int j=0; j<inside_nodes[i].size(); ++j)
                num_nodes += inside_nodes[i][j].size();
        
        Box nBox(IntVect::TheZeroVector(),(num_nodes-1)*amrex::BASISV(0));
        altSurfNodes.resize(nBox,nCompSurf);
        Real** csd = new Real*[nCompSurf];
        for (int d=0; d<nCompSurf; ++d)
            csd[d] = altSurfNodes.dataPtr(d);
        
        // De-scramble data so it looks like a isosurface node fab
        for (int lev=0; lev<Nlev; ++lev)
        {
            for (int i=0; i<paths[lev]->size(); ++i)
            {
                int proc = paths[lev]->DistributionMap()[i];
                int nPaths = inside_nodes[lev][i].size();
                for (int j=0; j<nPaths; ++j)
                {
                    int idx = inside_nodes[lev][i][j] - 1;

                    for (int d=0; d<nCompSurf; ++d)
                        csd[d][idx] = my_altSurfData[offset[proc] + d];

                    offset[proc]+= nCompSurf; // point to next box-ints on this proc
                }
            }
        }
    }
    if (ParallelDescriptor::IOProcessor())
        std::cerr << "New iso built on ioproc " << std::endl;
}

void
advect_alt_surf(const Vector<MultiFab*>&  paths,
                const vector<string>&     pathNames,
                int                       xComp,
                int                       vComp,
                Real                      dt,
                int                       isoComp,
                Real                      isoVal,
                const Vector<Vector<Vector<int> > >& inside_nodes,
                FArrayBox&                altSurfNodes,
                vector<string>&           altSurfNames)
{
}

void
write_ml_streamline_data(const std::string&       outfile,
                         const Vector<MultiFab*>& data,
                         int                      sComp,
                         const vector<string>&    names,
                         const Vector<int>&       faceData,
                         int                      nElts,
                         const Vector<Vector<Vector<int> > >& inside_nodes,
                         const AmrData&           amrdata)
{
    const std::string LevelDirName("Level");
    const std::string StreamDataFileName("Str");

#undef OLDFORMAT
#define OLDFORMAT 1

#if defined(OLDFORMAT)
    const std::string FileFormatName("Oddball-multilevel-connected-data-format");
#else
    const std::string FileFormatName("Oddball-multilevel-connected-data-format-1.0");
#endif

    int Nlevels = data.size();
    Vector<Box> probDomain(Nlevels);
    Vector<BoxArray> stateBoxArray(Nlevels);
    for (int i=0; i<Nlevels; ++i) {
        probDomain[i] = amrdata.ProbDomain()[i];
        stateBoxArray[i] = amrdata.boxArray(i);
    }

    if (ParallelDescriptor::IOProcessor())
    {
        if (!amrex::UtilCreateDirectory(outfile, 0755))
            amrex::CreateDirectoryFailed(outfile);

        const std::string HeaderFileName("Header");
        const std::string ElementFileName("Elements");
        const std::string ElementFileFormat("ELEMENT_DATA_ASCII");

        // Write Header
        const string FullHeaderFileName = outfile + "/" + HeaderFileName;
        std::ofstream ofh;
        ofh.open(FullHeaderFileName.c_str());
        ofh << FileFormatName << '\n';
        ofh << data.size() << '\n'; // number of levels
        ofh << names.size() << '\n'; // number of variables per node
        for (int i=0; i<names.size(); ++i) {
            ofh << names[i] << '\n'; 
        }

#if ! defined(OLDFORMAT)
        ofh << ElementFileName << '\n'; 
        ofh << ElementFileFormat << '\n'; 

        for (int i=0; i<BL_SPACEDIM; ++i)  ofh << probLo[i] << " ";
        ofh << '\n';
        for (int i=0; i<BL_SPACEDIM; ++i)  ofh << probHi[i] << " ";
        ofh << '\n';
        for (int lev=0; lev<data.size(); ++lev) {
            ofh << probDomain[lev] << '\n';
            ofh << stateBoxArray[lev] << '\n';
        }
#endif
        ofh.close();

        // Write elements
        const string FullElementFileName = outfile + "/" + ElementFileName;
        int nodesPerElt = faceData.size() / nElts;
        std::ofstream ofe;
        ofe.open(FullElementFileName.c_str());
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
            if (!amrex::UtilCreateDirectory(DataFile, 0755))
                amrex::CreateDirectoryFailed(DataFile);
        //
        // Force other processors to wait till directory is built.
        //
        ParallelDescriptor::Barrier();
        
        std::string FabDataFile = DataFile + "/Str";

        if (names.size() == data[lev]->size())
        {
            VisMF::Write(*data[lev],FabDataFile);
        }
        else
        {
            MultiFab tmp(data[lev]->boxArray(),data[lev]->DistributionMap(),names.size(),data[lev]->nGrow());
            MultiFab::Copy(tmp,*data[lev],sComp,0,names.size(),data[lev]->nGrow());
            VisMF::Write(tmp,FabDataFile);
        }
    }
}

void
dump_ml_streamline_data(const std::string&       outFile,
                        const Vector<MultiFab*>& data,
                        int                      sComp,
                        const vector<string>&    names,
                        const Vector<int>&       faceData,
                        int                      nElts,
                        const Vector<Vector<Vector<int> > >& inside_nodes)
{
  // Create a folder and have each processor write their own data, one file per streamline
  auto myProc = ParallelDescriptor::MyProc();
  auto nProcs = ParallelDescriptor::NProcs();
  int cnt = 0;
  
  if (!amrex::UtilCreateDirectory(outFile, 0755))
    amrex::CreateDirectoryFailed(outFile);
  ParallelDescriptor::Barrier();
  
  const Box nullBox(IntVect::TheZeroVector(),IntVect::TheZeroVector());
  for (int lev=0; lev<data.size(); ++lev)
  {
    const auto& ba = data[lev]->boxArray();
    const auto& dm = data[lev]->DistributionMap();
    int nComp = data[lev]->nComp();
    for (int j=0; j<ba.size(); ++j)
    {
      if (dm[j] == myProc)
      {
        const auto& b = ba[j];
        if (b!=nullBox)
        {
          std::string fileName = outFile + "/str_";
          fileName = Concatenate(fileName,myProc) + "_";
          fileName = Concatenate(fileName,cnt++);
          std::ofstream ofs(fileName.c_str());

          for (int i=0; i<names.size(); ++i) {
            ofs << names[i] << " ";
          }
          ofs << '\n';
          
          for (auto iv=b.smallEnd(); iv<=b.bigEnd(); b.next(iv))
          {
            for (int i=sComp; i<nComp; ++i)
            {
              ofs << (*data[lev])[j](iv,i) << " ";
            }
            ofs << '\n';
          }
          ofs.close();
        }
      }
      ParallelDescriptor::Barrier();
    }
  }
}
