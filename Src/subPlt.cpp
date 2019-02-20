#include "AMReX_ParmParse.H"
#include "AMReX_ParallelDescriptor.H"
#include "AMReX_DataServices.H"
#include "AMReX_Utility.H"
#include "WritePlotFile.H"

using std::cout;
using std::endl;
using std::cerr;
using namespace amrex;

const bool verbose_DEF   = false;

static
void
PrintUsage (const char* progName)
{
    cout << '\n';
    cout << "Usage:" << '\n';
    cout << progName << '\n';
    cout << "    infile=inFileName" << '\n';
    cout << "   [outfile=outFileName <defaults to <inFileName>_section>]" << '\n';
    cout << "   [-sComp=N <defaults to 0, unless \"comps\" used]" << '\n';
    cout << "   [-nComp=N <defaults to all in <inFileName>, unless comps used]" << '\n';
    cout << "   [-comps=\"N1 N2 N3...\"]" << '\n';
    cout << "   [-help]" << '\n';
    cout << "   [-verbose]" << '\n';
    cout << '\n';
    exit(1);
}

static Vector< Vector<int> > contigLists(const Vector<int> orig);

int
main (int   argc,
      char* argv[])
{
    if (argc == 1)
        PrintUsage(argv[0]);

    Initialize(argc,argv);    

    ParmParse pp;

    if (pp.contains("help"))
        PrintUsage(argv[0]);

    FArrayBox::setFormat(FABio::FAB_IEEE_32);
    //
    // Scan the arguments.
    //
    std::string infile;

    bool verbose = verbose_DEF;
    verbose = (pp.contains("verbose") ? true : false);
    if (verbose)
        AmrData::SetVerbose(true);

    pp.get("infile",infile);

    vector<std::string> pieces = Tokenize(infile,std::string("/"));
    std::string outfile = pieces[pieces.size()-1] + std::string("_section");
    pp.query("outfile",outfile);

    // Read in pltfile and get amrData ref
    DataServices::SetBatchMode();
    Amrvis::FileType fileType(Amrvis::NEWPLT);
    DataServices dataServices(infile, fileType);
    if (!dataServices.AmrDataOk())
        DataServices::Dispatch(DataServices::ExitRequest, NULL);
    AmrData& amrData = dataServices.AmrDataRef();

    Vector<int> comps;
    if (int nc = pp.countval("comps"))
    {
        comps.resize(nc);
        pp.getarr("comps",comps,0,nc);
    }
    else
    {
        int sComp = 0;
        pp.query("sComp",sComp);
        int nComp = amrData.NComp();
        pp.query("nComp",nComp);
        BL_ASSERT(sComp+nComp <= amrData.NComp());
        comps.resize(nComp);
        for (int i=0; i<nComp; ++i)
            comps[i] = sComp + i;
    }

    Vector< Vector<int> > compsNEW = contigLists(comps);
    int finestLevel = amrData.FinestLevel(); pp.query("finestLevel",finestLevel);
    Box subbox = amrData.ProbDomain()[finestLevel];
    Vector<int> inBox;

    if (int nx=pp.countval("box"))
    {
        pp.getarr("box",inBox,0,nx);
        int d=BL_SPACEDIM;
        BL_ASSERT(inBox.size()==2*d);
        subbox=Box(IntVect(D_DECL(inBox[0],inBox[1],inBox[2])),
                   IntVect(D_DECL(inBox[d],inBox[d+1],inBox[d+2])),
                   IndexType::TheCellType());
    }

    Vector<Box> subboxes(finestLevel+1,subbox);
    for (int iLevel = finestLevel-1; iLevel>=0; --iLevel)
        subboxes[iLevel] = coarsen(subboxes[iLevel+1],amrData.RefRatio()[iLevel]);
    for (int iLevel = 1; iLevel<=finestLevel; ++iLevel)
        subboxes[iLevel] = refine(subboxes[iLevel-1],amrData.RefRatio()[iLevel-1]);

    Vector<MultiFab*> data_sub(finestLevel+1);
    Vector<std::string> names(comps.size());
    Vector<std::string> subNames;
    Vector<int> fillComps;
   
    Vector<Real> plo(BL_SPACEDIM), phi(BL_SPACEDIM);
    Vector<Box> psize(finestLevel+1);
    const IntVect ilo = subbox.smallEnd();
    const IntVect ihi = subbox.bigEnd();

   for (int i =0 ; i< BL_SPACEDIM; i++) {
       
       plo[i] = amrData.ProbLo()[i]+(ilo[i])*amrData.DxLevel()[finestLevel][i];
       phi[i] = amrData.ProbLo()[i]+(ihi[i]+1)*amrData.DxLevel()[finestLevel][i];
       
     }
        
    int actual_lev = -1;
    for (int iLevel = 0; iLevel <= finestLevel; ++iLevel)
    {
        // Build the BoxArray for the result
        BoxArray ba_sub = intersect(amrData.boxArray(iLevel),subboxes[iLevel]);

        if (ba_sub.size() > 0)
        {
            actual_lev = iLevel;
            DistributionMapping dmap_sub(ba_sub);
            data_sub[iLevel] = new MultiFab(ba_sub,dmap_sub,comps.size(),0);

            for (int i=0; i<comps.size(); ++i)
            {
                data_sub[iLevel]->copy(amrData.GetGrids(iLevel,comps[i],subboxes[iLevel]),0,i,1);
                amrData.FlushGrids(comps[i]);                
                names[i] = amrData.PlotVarNames()[comps[i]];
                
                if (ParallelDescriptor::IOProcessor()) 
                    cout << "Filling " << names[i] << " on level " << iLevel << endl;
                
            }
        }
    }

    // Write out the subregion pltfile
    WritePlotFile(data_sub,subboxes,amrData,outfile,verbose,names);

    Finalize();
    return 0;
}

static Vector< Vector<int> > contigLists(const Vector<int> orig)
{
    BL_ASSERT(orig.size() > 0);
    Vector< Vector<int> > res(1);
    int mySet = 0;
    res[mySet].resize(1);
    int myCnt = 0;
    res[mySet][myCnt] = orig[0];

    for (int i=1; i<orig.size(); ++i)
    {
        if (res[mySet][myCnt]+1 == orig[i])
        {
            res[mySet].resize(++myCnt + 1);
            res[mySet][myCnt] = orig[i];
        }
        else
        {
            myCnt = 0;
            res.resize(++mySet+1);
            res[mySet].resize(1);
            res[mySet][myCnt] = orig[i];
        }
    }
    return res;
}


