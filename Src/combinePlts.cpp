#include <string>
#include <iostream>

#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include <AMReX_DataServices.H>
#include <AMReX_WritePlotFile.H>

using namespace amrex;

#define Y_IN_PLOTFILE
#undef Y_IN_PLOTFILE

static
void 
print_usage (int,
             char* argv[])
{
	std::cerr << "usage:\n";
    std::cerr << argv[0] << " infileL=<name> infileR=<name> outfile=<> [options] \n\tOptions:\n";
    exit(1);
}

std::string
getFileRoot(const std::string& infile)
{
    vector<std::string> tokens = Tokenize(infile,std::string("/"));
    return tokens[tokens.size()-1];
}

int
main (int   argc,
      char* argv[])
{
    amrex::Initialize(argc,argv);
    {
    if (argc < 2)
        print_usage(argc,argv);

    ParmParse pp;

    if (pp.contains("help"))
        print_usage(argc,argv);

    if (pp.contains("verbose"))
        AmrData::SetVerbose(true);

    std::string infileL;
    pp.get("infileL",infileL);

    std::string infileR;
    pp.get("infileR",infileR);

    std::string outfile;
    pp.get("outfile",outfile);

    DataServices::SetBatchMode();
    Amrvis::FileType fileType(Amrvis::NEWPLT);

    // Set up for reading left pltfile
    DataServices dataServicesL(infileL, fileType);

    if (!dataServicesL.AmrDataOk())
        //
        // This calls ParallelDescriptor::EndParallel() and exit()
        //
        DataServices::Dispatch(DataServices::ExitRequest, NULL);
    
    AmrData& amrDataL = dataServicesL.AmrDataRef();

    Vector<int> compsL;
    if (int nc = pp.countval("compsL"))
    {
        compsL.resize(nc);
        pp.getarr("compsL",compsL,0,nc);
    }
    else
    {
        int sCompL = 0;
        pp.query("sCompL",sCompL);
        int nCompL = amrDataL.NComp();
        pp.query("nCompL",nCompL);
        BL_ASSERT(sCompL+nCompL <= amrDataL.NComp());
        compsL.resize(nCompL);
        for (int i=0; i<nCompL; ++i)
            compsL[i] = sCompL + i;
    }

    // Set up for reading right pltfile
    DataServices dataServicesR(infileR, fileType);

    if (!dataServicesR.AmrDataOk())
        //
        // This calls ParallelDescriptor::EndParallel() and exit()
        //
        DataServices::Dispatch(DataServices::ExitRequest, NULL);
    
    AmrData& amrDataR = dataServicesR.AmrDataRef();

    Vector<int> compsR;
    if (int nc = pp.countval("compsR"))
    {
        compsR.resize(nc);
        pp.getarr("compsR",compsR,0,nc);
    }
    else
    {
        int sCompR = 0;
        pp.query("sCompR",sCompR);
        int nCompR = amrDataR.NComp();
        pp.query("nCompR",nCompR);
        BL_ASSERT(sCompR+nCompR <= amrDataR.NComp());
        compsR.resize(nCompR);
        for (int i=0; i<nCompR; ++i)
            compsR[i] = sCompR + i;
    }

    int Nlev = amrDataL.FinestLevel() + 1;
    BL_ASSERT(Nlev == amrDataR.FinestLevel() + 1);

    const int nComp = compsL.size() + compsR.size();
    Vector<MultiFab*> fileData(Nlev);
    for (int lev=0; lev<Nlev; ++lev)
    {
        BL_ASSERT(amrDataL.boxArray(lev) == amrDataR.boxArray(lev));
        const DistributionMapping dm(amrDataL.boxArray(lev));
        fileData[lev] = new MultiFab(amrDataL.boxArray(lev),dm,nComp,0);
    }
    if (ParallelDescriptor::IOProcessor())
        std::cerr << "Full MultiFab allocated " << std::endl;

    for (int lev=0; lev<Nlev; ++lev)
    {
        for (int i=0; i<compsL.size(); ++i)
        {
            fileData[lev]->copy(amrDataL.GetGrids(lev,compsL[i]),0,i,1);
            if (ParallelDescriptor::IOProcessor())
                std::cerr << "After GetGrids (L): " << amrDataL.PlotVarNames()[compsL[i]] << std::endl;
            //amrDataL.FlushGrids(compsL[i]);
            if (ParallelDescriptor::IOProcessor())
                std::cerr << "AmrData flushed (L): " << amrDataL.PlotVarNames()[compsL[i]] << std::endl;
        }
        for (int i=0; i<compsR.size(); ++i)
        {
            fileData[lev]->copy(amrDataR.GetGrids(lev,compsR[i]),0,compsL.size()+i,1);
            if (ParallelDescriptor::IOProcessor())
                std::cerr << "After GetGrids (R): " << amrDataR.PlotVarNames()[compsR[i]] << std::endl;
            //amrDataR.FlushGrids(compsR[i]);
            if (ParallelDescriptor::IOProcessor())
                std::cerr << "AmrData flushed (R): " << amrDataR.PlotVarNames()[compsR[i]] << std::endl;
        }
    }
    if (ParallelDescriptor::IOProcessor())
        std::cerr << "File data loaded" << std::endl;

    Real progMin, progMax;
    for (int i=0; i<compsL.size(); ++i) {
        amrDataL.MinMax(amrDataL.ProbDomain()[Nlev-1], amrDataL.PlotVarNames()[compsL[i]], Nlev-1, progMin, progMax);
        if (ParallelDescriptor::IOProcessor()) {
            std::cout << amrDataL.PlotVarNames()[compsL[i]] << " min/max: " << progMin << ", " << progMax << std::endl;
        }
    }
    for (int i=0; i<compsR.size(); ++i) {
        amrDataR.MinMax(amrDataR.ProbDomain()[Nlev-1], amrDataR.PlotVarNames()[compsR[i]], Nlev-1, progMin, progMax);
        if (ParallelDescriptor::IOProcessor()) {
            std::cout << amrDataR.PlotVarNames()[compsR[i]] << " min/max: " << progMin << ", " << progMax << std::endl;
        }
    }

    Vector<std::string> names(nComp);
    for (int i=0; i<compsL.size(); ++i)
        names[i] = amrDataL.PlotVarNames()[compsL[i]];
    for (int i=0; i<compsR.size(); ++i)
        names[compsL.size()+i] = amrDataR.PlotVarNames()[compsR[i]];
    
    bool verb = false;
    WritePlotFile(fileData,amrDataL,outfile,verb,names);
    }
    amrex::Finalize();
    return 0;
}
