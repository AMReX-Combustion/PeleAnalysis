#include <string>
#include <iostream>

#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include <AMReX_DataServices.H>
#include <WritePlotFile.H>
#include <AppendToPlotFile.H>

using namespace amrex;

#define Y_IN_PLOTFILE
#undef Y_IN_PLOTFILE

static
void 
print_usage (int,
             char* argv[])
{
	std::cerr << "usage:\n";
    std::cerr << argv[0] << " infile=<name> outfile=<> [options] \n\tOptions:\n";
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

    std::string infile;
    pp.get("infile",infile);

    std::string outfile;
    pp.get("outfile",outfile);

    DataServices::SetBatchMode();
    Amrvis::FileType fileType(Amrvis::NEWPLT);

    // Set up for reading pltfile
    DataServices dataServices(infile, fileType);

    if (!dataServices.AmrDataOk())
        //
        // This calls ParallelDescriptor::EndParallel() and exit()
        //
        DataServices::Dispatch(DataServices::ExitRequest, NULL);
    
    AmrData& amrData = dataServices.AmrDataRef();

    int nComp = 0;
    Vector<int> comps;
    if (int nc = pp.countval("comps"))
    {
        comps.resize(nc);
        pp.getarr("comps",comps,0,nc);
        nComp = comps.size();
    }
    else
    {
        int sComp = 0;
        pp.query("sComp",sComp);
        nComp = amrData.NComp();
        pp.query("nComp",nComp);
        BL_ASSERT(sComp+nComp <= amrData.NComp());
        comps.resize(nComp);
        for (int i=0; i<nComp; ++i)
            comps[i] = sComp + i;
    }

    int Nlev = amrData.FinestLevel() + 1;

    int max_grid_size = 128;
    pp.query("max_grid_size",max_grid_size);
    Vector<MultiFab*> fileData(Nlev);
    for (int lev=0; lev<Nlev; ++lev)
    {
        BoxArray newba = amrData.boxArray(lev);
        newba.maxSize(max_grid_size);
        const DistributionMapping dm(newba);
        fileData[lev] = new MultiFab(newba,dm,nComp,0);
    }
    if (ParallelDescriptor::IOProcessor())
        std::cerr << "Full MultiFab allocated " << std::endl;

    for (int lev=0; lev<Nlev; ++lev)
    {
        for (int i=0; i<comps.size(); ++i)
        {
            fileData[lev]->copy(amrData.GetGrids(lev,comps[i]),0,i,1);
            if (ParallelDescriptor::IOProcessor())
                std::cerr << "After GetGrids: " << amrData.PlotVarNames()[comps[i]] << std::endl;
            amrData.FlushGrids(comps[i]);
            if (ParallelDescriptor::IOProcessor())
                std::cerr << "AmrData flushed: " << amrData.PlotVarNames()[comps[i]] << std::endl;
        }
    }
    if (ParallelDescriptor::IOProcessor())
        std::cerr << "File data loaded" << std::endl;

    Vector<std::string> names(nComp);
    for (int i=0; i<comps.size(); ++i)
        names[i] = amrData.PlotVarNames()[comps[i]];
    
    bool verb = false;
    WritePlotFile(fileData,amrData,outfile,verb,names);
    }
    amrex::Finalize();
    return 0;
}
