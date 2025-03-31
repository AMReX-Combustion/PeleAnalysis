#include <string>
#include <iostream>
#include <set>

#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include <AMReX_DataServices.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_PlotFileUtilHDF5.H>

using namespace amrex;

static
void 
print_usage (int,
             char* argv[])
{
  std::cerr << "usage:\n";
  std::cerr << argv[0] << " infile infile=f1 [options] \n\tOptions:\n";
  exit(1);
}

std::string
getFileRoot(const std::string& infile)
{
  std::vector<std::string> tokens = Tokenize(infile,std::string("/"));
  return tokens[tokens.size()-1];
}

int
main (int   argc,
      char* argv[])
{
  amrex::Initialize(argc,argv);
  {
    if (argc < 2) {
      print_usage(argc,argv);
    }

    // ---------------------------------------------------------------------
    // Set defaults input values
    // ---------------------------------------------------------------------
    int finestLevel           = 1000;

    // ---------------------------------------------------------------------
    // ParmParse
    // ---------------------------------------------------------------------
    ParmParse pp;

    if (pp.contains("help")) {
      print_usage(argc,argv);
    }

    pp.query("finestLevel",finestLevel);

    std::string plotFileName; pp.get("infile",plotFileName);
    
   
    // Initialize DataService
    DataServices::SetBatchMode();
    Amrvis::FileType fileType(Amrvis::NEWPLT);
    DataServices dataServices(plotFileName, fileType);
    if( ! dataServices.AmrDataOk()) {
      DataServices::Dispatch(DataServices::ExitRequest, NULL);
    }
    AmrData& amrData = dataServices.AmrDataRef();

    // Plotfile global infos
    finestLevel = std::min(finestLevel,amrData.FinestLevel());
    int Nlev = finestLevel + 1;
    const Vector<std::string>& plotVarNames = amrData.PlotVarNames();
    int nvars = plotVarNames.size();
    int id_comp_last = 0;    
    RealBox rb(&(amrData.ProbLo()[0]), 
               &(amrData.ProbHi()[0]));


    // Check symmetry/periodicity in given coordinate direction
    Vector<int> sym_dir(AMREX_SPACEDIM,0);
    pp.queryarr("sym_dir",sym_dir,0,AMREX_SPACEDIM);  

    Vector<int> is_per(AMREX_SPACEDIM,0);
    pp.queryarr("is_per",is_per,0,AMREX_SPACEDIM);
    Print() << "Periodicity assumed for this case: ";
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        Print() << is_per[idim] << " ";
    }
    Print() << "\n";

    int coord = 0;

    // ---------------------------------------------------------------------
    // Let's start the real work
    // ---------------------------------------------------------------------
    Vector<MultiFab*> fileData(Nlev);
    Vector<Geometry> geoms(Nlev);
    const int nGrow = 1;

    // Read data on all the levels
    for (int lev=0; lev<Nlev; ++lev) {
      const DistributionMapping dm(amrData.boxArray(lev));
      geoms[lev] = Geometry(amrData.ProbDomain()[lev],&rb,coord,&(is_per[0]));
      fileData[lev] = new MultiFab(amrData.boxArray(lev),dm,nvars,0);
    }
 
    Vector<int> idcomp;
    for (int i = 0; i < plotVarNames.size(); ++i) {  // loop though current pltfile variable names
	  idcomp.push_back(i);                       // sets index location of plotVarName that matches
      }
 
        
    for (int lev=0; lev<Nlev; ++lev) {
	for (int i=0; i<nvars; ++i) {
	  fileData[lev]->ParallelCopy(amrData.GetGrids(lev,idcomp[i]),0,id_comp_last+i,1);	    
	}
    }   
  

    // ---------------------------------------------------------------------
    // Write the results
    // ---------------------------------------------------------------------    
    
    std::string hdf5_compression{"ZLIB@9"};
    std::string outfile(getFileRoot(plotFileName) + "_hdf5");
    Print() << "Writing new data to " << outfile << std::endl;
    Vector<int> isteps(Nlev, 0);
    Vector<IntVect> refRatios(Nlev-1,{AMREX_D_DECL(2, 2, 2)});
    amrex::WriteMultiLevelPlotfileHDF5SingleDset(outfile, Nlev, GetVecOfConstPtrs(fileData), plotVarNames,
                                   geoms, amrData.Time(), isteps, refRatios, hdf5_compression);

  }
  Finalize();
  return 0;
}
