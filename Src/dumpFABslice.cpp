#include <string>
#include <iostream>

#include <AMReX_ParmParse.H>
#include <AMReX_DataServices.H>

using namespace amrex;

static
void
print_usage (int,
             char* argv[])
{
  std::cerr << "usage:\n";
  std::cerr << argv[0] << " infile infile=<s> [options] \n\tOptions:\n";
  std::cerr << "\t     varNames=<s,s,...> variables to slice (will output individual planes), (OPT, will output all by default)\n";
  std::cerr << "\t     dir=<i> direction to take slice (only needed for 3D)\n";
  std::cerr << "\t     num=<i> slice number to take (only needed for 3D)\n";

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
  Initialize(argc,argv);
  {
    if (argc < 2)
      print_usage(argc,argv);
    //declare ParmParse
    ParmParse pp;

    if (pp.contains("help"))
      print_usage(argc,argv);

    if (pp.contains("verbose"))
      AmrData::SetVerbose(true);

    std::string plotFileName;
    pp.get("infile",plotFileName);
    DataServices::SetBatchMode();
    Amrvis::FileType fileType(Amrvis::NEWPLT);
    DataServices dataServices(plotFileName, fileType);
    if( ! dataServices.AmrDataOk()) {
      DataServices::Dispatch(DataServices::ExitRequest, NULL);
    }
#if AMREX_SPACEDIM == 2
    constexpr int dir = 2;
    constexpr int num = 0;
#else
    int dir,num;
    pp.get("dir",dir);
    pp.get("num",num);
#endif

    AmrData& amrData = dataServices.AmrDataRef();
    int finestLevel = amrData.FinestLevel();
    pp.query("finestLevel",finestLevel);
    Box domainBox = amrData.ProbDomain()[finestLevel];
    std::string fabName;
    int numVars = pp.countval("varNames");
    Vector<std::string> varNames;
    if (numVars > 0) {
      varNames.resize(numVars);
      pp.getarr("varNames",varNames);
      for (int n = 0; n < numVars; n++) {
	DataServices::Dispatch(DataServices::DumpSlicePlaneOneVar,&dataServices,dir,num,varNames[n]);	
      }
    } else {
      DataServices::Dispatch(DataServices::DumpSlicePlaneAllVars,&dataServices,dir,num);	
    }    
  }  
  amrex::Finalize();
  return 0;
}
