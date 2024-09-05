#include <string>
#include <iostream>

#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include <AMReX_DataServices.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_PlotFileUtil.H>

using namespace amrex;

static void print_usage (int, char* argv[])
{
  std::cerr << "usage:\n";
  std::cerr << argv[0] << " infiles=<> outfile=<> vars=<>" ;
  exit(1);
}

int main (int argc, char* argv[])
{
  amrex::Initialize(argc,argv);

  if (argc < 2)
    print_usage(argc,argv);
  
  ParmParse pp;

  // get infile names and count
  int nfiles(pp.countval("infiles"));
  Vector<std::string> infiles; infiles.resize(nfiles);
  pp.getarr("infiles",infiles);
  
  // get and count variables to copy
  int nvars(pp.countval("vars"));
  int id_comp_last = 0;
  Vector<std::string> vars; vars.resize(nvars);
  pp.getarr("vars",vars);  
  Vector<std::string> new_vars = vars;
  Vector<std::string> names;
  
  // outfile name
  std::string outfile;
  pp.get("outfile",outfile);
  
  DataServices::SetBatchMode();
  Amrvis::FileType fileType(Amrvis::NEWPLT);
  
    // setting up reading for pltfile
  DataServices dataServices0(infiles[0], fileType);
  
  // check if plotfile is good
  if (!dataServices0.AmrDataOk())
    DataServices::Dispatch(DataServices::ExitRequest, NULL);
  
  AmrData& amrData0 = dataServices0.AmrDataRef();
  
  // getting the finest level (or whatever user sets)
  int finestLevel = -1;
  finestLevel = amrData0.FinestLevel();
  pp.query("finestLevel",finestLevel);
  int Nlev = finestLevel + 1;

  // setting up periodicity
  Vector<int> is_per(AMREX_SPACEDIM,1);
  pp.queryarr("is_per",is_per,0,AMREX_SPACEDIM);
  Print() << "Periodicity assumed for this case: ";
  for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
    Print() << is_per[idim] << " ";
  }
  Print() << "\n";
  
  // setting up pltfile boxes
  RealBox rb(&(amrData0.ProbLo()[0]), 
	     &(amrData0.ProbHi()[0]));
  Vector<MultiFab*> fileData(Nlev);
  Vector<Geometry> geoms(Nlev);
  int coord = 0;
  for (int lev=0; lev<Nlev; ++lev)
    {
      const DistributionMapping dm(amrData0.boxArray(lev));
      geoms[lev] = Geometry(amrData0.ProbDomain()[lev],&rb,coord,&(is_per[0]));
      fileData[lev] = new MultiFab(amrData0.boxArray(lev),dm,nvars,0);
    }
  
  // loop over all files
  for (int iFile = 0; iFile < nfiles; iFile++) {
    Print() << "Loading plotfile: " << infiles[iFile] << std::endl;
    
    // set up for reading pltfile
    DataServices dataServices(infiles[iFile], fileType);
    
    // check if file is okay
    if (!dataServices.AmrDataOk())
      DataServices::Dispatch(DataServices::ExitRequest, NULL);
    
    AmrData& amrData = dataServices.AmrDataRef();
    
    // check if finest level is all the same
    int usrfinestLevel = -1;
    pp.query("finestLevel",usrfinestLevel);
    if (amrData.FinestLevel() != amrData0.FinestLevel() && usrfinestLevel != -1) {
      Print() << "Warning, finest levels for all pltfiles are not the same \n";
      DataServices::Dispatch(DataServices::ExitRequest, NULL);
    }

    const Vector<std::string>& plotVarNames = amrData.PlotVarNames();
    
    Vector<int> idcomp;
    for (int i = 0; i < plotVarNames.size(); ++i) {  // loop though current pltfile variable names
      for (int j = 0; j < new_vars.size(); ++j) {   // loop though remaining search variable names
	if (plotVarNames[i] == new_vars[j]) {       // if search name = current pltfile variable name
	  idcomp.push_back(i);                       // sets index location of plotVarName that matches
	  new_vars.erase(new_vars.begin() + j);    // delets search variable name so we do not get repeats
	}
      }
    }
    
    // copys data
    if (idcomp.size() > 0) {   // makes sure there a variable we want to copy
      for (int lev=0; lev<Nlev; ++lev) {
	for (int i=0; i<idcomp.size(); ++i) {
	  fileData[lev]->ParallelCopy(amrData.GetGrids(lev,idcomp[i]),0,id_comp_last+i,1);	    
	}
      }
      // set the correct name with the correct location
      for (int i=0; i<idcomp.size(); ++i)
	names.push_back(plotVarNames[idcomp[i]]);
    }
    
    // updates ids of comps
    id_comp_last += idcomp.size();
    idcomp.clear();
  }
  
  if (new_vars.size() > 0) {
    Print() << "Error the following comps were not found: \n";
    for (int i = 0; i < new_vars.size(); ++i)
      Print() << new_vars[i] << std::endl;
    DataServices::Dispatch(DataServices::ExitRequest, NULL);	      
  }
  
  // write pltfile
  Vector<int> isteps(Nlev, 0);
  Vector<IntVect> refRatios(Nlev-1,{AMREX_D_DECL(2, 2, 2)});
  amrex::WriteMultiLevelPlotfile(outfile, Nlev, GetVecOfConstPtrs(fileData), names, geoms, 0.0, isteps, refRatios);
  
  amrex::Finalize();
  return 0;
}
