#include <string>
#include <iostream>

#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include <AMReX_DataServices.H>
#include <AMReX_PlotFileUtil.H>

#include <AMReX_BLFort.H>

using namespace amrex;

int main(int argc, char *argv[])
{
  amrex::Initialize(argc, argv);

  int nProcs = ParallelDescriptor::NProcs();
  //if (nProcs>1)
  //  amrex::Error("Not yet fully-implemented in parallel!");

  int myProc = ParallelDescriptor::MyProc();

  ParmParse pp;
  
  bool verbose(false);
  if(pp.contains("verbose") || (pp.contains("v"))) {
    verbose = true;
    AmrData::SetVerbose(true);
  }

  if (ParallelDescriptor::IOProcessor()) {
    verbose = true;
    AmrData::SetVerbose(true);
  }
  
  DataServices::SetBatchMode();
  Amrvis::FileType fileType(Amrvis::NEWPLT);

  // read in list of fabs to load
  //int nFiles=pp.countval("infile");
  int start,end,interval,nfiles;  
  std::string infile_head, infile_tail;

  /*
    std::string infile[nFiles];
    if (verbose) std::cout << "infiles =";
    for (int iFile=0; iFile<nFiles; iFile++) {
      pp.get("infile",infile[iFile],iFile);
      if (verbose) std::cout << " " << infile[iFile]; 
    }
    if (verbose) std::cout << std::endl;
     */
  pp.get("infile_head",infile_head);
  pp.get("infile_tail",infile_tail);
  pp.get("start",start);
  pp.get("end",end);
  pp.get("interval",interval);
  nfiles = (end-start)/interval;
  // read plot file name to output
  std::string outfile;
  pp.get("outfile",outfile);
  if (verbose) std::cout << "outfile = " << outfile << std::endl;

  // read a time to assign
  Real time(0.);
  pp.query("time",time);
  if (verbose) std::cout << "time = " << time << std::endl;

  // read in the variable names to use
  int nVars=pp.countval("names");
  Vector<std::string> names; names.resize(nVars);
  if (verbose) std::cout << "variables =";
  for (int iVar=0; iVar<nVars; iVar++) {
    pp.get("names",names[iVar],iVar);
    if (verbose) std::cout << " " << names[iVar]; 
  }
  if (verbose) std::cout << std::endl;

  //
  // Read fab data (only using header at this stage)
  //
  if (verbose) std::cout << "Reading fab..." << std::endl;
  std::ifstream ifs;
  std::string start_str = std::to_string(start);
  if (start < 100000) {
    start_str = "0"+start_str;
  }
		      
  std::string infiletest = infile_head + start_str + infile_tail;
  amrex::Print() << "Getting geometry from file: "+infiletest << std::endl;
  FArrayBox fab;

  //Box* boxptr = new Box[1];
  /*Real* probLo = malloc(3*sizeof(Real));
  probLo[0] = 0.0; probLo[1] = 0.0; probLo[2] = 0.0;
  Real* probHi = malloc(3*sizeof(Real));
  probHi[0] = 0.0; probHi[1] = 0.0; probHi[2] = 0.0;
  Real* nx = malloc(3*sizeof(int));
  nx[0] = 0; nx[1] = 0; nx[2] = 0;*/
  Vector<Real> probLo = {0.0,0.0,0.0};
  Vector<Real> probHi = {0.0,0.0,0.0};
  Vector<int> nx = {0,0,0};
  
  if (ParallelDescriptor::IOProcessor()) {
    ifs.open(infiletest);
    fab.readFrom(ifs);
    ifs.close();
    if (verbose) std::cout << "   ... done." << std::endl;


    std::cout << "fab.nComp = " << fab.nComp() << std::endl;
    if (fab.nComp()!=nVars) {
      amrex::Error("Mismatch fab.nComp != nVars");
    }
    
    // figure out how big the full domain box needs to be
    const Box& box = fab.box();
    std::cout << "fab.box().length = " << box.length() << std::endl;
    //boxptr[0] = boxloc;
    //ParallelDescriptor::Bcast(boxptr,1);
    //}
  //Box box = boxptr[0];

  
  // number of cells in each direction
  // (assume stacking in z)
    //int nx[3];
    nx[0] = box.length(0);
    nx[1] = box.length(1);
    nx[2] = nfiles;
    int pCells=nx[0]*nx[1];
    long nCells=nx[0]*nx[1]*nx[2];
    
    if (verbose) std::cout << "cells = "
			   << nx[0] << " " << nx[1] << " " << nx[2] << " "
			   << " (" << nCells << ")" << std::endl;
    
    // assign a physical size
    int calcz;
    if (pp.countval("probLo")==3) {
      for (int i=0; i<3; i++) {
	pp.get("probLo",probLo[i],i);
      }
      calcz = 0; 
    } else if (pp.countval("probLo") == 2) {
      for (int i=0; i<2; i++) {
	pp.get("probLo",probLo[i],i);
      }
      calcz = 1;
    } else {
      amrex::Error("Need 3 values for probLo");
    }
    
    if (pp.countval("probHi")==3) {
      for (int i=0; i<3; i++) {
	pp.get("probHi",probHi[i],i);
      }
    } else {
      amrex::Error("Need 3 values for probHi");
    }
    Real dx[3];
    if (calcz) {
      for (int i=0; i<2; i++) {
	dx[i] = (probHi[i]-probLo[i])/(Real)nx[i];
      }
      dx[2] = dx[0];
      probLo[2] = -dx[2]*nx[2];
    } else {
      for (int i=0; i<3; i++) {
	dx[i] = (probHi[i]-probLo[i])/(Real)nx[i];
      }
    }
  }

  ParallelDescriptor::ReduceRealSum(probLo.data(),3);
  ParallelDescriptor::ReduceRealSum(probHi.data(),3);
  ParallelDescriptor::ReduceIntSum(nx.data(),3);

  //if (verbose)
  Print() << "probLo = "  << probLo[0] << " " << probLo[1] << " " << probLo[2] << std::endl;
  //if (verbose)
  Print() << "probHi = "  << probHi[0] << " " << probHi[1] << " " << probHi[2] << std::endl;
    
  //if (verbose)
  Print() << "nx = "   << nx[0] << " " << nx[1] << " " << nx[2] << std::endl;

  
  RealBox rb(probLo.data(),probHi.data()); // make real box for geometry 
  Vector<int> is_per(AMREX_SPACEDIM,0); //hard code to no periodicity for the minute
    
  
  //
  // Let's try a slab domain decomposition
  //

  Vector<int> plovec = {0,0,0};
  Vector<int> phivec = {nx[0]-1,nx[1]-1,nx[2]-1};
  
  IntVect     pdLo(plovec);
  IntVect     pdHi(phivec);
  Box         probDomain(pdLo,pdHi);
  int coord = 0; //hard code cartesian
  Geometry geoms(probDomain, &rb, coord, &(is_per[0]));
  // make a box for each file
  int         nBoxes(nfiles);
  BoxArray    domainBoxArray(nBoxes);

  // Make some the slabs
  Box         tempBox(probDomain);    
  Vector<int> tempBoxSmall(nBoxes,0);
  Vector<int> tempBoxBig(nBoxes,0);

  if (verbose)
    std::cout << "Slab decomposition:" << std::endl;
  for (int iBox=0; iBox<nBoxes; iBox++) {
    tempBoxSmall[iBox] = probDomain.smallEnd(2) + iBox;
    tempBoxBig[iBox]   = tempBoxSmall[iBox];
    if (verbose) {
      std::cout << "   iBox / small / big: "
		<< iBox << " / "
		<< tempBoxSmall[iBox] << " / "
		<< tempBoxBig[iBox] << std::endl;
    }
  }    
  for (int iBox=0; iBox<nBoxes; iBox++) {
    tempBox.setSmall(Amrvis::ZDIR, tempBoxSmall[iBox]);
    tempBox.setBig(Amrvis::ZDIR, tempBoxBig[iBox]); 
    domainBoxArray.set(iBox, tempBox);
  }
  
  // And now the distibution mapping
  /*
  Vector<int> pmap(nBoxes);
  for (int iBox=0; iBox<nBoxes; iBox++)
    pmap[iBox] = iBox%nProcs; // maybe think about hacking this to parallelise this step
  DistributionMapping domainDistMap(pmap);
  */
  DistributionMapping domainDistMap(domainBoxArray);
  
  MultiFab *mf;
  mf = new MultiFab(domainBoxArray,domainDistMap,nVars,0);
  //
  // Populate data
  //
  if (verbose)
    std::cout << "Populating data:" << std::endl;

  // how do we loop over the boxes in the multifab properly?
  //int iFile=0;
  for (MFIter mfi(*mf); mfi.isValid(); ++mfi) {
    
    // destination fab
    FArrayBox& myFab = (*mf)[mfi];
    
    // load data
    //std::cout << "MF idx = "+std::to_string(mfi.index()) << std::endl;
    int iFile=nfiles-myFab.smallEnd()[2]-1;
    int fileNum = start+iFile*interval;
    std::string filenum_str = std::to_string(fileNum);
    if (fileNum < 100000) {
      filenum_str = "0"+filenum_str;
    }
    
    std::string infile_local = infile_head+filenum_str+infile_tail;//infile[nFiles-iFile-1].c_str()
    std::cout << "iFile = " << iFile << ": " <<  infile_local << std::endl;
    ifs.open(infile_local);
    fab.readFrom(ifs);
    ifs.close();
    Vector<int> shiftvec = {0,0,myFab.smallEnd()[2]};
    IntVect shift(shiftvec);
    fab.shift(shift);
    const Box& inBox = fab.box();
    //amrex::Print() << "inBox = " << inBox << std::endl;
    //amrex::Print() << "myFab box= " << myFab.box() << std::endl;
    for (int dir=0; dir<2; dir++) {
      if (inBox.length(dir)!=nx[dir]) {
	std::cerr << "file = " << infile_local << std::endl;
	std::cerr << "inBox.length(" << dir << ") = " << inBox.length(dir) << std::endl;
	amrex::Error("inBox.length mismatch!");
      }
    }

    // copy data
    myFab.copy(fab);
    
    //iFile++;
  }

  // write the output plotfile
  // should be able to replace with modern call to writeplotfile

  if (verbose) {
    std::cout << "*** writing plotfile " << std::endl;
  }
  int levelSteps;
  WriteSingleLevelPlotfile(outfile,*mf,names, geoms,time,levelSteps);
  amrex::Finalize();
  return 0;
}
















