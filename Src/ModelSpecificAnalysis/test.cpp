#include <string>
#include <iostream>
#include <set>

#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include <AMReX_DataServices.H>
#include <AMReX_BCRec.H>
#include <AMReX_Interpolater.H>
#include <WritePlotFile.H>

#include <AMReX_BLFort.H>
#include <mechanism.h>
#include <chemistry_file.H>
#include <util.H>
#include <util_F.H>

using namespace amrex;
using namespace analysis_util;

int
main (int   argc,
      char* argv[])
{
  Initialize(argc,argv);
  {
    init_mech();
    
    ParmParse pp;

    Vector<std::string> elem_names = GetElemNames();
    int nelements = NumElements();
    Print() << "{ ";
    for (int i=0; i<nelements; ++i) {
      Print() << elem_names[i] << " ";
    }
    Print() << "}" << std::endl;

    Vector<std::string> spec_names = GetSpecNames();
    int nspecies = NumSpecies();
    Print() << "{ ";
    for (int i=0; i<nspecies; ++i) {
      Print() << spec_names[i] << " ";
    }
    Print() << "}" << std::endl;
  }
  Finalize();
  return 0;
}
