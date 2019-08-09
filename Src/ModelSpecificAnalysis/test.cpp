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
    for (int i=0; i<nelements; ++i) {
      Print() << i << " " << elem_names[i] << std::endl;
    }

    Vector<std::string> spec_names = GetSpecNames();
    int nspecies = NumSpecies();
    for (int i=0; i<nspecies; ++i) {
      Print() << i << " " << spec_names[i] << std::endl;
    }
  }
  Finalize();
  return 0;
}
