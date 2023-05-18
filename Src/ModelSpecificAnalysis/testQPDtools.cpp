#include <string>
#include <iostream>
#include <set>

#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include <AMReX_DataServices.H>
#include <AMReX_BCRec.H>
#include <AMReX_Interpolater.H>
#include <AMReX_WritePlotFile.H>

#include <AMReX_BLFort.H>
#include <PelePhysics.H>
#include <util.H>

using namespace amrex;
using namespace analysis_util;

int
main (int   argc,
      char* argv[])
{
  Initialize(argc,argv);
  {
    ParmParse pp;

    Print() << " ==> Element list \n";
    Vector<std::string> elem_names;
    CKSYME_STR(elem_names);
    Print() << "{ ";
    for (int i=0; i<NUM_ELEMENTS; ++i) {
      Print() << elem_names[i] << " ";
    }
    Print() << "}" << std::endl;

    Print() << "\n ==> Species list \n";
    Vector<std::string> spec_names;
    CKSYMS_STR(spec_names);
    Print() << "{ ";
    for (int i=0; i<NUM_SPECIES; ++i) {
      Print() << spec_names[i] << " ";
    }
    Print() << "}" << std::endl;

    Print() << "\n ==> Species composition \n";
    for (int i=0; i<NUM_SPECIES; ++i) {
      Print() << spec_names[i] << " = { ";
      for (int j=0; j<NUM_ELEMENTS; ++j) {
        int ne = NumElemXinSpecY(j,i);
        if (ne > 0) {
          Print() << elem_names[j] << ":" << ne << " "; 
        }
      }
      Print() << "}\n";
    }

    Print() << "\n ==> NumReactions: " << NUM_REACTIONS << std::endl;

    Print() << "\n ==> Species in reactions \n";
    for (int i=0; i<NUM_SPECIES; ++i) {
      Vector<int> rnsL;
      getReactWithXOnL(rnsL,i);
      Print() << "Spec, rnsL, rnsR: " << spec_names[i] <<  ": ";
      Print() << "{ ";
      for (int j=0; j<rnsL.size(); ++j) {
        Print() << rnsL[j] << " ";
      }
      Print() << "}";
      Vector<int> rnsR;
      getReactWithXOnR(rnsR,i);
      Print() << " { ";
      for (int j=0; j<rnsR.size(); ++j) {
        Print() << rnsR[j] << " ";
      }
      Print() << "}" << std::endl;
    }

    Print() << "\n ==> RMAP \n";
    Vector<int> rmap = GetReactionMap();
    Print() << "{ ";
    for (int j=0; j<rmap.size(); ++j) {
      Print() << rmap[j] << " ";
    }
    Print() << "}\n";

    // Build revese reaction map
    Vector<int> rrmap(NUM_REACTIONS);
    for (int i=0; i<NUM_REACTIONS; ++i) {
      rrmap[rmap[i]] = i;
    }
    Print() << "\n ==> RRMAP \n";
    Print() << "{ ";
    for (int j=0; j<rrmap.size(); ++j) {
      Print() << rrmap[j] << " ";
    }
    Print() << "}\n";

    Print() << "\n ==> Stoich coeffs \n";
    for (int j=0; j<NUM_REACTIONS; ++j) {
      auto sc = specCoeffsInReactions(rmap[j]);
      Print() << "rn[" << j << "] (" << rmap[j] << ") = { ";
      for (int i=0; i<sc.size(); ++i) {
        Print() << sc[i].first << ":" << sc[i].second << " ";
      }
      Print() << "}\n";
    }

    Print() << "\n ==> Edges \n";
    std::string trElem = "C";
    pp.query("trElem",trElem);
    AMREX_ALWAYS_ASSERT(IndexElem(trElem)>=0 && IndexElem(trElem)<NUM_ELEMENTS);
    auto edges = getEdges(trElem,1,1);
    for (auto edge : edges) {
      Print() << edge << std::endl;
    }
  }
  Finalize();
  return 0;
}
