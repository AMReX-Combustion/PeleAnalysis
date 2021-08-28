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
#include <PelePhysics.H>
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

    for (int i=0; i<nspecies; ++i) {
      Print() << spec_names[i] << " = { ";
      for (int j=0; j<nelements; ++j) {
        int ne = NumElemXinSpecY(elem_names[j],spec_names[i]);
        if (ne > 0) {
          Print() << elem_names[j] << ":" << ne << " "; 
        }
      }
      Print() << "}\n";
    }

    int nreactions = NumReactions();
    Print() << "NumReactions: " << nreactions << std::endl;

    Vector<int> rns(nreactions);
    for (int i=0; i<nspecies; ++i) {
      rns = ReactionsWithXonL(i);
      Print() << "Spec, rnsL, rnsR: " << spec_names[i] <<  ": ";
      Print() << "{ ";
      for (int j=0; j<rns.size(); ++j) {
        Print() << rns[j] << " ";
      }
      Print() << "}";
      rns = ReactionsWithXonR(i);
      Print() << " { ";
      for (int j=0; j<rns.size(); ++j) {
        Print() << rns[j] << " ";
      }
      Print() << "}" << std::endl;
    }

    Vector<int> rmap = GetReactionMap();
    Print() << "rmap = { ";
    for (int j=0; j<rmap.size(); ++j) {
      Print() << rmap[j] << " ";
    }
    Print() << "}\n";

    // Build revese reaction map
    Vector<int> rrmap(nreactions);
    for (int i=0; i<nreactions; ++i) {
      rrmap[rmap[i]] = i;
    }
    Print() << "rrmap = { ";
    for (int j=0; j<rrmap.size(); ++j) {
      Print() << rrmap[j] << " ";
    }
    Print() << "}\n";

    for (int j=0; j<nreactions; ++j) {
      auto sc = specCoeffsInReactions(rmap[j]);
      Print() << "rn[" << j << "] (" << rmap[j] << ") = { ";
      for (int i=0; i<sc.size(); ++i) {
        Print() << sc[i].first << ":" << sc[i].second << " ";
      }
      Print() << "}\n";
    }

    std::string trElem = "C";
    pp.query("trElem",trElem);
    AMREX_ALWAYS_ASSERT(IndexElem(trElem)>=0 && IndexElem(trElem)<NumElements());
    auto edges = getEdges(trElem,1,1);
    for (auto edge : edges) {
      Print() << edge << std::endl;
    }


  }
  Finalize();
  return 0;
}
