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

    Vector<std::string> elem_names = GetElemNames();
    int nelements = elem_names.size();
    Print() << "Elements = { ";
    for (int i=0; i<nelements; ++i) {
      Print() << elem_names[i] << " ";
    }
    Print() << "}" << std::endl;

    Vector<std::string> spec_names = GetSpecNames();
    int nspecies = spec_names.size();
    Print() << "Species = { ";
    for (int i=0; i<nspecies; ++i) {
      Print() << spec_names[i] << " ";
    }
    Print() << "}" << std::endl;
    int nreactions = NumReactions();
    Print() << "NumReactions: " << nreactions << std::endl;
    Print() << "\n";

    Print() << "Elemental compositions: " << std::endl;
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
    Print() << "\n";

    Vector<int> rmap = GetReactionMap();
    Print() << "Reaction map = { ";
    for (int j=0; j<rmap.size(); ++j) {
      Print() << rmap[j] << " ";
    }
    Print() << "}\n\n";

    Print() << "For each species, reactions on left and right:" << std::endl;

    for (int i=0; i<nspecies; ++i) {
      Print() << spec_names[i] <<  ": ";
      Vector<int> rns = ReactionsWithXonL(i);
      Print() << "{ ";
      for (int j=0; j<rns.size(); ++j) {
        Print() << rns[j] << " (orig: " << rmap[rns[j]]+1 << ") ";
      }
      Print() << "}";
      rns = ReactionsWithXonR(i);
      Print() << " { ";
      for (int j=0; j<rns.size(); ++j) {
        Print() << rns[j] << " ";
      }
      Print() << "}" << std::endl;
    }
    Print() << "\n";

    std::string trElem = "C";
    pp.query("trElem",trElem);
    AMREX_ALWAYS_ASSERT(IndexElem(trElem)>=0 && IndexElem(trElem)<NumElements());
    Print() << "Edges: " << std::endl;
    auto edges = getEdges(trElem,1,1);
    for (auto edge : edges) {
      Print() << edge << std::endl;
    }
  }
  Finalize();
  return 0;
}
