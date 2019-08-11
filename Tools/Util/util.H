#ifndef util_H
#define util_H

#include <AMReX_Vector.H>
#include <map>
#include <list>

namespace analysis_util
{
  int NumElements();
  int NumSpecies();
  int NumReactions();
  int IndexSpec (const std::string& specName);
  int IndexElem (const std::string& elemName);

  amrex::Vector<int> GetReactionMap();
  
  amrex::Vector<std::string> GetElemNames();
  amrex::Vector<std::string> GetSpecNames();

  amrex::Vector<int> ReactionsWithXonL(int ispec);
  amrex::Vector<int> ReactionsWithXonR(int ispec);
  amrex::Vector<std::pair<std::string,int> > specCoeffsInReactions(int ireac);
  int NumElemXinSpecY(const std::string& elem,
                      const std::string& spec);

  //
  // Data structure for generating chemical path diagrams,
  //   an Edge represents the transfer of an atom from one chemical
  //   species, sp1, to another, sp2, which occurs through a list
  //   of reactions.  For each reaction, n of these atoms are
  //   transferred (RWL is a list of pairs of rxn ID to atom transfer cnt)
  //
  class Edge
  {
  public:
    friend std::ostream& operator<< (std::ostream& os, const Edge& e);
        
    Edge (const std::string& n1,
          const std::string& n2,
          const amrex::Vector<std::pair<int,amrex::Real> > rwl);
        
    Edge (const std::string& n1,
          const std::string& n2,
          int reac,
          amrex::Real weight );
        
    int equivSign (const Edge& rhs) const;
    void combine (const Edge& rhs, int sgn);
    bool touchesSp(const std::string& rhs) const;
    void reverse();
    const amrex::Vector<std::pair<int,amrex::Real> >& rwl () const;
    const std::string left() const;
    const std::string right() const;
    bool operator< (const Edge& rhs) const;
  private:
    std::string sp1, sp2;
    amrex::Vector<std::pair<int,amrex::Real> > RWL; // RateWeightList, each pair is (reac,coeff)
  };

  // 
  // Helper class for building chem edges.  A group is a list of constituent 
  // atoms, and this class allows a coupleof useful math operations on these
  // groups.  
  //
  class Group
  {
  public:
    friend std::ostream& operator<< (std::ostream& os, const Group& g);
        
    Group () {}
    Group (const std::map<std::string,int>& eltCnts);
    Group (const Group& rhs);
    Group operator- (const Group& rhs) const;
    Group operator* (int rhs) const;
    int operator[](const std::string& id) const;
    bool sameSign() const;
    bool contains(const std::string& id) const;
    amrex::Real awt() const;
    int size() const;
  private:
    void FillAtomicWeights () const;
    std::map<std::string,int> mEltCnts;
    mutable std::map<std::string,amrex::Real> AtomicWeight;
  };

  // 
  // Compute edges from chem mechanism
  //
  std::list<Edge> getEdges (const std::string& trElt,
                            int PrintVerbose=0,
                            int HackSplitting=1);
}
#endif
