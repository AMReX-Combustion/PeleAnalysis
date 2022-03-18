#include <util.H>

#include <AMReX_Print.H>

using namespace amrex;

static std::string
decodeStringFromFortran(const Vector<int>& coded,
                        int                length)
{
    std::string result;
    for (int i = 0; i < length; ++i)
        result += coded[i];
    return result;
}

namespace analysis_util {

  int
  NumSpecies()
  {
    return NUM_SPECIES;
  }
  
  int
  NumElements()
  {
    return NUM_ELEMENTS;
  }
  
  int
  NumReactions()
  {
    return NUM_REACTIONS;
  }
  
  Vector<std::string>
  GetSpecNames ()
  {
    Vector<std::string> names;
    CKSYMS_STR(names);
    return names;
  }

  Vector<std::string>
  GetElemNames ()
  {
    Vector<std::string> names;
    CKSYME_STR(names);
    return names;
  }

  amrex::Vector<int>
  ReactionsWithXonL(int ispec)
  {
    int max_num_spec, num_spec;
    Vector<int> specIDs, specNUs;
    int rid = 0; // Special flag to find max_num_spec
    CKINU(&rid,&max_num_spec,specIDs.data(),specNUs.data());
    int nreactions = analysis_util::NumReactions();
    Vector<int> rns;
    specIDs.resize(max_num_spec);
    specNUs.resize(max_num_spec);
    for (int j=1; j<=nreactions; ++j) {
      CKINU(&j,&num_spec,specIDs.data(),specNUs.data());
      if (num_spec > 0) {
        for (int is=0; is<num_spec; ++is) {
          if (specIDs[is]==ispec+1 && specNUs[is]<0) {
            rns.push_back(j-1);
          }
        }
      }
    }
    return rns;
  }

  amrex::Vector<int>
  ReactionsWithXonR(int ispec)
  {
    int max_num_spec, num_spec;
    Vector<int> specIDs, specNUs;
    int rid = 0; // Special flag to find max_num_spec
    CKINU(&rid,&max_num_spec,specIDs.data(),specNUs.data());
    int nreactions = analysis_util::NumReactions();
    Vector<int> rns;
    specIDs.resize(max_num_spec);
    specNUs.resize(max_num_spec);
    for (int j=1; j<=nreactions; ++j) {
      CKINU(&j,&num_spec,specIDs.data(),specNUs.data());
      if (num_spec > 0) {
        for (int is=0; is<num_spec; ++is) {
          if (specIDs[is]==ispec+1 && specNUs[is]>0) {
            rns.push_back(j-1);
          }
        }
      }
    }
    return rns;
  }

  amrex::Vector<std::pair<std::string,int> >
  specCoeffsInReactions(int ireac)
  {
    amrex::Vector<std::pair<std::string,int> > coeffs;

    int max_num_spec, num_spec;
    Vector<int> specIDs, specNUs;
    Vector<std::string> spec_names = analysis_util::GetSpecNames();
    int rid = 0; // Special flag to find max_num_spec
    CKINU(&rid,&max_num_spec,specIDs.data(),specNUs.data());
    specIDs.resize(max_num_spec);
    specNUs.resize(max_num_spec);
    int freac = ireac + 1;
    CKINU(&freac,&num_spec,specIDs.data(),specNUs.data());
    for (int is=0; is<num_spec; ++is) {
      coeffs.push_back(std::make_pair(spec_names[specIDs[is]-1],specNUs[is]));
    }
    return coeffs;
  }

  int
  IndexSpec (const std::string& specName)
  {
    Vector<std::string> speciesNames = GetSpecNames();
    for (int i=0; i<speciesNames.size(); i++) {
      if (specName == speciesNames[i]) {
        return i;
      }
    }
    return -1;
  }

  int
  IndexElem (const std::string& elemName)
  {
    Vector<std::string> elemNames = GetElemNames();
    for (int i=0; i<elemNames.size(); i++) {
      if (elemName == elemNames[i]) {
        return i;
      }
    }
    return -1;
  }
  
  amrex::Vector<int>
  GetReactionMap()
  {
    int nreactions = analysis_util::NumReactions();    
    Vector<int> rmap(nreactions);
    GET_RMAP(rmap.data());
    return rmap;
  }

  int
  NumElemXinSpecY(const std::string& elem,
                  const std::string& spec)
  {
    int nelements = analysis_util::NumElements();
    int nspecies = analysis_util::NumSpecies();
    Vector<int> ncf(nspecies * nelements);
    CKNCF(ncf.data());
    int ielem = analysis_util::IndexElem(elem);
    int ispec = analysis_util::IndexSpec(spec);
    return ncf[ispec*nelements + ielem];
  }

  Edge::Edge (const std::string& n1,
              const std::string& n2,
              const amrex::Vector<std::pair<int,amrex::Real> > rwl)
    : sp1(n1), sp2(n2), RWL(rwl) {}

  Edge::Edge (const std::string& n1,
              const std::string& n2,
              int reac,
              amrex::Real weight )
    : sp1(n1), sp2(n2)
  {
    RWL.push_back(std::pair<int,amrex::Real>(reac,weight));
  }

  int Edge::equivSign (const Edge& rhs) const
  {   
    if ( (sp1 == rhs.sp1) && (sp2 == rhs.sp2) ) {
      return 1;
    } else {
      if ( (sp1 == rhs.sp2) && (sp2 == rhs.sp1) )
        return -1;
    }
    return 0;
  }

  void Edge::combine (const Edge& rhs, const int   sgn)
  {
    if (sgn!=0)
    {            
      int oldSize = RWL.size();
      int delSize = rhs.RWL.size();
      RWL.resize(oldSize + delSize);
      for (int i=0; i<delSize; ++i) {
        RWL[oldSize+i] = std::pair<int,amrex::Real>(rhs.RWL[i].first,sgn*rhs.RWL[i].second);
      }
    }
  }

  bool
  Edge::touchesSp(const std::string& rhs) const
  {
    return ( (sp1 == rhs) || (sp2 == rhs) );
  }

  void Edge::reverse()
  {
    for (int i=0; i<RWL.size(); ++i) {
      RWL[i].second = - RWL[i].second;
    }
    std::swap(sp1,sp2);
  }

  const amrex::Vector<std::pair<int,amrex::Real> >&
  Edge::rwl () const
  {
    return RWL;
  }

  const std::string
  Edge::left() const
  {
    return sp1;
  }

  const std::string
  Edge::right() const
  {
    return sp2;
  }

  bool
  Edge::operator<(const Edge& rhs) const
  {
    return (this->left() == rhs.left() ?  this->right() < rhs.right() : this->left() < rhs.left());
  }

  std::ostream& operator<< (std::ostream& os, const Edge& e)
  {
    os << e.sp1 << " <=> " << e.sp2 << "  ( ";
    for (amrex::Vector<std::pair<int,amrex::Real> >::const_iterator it=e.RWL.begin(); it!=e.RWL.end(); ++it) {
      os << it->first << ":" << it->second << " ";
    }
    os << ")";
    return os;
  }

  Group::Group (const std::map<std::string,int>& eltCnts)
  {
    mEltCnts = eltCnts;
  }

  Group::Group (const Group& rhs)
  {
    mEltCnts = rhs.mEltCnts;
  }

  Group
  Group::operator- (const Group& rhs) const
  {
    Group val(*this);
    for (std::map<std::string,int>::const_iterator rit = rhs.mEltCnts.begin();rit!=rhs.mEltCnts.end();++rit)
    {
      std::map<std::string,int>::iterator fit = val.mEltCnts.find( rit->first );
      if ( fit != val.mEltCnts.end() ) {
        int res = fit->second - rit->second;
        if (res != 0) {
          fit->second = res;
        } else {
          val.mEltCnts.erase(fit);
        }
      }
      else {
        val.mEltCnts[rit->first] = -rit->second;
      }
    }
    return val;
  }

  Group
  Group::operator* (int rhs) const
  {
    Group val(*this);
    for (std::map<std::string,int>::iterator vit = val.mEltCnts.begin();vit!=val.mEltCnts.end();++vit) {
      vit->second = rhs*vit->second;
    }
    return val;
  }

  int
  Group::operator[](const std::string& id) const
  {
    std::map<std::string,int>::const_iterator fit = mEltCnts.find( id );
    if (fit!=mEltCnts.end()) {
      return fit->second;
    }
    return 0;
  }

  bool
  Group::sameSign () const
  {
    if (mEltCnts.size()>0)
    {
      std::map<std::string,int>::const_iterator it = mEltCnts.begin();
      if (it->second < 0) {
        for (++it; it!=mEltCnts.end(); ++it)
          if (it->second > 0)
            return false;
      } else {
        for (++it; it!=mEltCnts.end(); ++it) {
          if (it->second < 0) {
            return false;
          }
        }
      }
    }
    return true;
  }

  bool
  Group::contains (const std::string& id) const
  {
    return mEltCnts.find(id)!=mEltCnts.end();
  }
    
  amrex::Real
  Group::awt () const
  {
    if (AtomicWeight.size() == 0) {
      FillAtomicWeights(); // read up static atomic weight data
    }
    amrex::Real myWt = 0;
    for (std::map<std::string,int>::const_iterator it = mEltCnts.begin(); it!=mEltCnts.end(); ++it) {
      amrex::Real thisAwt = AtomicWeight[it->first];
      myWt += std::abs(it->second) * thisAwt;
    }
    return myWt;
  }

  int
  Group::size() const
  {
    int mySize = 0;
    for (std::map<std::string,int>::const_iterator it = mEltCnts.begin(); it!=mEltCnts.end(); ++it) {
      mySize += std::abs(it->second);
    }
    return mySize;
  }

  Group operator* (int n, const Group& g)
  {
    return g.operator*(n);
  }

  void
  Group::FillAtomicWeights () const
  {
    AtomicWeight["H"] = 1.00797;
    AtomicWeight["HE"] = 4.00260;
    AtomicWeight["LI"] = 6.93900;
    AtomicWeight["BE"] = 9.01220;
    AtomicWeight["B"] = 10.81100;
    AtomicWeight["C"] = 12.01115;
    AtomicWeight["N"] = 14.00670;
    AtomicWeight["O"] = 15.99940;
    AtomicWeight["F"] = 18.99840;
    AtomicWeight["NE"] = 20.18300;
    AtomicWeight["NA"] = 22.98980;
    AtomicWeight["MG"] = 24.31200;
    AtomicWeight["AL"] = 26.98150;
    AtomicWeight["SI"] = 28.08600;
    AtomicWeight["P"] = 30.97380;
    AtomicWeight["S"] = 32.06400;
    AtomicWeight["CL"] = 35.45300;
    AtomicWeight["AR"] = 39.94800;
    AtomicWeight["K"] = 39.10200;
    AtomicWeight["CA"] = 40.08000;
    AtomicWeight["SC"] = 44.95600;
    AtomicWeight["TI"] = 47.90000;
    AtomicWeight["V"] = 50.94200;
    AtomicWeight["CR"] = 51.99600;
    AtomicWeight["MN"] = 54.93800;
    AtomicWeight["FE"] = 55.84700;
    AtomicWeight["CO"] = 58.93320;
    AtomicWeight["NI"] = 58.71000;
    AtomicWeight["CU"] = 63.54000;
    AtomicWeight["ZN"] = 65.37000;
    AtomicWeight["GA"] = 69.72000;
    AtomicWeight["GE"] = 72.59000;
    AtomicWeight["AS"] = 74.92160;
    AtomicWeight["SE"] = 78.96000;
    AtomicWeight["BR"] = 79.90090;
    AtomicWeight["KR"] = 83.80000;
    AtomicWeight["RB"] = 85.47000;
    AtomicWeight["SR"] = 87.62000;
    AtomicWeight["Y"] = 88.90500;
    AtomicWeight["ZR"] = 91.22000;
    AtomicWeight["NB"] = 92.90600;
    AtomicWeight["MO"] = 95.94000;
    AtomicWeight["TC"] = 99.00000;
    AtomicWeight["RU"] = 101.07000;
    AtomicWeight["RH"] = 102.90500;
    AtomicWeight["PD"] = 106.40000;
    AtomicWeight["AG"] = 107.87000;
    AtomicWeight["CD"] = 112.40000;
    AtomicWeight["IN"] = 114.82000;
    AtomicWeight["SN"] = 118.69000;
    AtomicWeight["SB"] = 121.75000;
    AtomicWeight["TE"] = 127.60000;
    AtomicWeight["I"] = 126.90440;
    AtomicWeight["XE"] = 131.30000;
    AtomicWeight["CS"] = 132.90500;
    AtomicWeight["BA"] = 137.34000;
    AtomicWeight["LA"] = 138.91000;
    AtomicWeight["CE"] = 140.12000;
    AtomicWeight["PR"] = 140.90700;
    AtomicWeight["ND"] = 144.24000;
    AtomicWeight["PM"] = 145.00000;
    AtomicWeight["SM"] = 150.35000;
    AtomicWeight["EU"] = 151.96000;
    AtomicWeight["GD"] = 157.25000;
    AtomicWeight["TB"] = 158.92400;
    AtomicWeight["DY"] = 162.50000;
    AtomicWeight["HO"] = 164.93000;
    AtomicWeight["ER"] = 167.26000;
    AtomicWeight["TM"] = 168.93400;
    AtomicWeight["YB"] = 173.04000;
    AtomicWeight["LU"] = 174.99700;
    AtomicWeight["HF"] = 178.49000;
    AtomicWeight["TA"] = 180.94800;
    AtomicWeight["W"] = 183.85000;
    AtomicWeight["RE"] = 186.20000;
    AtomicWeight["OS"] = 190.20000;
    AtomicWeight["IR"] = 192.20000;
    AtomicWeight["PT"] = 195.09000;
    AtomicWeight["AU"] = 196.96700;
    AtomicWeight["HG"] = 200.59000;
    AtomicWeight["TL"] = 204.37000;
    AtomicWeight["PB"] = 207.19000;
    AtomicWeight["BI"] = 208.98000;
    AtomicWeight["PO"] = 210.00000;
    AtomicWeight["AT"] = 210.00000;
    AtomicWeight["RN"] = 222.00000;
    AtomicWeight["FR"] = 223.00000;
    AtomicWeight["RA"] = 226.00000;
    AtomicWeight["AC"] = 227.00000;
    AtomicWeight["TH"] = 232.03800;
    AtomicWeight["PA"] = 231.00000;
    AtomicWeight["U"] = 238.03000;
    AtomicWeight["NP"] = 237.00000;
    AtomicWeight["PU"] = 242.00000;
    AtomicWeight["AM"] = 243.00000;
    AtomicWeight["CM"] = 247.00000;
    AtomicWeight["BK"] = 249.00000;
    AtomicWeight["CF"] = 251.00000;
    AtomicWeight["ES"] = 254.00000;
    AtomicWeight["FM"] = 253.00000;
    AtomicWeight["D"] = 002.01410;
    AtomicWeight["E"] = 5.48578E-4;
  }

  std::ostream& operator<< (std::ostream& os, const Group& g)
  {
    os << "Group < ";
    for (std::map<std::string,int>::const_iterator it=g.mEltCnts.begin(); it!=g.mEltCnts.end(); ++it) {
        os << it->first << ":" << it->second << " ";
    }
    os << ">";
    return os;
  }

  std::list<Edge>
  getEdges (const std::string& trElt, int PrintVerbose, int HackSplitting)
  {
    std::list<Edge> edges;
    std::map<std::string,Group> groups;
    const amrex::Vector<std::string>& spNames = GetSpecNames();
    const amrex::Vector<std::string>& eltNames = GetElemNames();
    for (int i=0; i<spNames.size(); ++i) {
      const std::string sp = spNames[i];
      std::map<std::string,int> eltCnts;
      for (int j=0; j<eltNames.size(); ++j) {
        const std::string elt = eltNames[j];
        int m = NumElemXinSpecY(elt,sp);
        if (m>0) {
          eltCnts[elt] = m;
        }
        groups[sp] = Group(eltCnts);
      }
    }
    for (int r=0; r<NumReactions(); ++r) {
      const amrex::Vector<std::pair<std::string,int> >& coeffs = specCoeffsInReactions(r);

      std::map<std::string,int> net;
      for (int i=0; i<coeffs.size(); ++i) {
        const std::string& sp=coeffs[i].first;
        const int co=coeffs[i].second;
        std::map<std::string,int>::iterator it = net.find(sp);
        if (it!=net.end()) {
          net[sp] += co;
          if (net[sp]==0) {
            net.erase(it);
          }
        } else {
          net[sp] = co;
        }
      }

      std::map<std::string,int> reac, prod;
      for (std::map<std::string,int>::const_iterator it = net.begin(); it!=net.end(); ++it) {
        const std::string& sp = it->first;
        const int n = it->second;
        if (n < 0) {
          if (NumElemXinSpecY(trElt,sp)>0) {
            reac[sp] = -n;
          }
        } else {
          if (NumElemXinSpecY(trElt,sp)>0) {
            prod[sp] = n;
          }
        }
      }

      int LR = reac.size();
      int LP = prod.size();
        
      if (LR==0 || LP==0) { // no tr-containing species in reac
        continue;
      }
        
      if ((LR == 1) || (LP == 1)) {// exactly one tr-containing species either side
        for (std::map<std::string,int>::const_iterator rit=reac.begin(); rit!=reac.end(); ++rit) {
          const std::string& spcr = rit->first;
          const int cor = rit->second;

          for (std::map<std::string,int>::const_iterator pit=prod.begin(); pit!=prod.end(); ++pit){
            const std::string& spcp = pit->first;
            const int cop = pit->second;
            int w = std::min(cor*groups[spcr][trElt],cop*groups[spcp][trElt]);
            edges.push_back(Edge(spcr,spcp,r,w));
          }
        }
        continue;
      }

      if ( (LR==2) && (LP==2) ) { // two tr-containing species each side
        std::map<std::string,int>::const_iterator r0 = reac.begin();
        std::map<std::string,int>::const_iterator r1 = r0; r1++;
        std::map<std::string,int>::const_iterator p0 = prod.begin();
        std::map<std::string,int>::const_iterator p1 = p0; p1++;
 
        const std::string& rs0 = r0->first;
        const std::string& rs1 = r1->first;
        const std::string& ps0 = p0->first;
        const std::string& ps1 = p1->first;

        int rc0 = r0->second;
        int rc1 = r1->second;
        int pc0 = p0->second;
        int pc1 = p1->second;

        Group b0(pc0 * groups[ps0] - rc0 * groups[rs0]);
        Group b1(pc1 * groups[ps1] - rc0 * groups[rs0]);
        int pick = 0;

        // HACK
        if ((HackSplitting==1) && (trElt=="H") && (r==61)) {
          pick = 0;
          if (PrintVerbose>0) {
            Print() << "...resorting to HACK!!" << std::endl;
          }
        } else {
          if ((HackSplitting==1) && (trElt=="H") && (r==310)) {
            pick = 0;
            if (PrintVerbose>0) {
              Print() << "...resorting to HACK!!" << std::endl;
            }
          } else {
            if (b0.sameSign() && b1.sameSign()) {

              if (b1.size() < b0.size()) {
                pick = 1;
                if (PrintVerbose>0) {
                  Print() << "...resorting to atom cnt" << std::endl;
                }
              }
              if (b1.size() == b0.size()) {
                if (b0.awt() > b1.awt()) {
                  pick = 1;
                  if (PrintVerbose>0) {
                    Print() << "...resorting to tot awt" << std::endl;
                  }
                } else if ( (PrintVerbose>0) && (b0.awt() < b1.awt()) ) {
                  Print() << "...resorting to tot awt" << std::endl;
                }
              }
            } else {
              if (b1.sameSign()) {
                pick = 1;
              }
            }
          }
        }
 
        int nR0 = rc0*groups[rs0][trElt];
        int nR1 = rc1*groups[rs1][trElt];
        int nP0 = pc0*groups[ps0][trElt];
        int nP1 = pc1*groups[ps1][trElt];

        if (pick == 0) {
          edges.push_back(Edge(rs0,ps0,r,std::min(nR0,nP0)));
          if (nP0<nR0) {
            edges.push_back(Edge(rs0,ps1,r,nR0-nP0));
          }
          
          edges.push_back(Edge(rs1,ps1,r,std::min(nR1,nP1)));
          if (nR0<nP0) {
            edges.push_back(Edge(rs1,ps0,r,nP0-nR0));
          }
        } else {
          edges.push_back(Edge(rs0,ps1,r,std::min(nR0,nP1)));
          if (nP1<nR0) {
            edges.push_back(Edge(rs0,ps0,r,nR0-nP1));
          }
 
          edges.push_back(Edge(rs1,ps0,r,std::min(nR1,nP0)));
          if (nR0<nP1) {
            edges.push_back(Edge(rs1,ps1,r,nP1-nR0));
          }
        }
        continue;
      }

      Print() << "Cannot decompose rxn: " << r << " " << LR << " " << LP << std::endl;
    }

    // Uniquify/combine edges
    std::list<Edge> uEdges;
    while (edges.size()!=0) {
      const Edge oe(edges.front()); edges.pop_front();
      bool foundIt = false;
      // There's probably a faster way to search this list...
      for (std::list<Edge>::iterator it=uEdges.begin(); !foundIt && it!=uEdges.end(); ++it) {
        int sgn = it->equivSign(oe);
        if (sgn!=0) {
          it->combine(oe,sgn);
          foundIt = true;
        }
      }
      if (!foundIt) {
        uEdges.push_back(oe);
      }
    }
    return uEdges;
  }
}


