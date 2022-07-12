#ifndef BL_STREAMDATA_H
#define BL_STREAMDATA_H

#include <AMReX_AmrData.H>


class StreamData
  : protected amrex::AmrData
{
public:
  struct MLloc
  {
    MLloc() :
      amr_lev(-1), box_idx(-1), pt_idx(-1) {}
    MLloc(int lev,int box,int pt) :
      amr_lev(lev), box_idx(box), pt_idx(pt) {}
    int amr_lev, box_idx, pt_idx;
  };

  StreamData(const std::string& file);
  bool ReadInfo();
  int NElts() const { return num_global_elements; }
  int NLines() const { return nLines; }
  const amrex::BoxArray &boxArray(int level) const;
  const amrex::Vector<std::string>& ComponentNames() const {return plotVars;}
  const std::string& StreamFileVersion() const {return plotFileVersion;}
  const amrex::Vector<int>& InsideNodes(int level, int box_id) const;
  int StreamIdxLo() const { return streamIdxLo; }
  int StreamIdxHi() const { return streamIdxHi; }
  static amrex::Box& invalidBox() {return ZBOX;}
  const amrex::FArrayBox* getFab (int level,
                                  int fabIdx,
                                  int compIdx);
  const amrex::DistributionMapping& DistributionMap (int level) const;
  void WriteFile(std::string&              outfile,
                 const amrex::Vector<amrex::MultiFab*>&  data,
                 int                       sComp,
                 int                       nComp,
                 const amrex::Vector<std::string>& names,
                 const std::string&        FileFormatName = StreamData::FileFormatName00()) const;

  void WriteFile(std::string&       outfile,
                 const amrex::Vector<int>&  comps,
                 const std::string& FileFormatName = StreamData::FileFormatName00()) const;

  int NodesPerElt() const {return nodes_per_element;}
  const amrex::Vector<int>& GlobalElements () const {return global_element_node_ids;}
  const amrex::Vector<int>& LocalElements (); // The list of globally-numbered nodes forming each element on this processor
  const amrex::Vector<int>& LocalElementIDs (); // The list of globally numbered elements local to this processor
  const std::vector<StreamData::MLloc>& GlobalNodeMap() const;
  const std::map<int,StreamData::MLloc>& LocalNodeMap();
  void PartitionElements (int nCompPass = 5);
  void PartitionElements (int nCompPass, const amrex::Vector<int>& comps);

  // Promote some AmrData functions
  using AmrData::NComp;
  using AmrData::GetFileName;
  using AmrData::StateNumber;
  using AmrData::GetGrids;
  using AmrData::FinestLevel;
  using AmrData::MinMax;
  using AmrData::ProbDomain;

  void FlushGrids(int componentIndex);

  const amrex::Vector<std::string> &StreamVarNames() const { return PlotVarNames(); }

  static std::string FileFormatName00() { return "Oddball-multilevel-connected-data-format";}
  static std::string FileFormatName10() { return "Oddball-multilevel-connected-data-format-1.0";}

  void SetStreamVerbose(int v) {stream_verbose = v;}
  int StreamVerbose() const {return stream_verbose;}

protected:
  void BuildGlobalNodeMap() const;

  // Actually does the read, note: can only flush the entire mf
  bool DefineFab(int level, int componentIndex, int fabIndex);

  std::string filename;
  amrex::Vector<int> global_element_node_ids;
  int num_global_elements;
  int nLines;
  int nodes_per_element;
  amrex::Vector<amrex::Vector<amrex::Vector<int> > > inside_nodes;
  int streamIdxLo, streamIdxHi;
  static amrex::Box ZBOX;
  bool is_version_0;
  std::string HeaderFileName;
  std::string ElementFileName;
  std::string ElementFileFormat;
  mutable std::vector<StreamData::MLloc> global_node_map;
  mutable std::map<int,StreamData::MLloc> local_node_map;

  // Used in the partitioning
  mutable std::map<int,std::map<int,const MLloc*> > remote_nodes;
  mutable std::map<int,std::map<int,const MLloc*> > tosend_nodes;
  mutable amrex::Vector<int> local_element_node_ids;
  mutable amrex::Vector<int> local_element_ids;
  mutable std::map<int,int> global_to_local_node_ids;
  mutable std::map<int,int> local_to_global_node_ids;
  mutable amrex::Vector<amrex::FArrayBox*> boundaryNodes;
  int stream_verbose;
};
#endif
