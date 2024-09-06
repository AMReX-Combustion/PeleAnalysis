#include <StreamData.H>
#include <AMReX_ParmParse.H>

using namespace amrex;
using std::cerr;
using std::endl;
using std::cout;
using std::ios;

#include <AMReX_Utility.H>

const std::string HeaderFileName_DEF("Header");
const std::string ElementFileName_DEF("Elements");
const std::string ElementFileFormat_DEF("ELEMENT_DATA_ASCII");


Box StreamData::ZBOX(IntVect::TheZeroVector(),IntVect::TheZeroVector());

StreamData::StreamData(const std::string& file)
    : AmrData(), filename(file)
{
    if (! (ReadInfo() ) )
        Abort(std::string(std::string("StreamData failed - bad filename: ")+file).c_str());
}

bool
StreamData::DefineFab(int level, int componentIndex, int fabIndex)
{   
    if( ! dataGridsDefined[level][componentIndex][fabIndex]) {
        BL_ASSERT(componentIndex < compIndexToVisMFMap.size());
        BL_ASSERT(componentIndex < compIndexToVisMFComponentMap.size());
        int whichVisMF(compIndexToVisMFMap[componentIndex]);
        int whichVisMFComponent(compIndexToVisMFComponentMap[componentIndex]);
        dataGrids[level][componentIndex]->setFab(
            fabIndex,
            std::unique_ptr<FArrayBox>(visMF[level][whichVisMF]->readFAB(fabIndex, whichVisMFComponent)));
        dataGridsDefined[level][componentIndex][fabIndex] = true;
    }
    return true;
}

const FArrayBox*
StreamData::getFab (int level,
                    int fabIdx,
                    int compIdx)
{
  if (fabIdx<0) {
      // Intends to access boundary/ghost data
      if (compIdx < boundaryNodes.size()) {
        if (!(boundaryNodes.size()>compIdx)
              || boundaryNodes[compIdx] == nullptr) {
              std::cout << "Boundary data not correctly partitioned for component " << compIdx << std::endl;
              Abort("Data not correctly partitioned");
          }
          return boundaryNodes[compIdx];
      }
      else {
          Abort("Data not correctly partitioned");
      }
  }
  DefineFab(level,compIdx,fabIdx);
  return &((*dataGrids[level][compIdx])[fabIdx]);
}

const BoxArray&
StreamData::boxArray(int level) const
{
    // use visMF[][0]:  all boxArrays are
    // guaranteed to be the same for each MultiFab
    return visMF[level][0]->boxArray();
}

const Vector<int>&
StreamData::InsideNodes(int level, int box_id) const
{
    BL_ASSERT(level<=finestLevel);
    BL_ASSERT(box_id<inside_nodes[level].size());
    return inside_nodes[level][box_id];
}

bool
StreamData::ReadInfo()
{
    VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);

    std::ifstream isPltIn;

    isPltIn.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());

    if(stream_verbose && ParallelDescriptor::IOProcessor()) {
            cout << "StreamData::opening file = " << filename << endl;
    }

    // Read header for path traces
    HeaderFileName = HeaderFileName_DEF;
    std::string HdrName = filename + "/" + HeaderFileName;

    isPltIn.open(HdrName.c_str(), ios::in);
    if(isPltIn.fail()) {
        if(ParallelDescriptor::IOProcessor()) {
            cerr << "Unable to open file: " << HdrName << endl;
        }
        return false;
    }

    isPltIn >> plotFileVersion;
    if (plotFileVersion == FileFormatName00()) {
        is_version_0 = true;
    } else if (plotFileVersion == FileFormatName10()) {
        is_version_0 = false;
    }
    else {
        Abort("Bad stream file version");
    }

    if(stream_verbose && ParallelDescriptor::IOProcessor()) {
        cout << "Stream file version:  " << plotFileVersion << endl;
    }

    int Nlev;
    isPltIn >> Nlev;
    finestLevel = Nlev - 1;

    isPltIn >> nComp;

    plotVars.resize(nComp);
    if(stream_verbose>1 && ParallelDescriptor::IOProcessor()) {
        cout << "Stream file variables (" << nComp << "): ";
    }
    std::set<std::string> unique_names;
    for (int i=0; i<nComp; ++i) {
        isPltIn >> plotVars[i];
        if (unique_names.count(plotVars[i])>0) {
          plotVars[i] += "_";
        }
        unique_names.insert(plotVars[i]);
        if(stream_verbose>1 && ParallelDescriptor::IOProcessor()) {
            cout << plotVars[i] << "[" << i << "] ";
        }
    }
    if(stream_verbose>1 && ParallelDescriptor::IOProcessor()) {
        cout << endl;
    }

    probDomain.resize(Nlev);

    if (is_version_0)
    {
        ElementFileName = ElementFileName_DEF;
        ElementFileFormat = ElementFileFormat_DEF;
    }
    else {

        isPltIn >> ElementFileName; 
        isPltIn >>  ElementFileFormat; 

        for (int i=0; i<BL_SPACEDIM; ++i)  isPltIn >> probLo[i];
        for (int i=0; i<BL_SPACEDIM; ++i)  isPltIn >> probHi[i];
        for (int lev=0; lev<Nlev; ++lev) {
            isPltIn >> probDomain[lev];
            // dont read ba
        }
    }

    ElementFileName = filename + "/" + ElementFileName;
    if(stream_verbose && ParallelDescriptor::IOProcessor()) {
      cout << "Reading streamline metadata" << endl;
    }
    dataGrids.resize(finestLevel + 1);
    dataGridsDefined.resize(finestLevel + 1);
    visMF.resize(finestLevel + 1);
    compIndexToVisMFMap.resize(nComp);
    compIndexToVisMFComponentMap.resize(nComp);

    MFInfo Fab_noallocate;
    Fab_noallocate.SetAlloc(false);

    bool found_valid_line_extent = false;
    for (int i = 0; i <= finestLevel; ++i) {

        // copied, but mostly bypassed code to account for multiple multifabs in a plot file
        int currentIndexComp(0);
        int currentVisMF(0);
        dataGrids[i].resize(nComp);
        dataGridsDefined[i].resize(nComp);
        
        while(currentIndexComp < nComp) {
            
          std::string ThisStreamFile = filename + "/";
            
          if (is_version_0) {
            ThisStreamFile += Concatenate("Level_",i,1) + "/Str";
          } else {
            std::string mfNameRelative;
            isPltIn >> mfNameRelative;
            ThisStreamFile += mfNameRelative;
          }

            visMF[i].resize(currentVisMF + 1);  // this preserves previous ones
            visMF[i][currentVisMF] = new VisMF(ThisStreamFile);
            int iComp(currentIndexComp);
            nGrow = visMF[i][currentVisMF]->nGrow();
            currentIndexComp += visMF[i][currentVisMF]->nComp();
            int currentVisMFComponent(0);
            for( ; iComp < currentIndexComp; ++iComp) {

                const BoxArray& ba = visMF[i][currentVisMF]->boxArray();
                // Find a valid box entry so that we can determine
                //    the lo,hi range of the pathlines
                for (int j=0; j<ba.size() && !found_valid_line_extent; ++j) {
                    const Box& bx = ba[j];
                    if (bx != ZBOX) {
                        streamIdxLo = bx.smallEnd(1);
                        streamIdxHi = bx.bigEnd(1);
                        found_valid_line_extent = true;
                    }
                }

                // make single component multifabs
                // defer reading the MultiFab data
                dataGrids[i][iComp] = new MultiFab(ba, DistributionMapping(ba), 1,
                                                   visMF[i][currentVisMF]->nGrow(),Fab_noallocate);
                dataGridsDefined[i][iComp].resize(visMF[i][currentVisMF]->size(),
                                                  false);
                compIndexToVisMFMap[iComp] = currentVisMF;
                compIndexToVisMFComponentMap[iComp] = currentVisMFComponent;
                ++currentVisMFComponent;
            }
            ++currentVisMF;
        }  // end while        
    }  // end for(i...finestLevel)

    isPltIn.close();

    // Read elements
    // FIXME: Should we read in serial and broadcast these??
    isPltIn.open(ElementFileName.c_str());
    if(isPltIn.fail()) {
        if(ParallelDescriptor::IOProcessor()) {
            cerr << "Unable to open file: " << ElementFileName << endl;
        }
        return false;
    }

    isPltIn >> num_global_elements;
    isPltIn >> nodes_per_element;
    if(stream_verbose && ParallelDescriptor::IOProcessor()) {
            cout << "Reading elements (num_global_elements,nodes_per_element): ("
                 << num_global_elements << ", " << nodes_per_element << ")" << endl;
    }
    global_element_node_ids.resize(num_global_elements*nodes_per_element);
    for (int i=0; i<num_global_elements*nodes_per_element; ++i) {
        isPltIn >> global_element_node_ids[i];
        global_element_node_ids[i]--; // Historical crud, these are written as 1-based
    }

    // Read element distribution after we know number of boxes/level
    if(stream_verbose && ParallelDescriptor::IOProcessor()) {
            cout << "Populating map of node distribution... " << endl;
    }
    inside_nodes.resize(finestLevel + 1);

    nLines = 0;
    for (int lev=0; lev<=finestLevel; ++lev)
    {
        inside_nodes[lev].resize(visMF[lev][0]->size());
        
        // Get number of nonempty boxes
        int num_non_zero;
        isPltIn >> num_non_zero;
        
        for (int j=0; j<num_non_zero; ++j)
        {
            // Get index of non-empty box, and number of entries
            int box_id, num_ids;
            isPltIn >> box_id >> num_ids;
            
            inside_nodes[lev][box_id].resize(num_ids);
            nLines += num_ids;
            
            for (int k=0; k<num_ids; ++k)
                isPltIn >> inside_nodes[lev][box_id][k];
        }
    }
    isPltIn.close();

    if(stream_verbose && ParallelDescriptor::IOProcessor()) {
        cout << "Finished reading headers for " << nLines << " stream lines and " << num_global_elements << " elements" << endl;
    }

    return true;

}

void
StreamData::WriteFile(std::string&       outfile,
                      const Vector<int>&  comps,
                      const std::string& FileFormatName) const
{
  Abort("Not implemented yet");
}

void
StreamData::WriteFile(std::string&               outfile,
                      const Vector<MultiFab*>&   data,
                      int                        sComp,
                      int                        _nComp,
                      const Vector<std::string>& names,
                      const std::string&         FileFormatName) const
{
    BL_ASSERT(names.size() == _nComp);

    const std::string LevelDirName("Level");
    const std::string StreamDataFileName("Str");

    if (ParallelDescriptor::IOProcessor())
    {
        if (!UtilCreateDirectory(outfile, 0755))
            CreateDirectoryFailed(outfile);

        // Write Header
        const std::string FullHeaderFileName = outfile + "/" + HeaderFileName_DEF;
        std::ofstream ofh;
        ofh.open(FullHeaderFileName.c_str());
        ofh << FileFormatName << '\n';
        ofh << data.size() << '\n'; // number of levels
        ofh << names.size() << '\n'; // number of variables per node
        for (int i=0; i<names.size(); ++i) {
            ofh << names[i] << '\n'; 
        }

        if (!is_version_0) {
            ofh << ElementFileName << '\n'; 
            ofh << ElementFileFormat << '\n'; 
            
            for (int i=0; i<BL_SPACEDIM; ++i)  ofh << probLo[i] << " ";
            ofh << '\n';
            for (int i=0; i<BL_SPACEDIM; ++i)  ofh << probHi[i] << " ";
            ofh << '\n';
            for (int lev=0; lev<data.size(); ++lev) {
                ofh << ProbDomain()[lev] << '\n';
                ofh << boxArray(lev) << '\n';
            }
        }
        ofh.close();

        // Write elements
        const std::string FullElementFileName = outfile + "/" + ElementFileName;
        int my_nodes_per_element = global_element_node_ids.size() / num_global_elements;
        std::ofstream ofe;
        ofe.open(FullElementFileName.c_str());
        ofe << num_global_elements << '\n';
        ofe << my_nodes_per_element << '\n';
        for (int i=0; i < global_element_node_ids.size(); ++i)
            ofe << global_element_node_ids[i] + 1 << " "; // Historical crud, these are written as 1-based
        ofe << '\n';

        // Write element distribution
        for (int i=0; i<inside_nodes.size(); ++i)
        {
            const Vector<Vector<int> >& inside_nodes_i = inside_nodes[i];

            // Count the non-zero-length entries
            int num_non_zero = 0;
            for (int j=0; j<inside_nodes_i.size(); ++j)
                if (inside_nodes_i[j].size() > 0)
                    num_non_zero++;

            ofe << num_non_zero << '\n';
            for (int j=0; j<inside_nodes_i.size(); ++j)
            {
                const Vector<int>& inside_nodes_i_j = inside_nodes[i][j];
                if (inside_nodes_i[j].size() > 0)
                {
                    ofe << j << " " << inside_nodes_i_j.size();
                    for (int k=0; k<inside_nodes_i_j.size(); ++k)
                        ofe << " " << inside_nodes_i_j[k];
                    ofe << '\n';
                }
            }
        }
        ofe.close();
    }
    //
    // Force other processors to wait till directory is built and Header is written.
    //
    ParallelDescriptor::Barrier();

    // Write data
    for (int lev=0; lev<data.size(); ++lev)
    {
        char buf[64];
        sprintf(buf, "Level_%d", lev);
        std::string DataFile = outfile + "/" + std::string(buf);
        
        if (ParallelDescriptor::IOProcessor())
            if (!UtilCreateDirectory(DataFile, 0755))
                CreateDirectoryFailed(DataFile);
        //
        // Force other processors to wait till directory is built.
        //
        ParallelDescriptor::Barrier();
        
        std::string FabDataFile = DataFile + "/Str";

        if (names.size() == data[lev]->size())
        {
            VisMF::Write(*data[lev],FabDataFile);
        }
        else
        {
          const DistributionMapping dml(data[lev]->boxArray());
          MultiFab tmp(data[lev]->boxArray(),dml,names.size(),data[lev]->nGrow());
          MultiFab::Copy(tmp,*data[lev],sComp,0,names.size(),data[lev]->nGrow());
          VisMF::Write(tmp,FabDataFile);
        }
    }
}

const DistributionMapping&
StreamData::DistributionMap (int level) const
{
    // All components will have the same distribution map, by assumption
    BL_ASSERT(level<=finestLevel);
    return dataGrids[level][0]->DistributionMap();
}

const std::vector<StreamData::MLloc>&
StreamData::GlobalNodeMap() const
{
    if (global_node_map.size()==0) {
        BuildGlobalNodeMap();
    }
    return global_node_map;
}

void
StreamData::FlushGrids(int componentIndex)
{
    AmrData::FlushGrids(componentIndex);
    if (boundaryNodes.size() > componentIndex && boundaryNodes[componentIndex]!=nullptr) {
        delete boundaryNodes[componentIndex];
        boundaryNodes[componentIndex] = new FArrayBox();
    }
}

void
StreamData::BuildGlobalNodeMap() const
{
    if (global_node_map.size() == 0) {
        global_node_map.resize(NLines());
        int myfinestLevel = FinestLevel();
        
        for (int iLev=0; iLev<=myfinestLevel; ++iLev)
        {
            int nBoxes = boxArray(iLev).size();
            for (int iBox=0; iBox<nBoxes; ++iBox)
            {
                const Vector<int>& nodes_this_box = InsideNodes(iLev,iBox);
                int num_nodes = nodes_this_box.size();
                for (int iNode=0; iNode<num_nodes; ++iNode)
                {
                    // Remember that inside_nodes is 1-based
                    global_node_map[nodes_this_box[iNode] - 1] = MLloc(iLev,iBox,iNode); 
                }
            }
        }

        // Get distribution on each level
        Vector<DistributionMapping> dm(myfinestLevel+1);
        for (int iLev = 0; iLev <= myfinestLevel; ++iLev) {
            dm[iLev] = DistributionMap(iLev);
        }
  
        int mynodes_per_element = NodesPerElt();
        int num_elements_global = NElts();
  
        int my_proc = ParallelDescriptor::MyProc();
  
        Vector<const MLloc*> n(mynodes_per_element);
        for (int i=0; i<num_elements_global; ++i)
        {
            const int offset = i * mynodes_per_element;
    
            for (int j=0; j<mynodes_per_element; ++j) {
                int global_node_id = global_element_node_ids[offset + j];
                BL_ASSERT(global_node_map.size() > global_node_id);
                n[j] = &(global_node_map[global_node_id]);
            }
    
            BL_ASSERT(n[0]->amr_lev <= myfinestLevel);
            BL_ASSERT(dm.size()>n[0]->amr_lev && dm[n[0]->amr_lev].size()>n[0]->box_idx);
            int element_owner = dm[n[0]->amr_lev][n[0]->box_idx];

            bool is_shared = false;
            for (int j=1; j<mynodes_per_element; ++j) {
      
                BL_ASSERT(n[j]->amr_lev <= myfinestLevel);
                BL_ASSERT(dm.size()>n[j]->amr_lev && dm[n[j]->amr_lev].size()>n[j]->box_idx);
                int node_owner = dm[n[j]->amr_lev][n[j]->box_idx];

                if (element_owner != node_owner) {
                    element_owner = rand() > 0.5  ?  element_owner :  node_owner;
                    is_shared = true;
                }
            }
    
            if (is_shared) {
                for (int j=0; j<mynodes_per_element; ++j) {
                    int node_owner = dm[n[j]->amr_lev][n[j]->box_idx];
                    int global_node_id = global_element_node_ids[offset + j];
      
                    if ( (element_owner==my_proc) && (element_owner!=node_owner)) {
                        remote_nodes[node_owner][global_node_id] = n[j];
                    }
        
                    if ( (node_owner==my_proc) && (element_owner!=my_proc) ) {
                        tosend_nodes[element_owner][global_node_id] = n[j];
                    }
                }
            }

            if (element_owner == my_proc) {
                local_element_ids.push_back(i);
            }
        }

        // Build map to new nodes
        int new_local_node_cnt = 0;
        for (std::map<int,std::map<int,const MLloc*> >::const_iterator it=remote_nodes.begin(),
                 End=remote_nodes.end(); it!=End ; ++it)
        {
            const std::map<int,const MLloc*>& nodes_recd = it->second;
            for (std::map<int,const MLloc*>::const_iterator it1=nodes_recd.begin(),
                     End1=nodes_recd.end(); it1!=End1; ++it1)
            {
                int global_node_id = it1->first;      
                global_to_local_node_ids[global_node_id] = new_local_node_cnt; // Assign new node number
                local_to_global_node_ids[new_local_node_cnt] = global_node_id; // Map new number back to global number
                new_local_node_cnt++;
            }
        }

#if 0
        for (int n=0; n<nprocs; ++n) {
            if (stream_verbose && my_proc==n) {
                std::cout << "proc " << n << " will build " << new_local_node_cnt << " (ghost) nodes." << std::endl;
            }
            ParallelDescriptor::Barrier();
        }
        ParallelDescriptor::Barrier();
#endif

        // Build new set of local elements 
        int num_local_elements = local_element_ids.size();
        local_element_node_ids.resize(mynodes_per_element * num_local_elements);
        int* ptr = local_element_node_ids.dataPtr();
        
        local_node_map.clear();
        for (int i=0, End=local_element_ids.size(); i<End; ++i)
        {
            int global_element_id = local_element_ids[i];
            const int offset_global = global_element_id * mynodes_per_element;
            for (int j=0; j<mynodes_per_element; ++j)
            {
                int global_node_id = global_element_node_ids[offset_global + j];
                int amr_lev, box_idx, pt_idx;
                
                // Local nodes are "ghost" nodes if they are in this map
                //  ghost nodes are stored in the boundaryNodes fab
                if (global_to_local_node_ids.count(global_node_id))
                {
                    amr_lev = -1;
                    box_idx = -1;
                    pt_idx = global_to_local_node_ids[global_node_id];
                }
                else
                {
                    // Otherwise, the node is in the usual data structure
                    BL_ASSERT(global_node_map.size() > global_node_id);
                    const MLloc& node = global_node_map[global_node_id];
                    amr_lev = node.amr_lev;
                    box_idx = node.box_idx;
                    pt_idx = node.pt_idx;
                    
                    BL_ASSERT(DistributionMap(amr_lev)[box_idx] == my_proc);
                }
                
                *ptr++ = global_node_id;
                
                // populate local node map
                if (local_node_map.count(global_node_id)==0) {
                    local_node_map[global_node_id] = MLloc(amr_lev,box_idx,pt_idx);
                }
            }
        }        

#if 0
        for (int n=0; n<nprocs; ++n) {
            if (stream_verbose && my_proc==n) {
                std::cout << "proc " << n << " now has " << local_node_map.size() << " nodes." << std::endl;
            }
            ParallelDescriptor::Barrier();
        }
        ParallelDescriptor::Barrier();
#endif
    }
}

const Vector<int>&
StreamData::LocalElements ()
{
  return local_element_node_ids;
}

const Vector<int>&
StreamData::LocalElementIDs ()
{
  return local_element_ids;
}

const std::map<int,StreamData::MLloc>&
StreamData::LocalNodeMap()
{
  return local_node_map;
}

typedef StreamData::MLloc MLloc;
static void DoIt(const Vector<int>& comps,
                 Vector<int>& sendcnts,
                 Vector<int>& recvcnts,
                 Vector<int>& sdispls,
                 Vector<int>& rdispls,
                 Vector<Real>& senddata,
                 Vector<Real>& recvdata,
                 const std::map<int,std::map<int,const MLloc*> >& remote_nodes,
                 const std::map<int,std::map<int,const MLloc*> >& tosend_nodes,
                 Vector<FArrayBox*>& dest,
                 StreamData& stream,
                 std::map<int,int>& global_to_local_node_ids)
{
  int num_procs = ParallelDescriptor::NProcs();

  // Build data structures for communicating node info
  sendcnts.resize(num_procs,0);
  recvcnts.resize(num_procs,0);
  sdispls.resize(num_procs,0);
  rdispls.resize(num_procs,0);

  int points_on_line = stream.StreamIdxHi()-stream.StreamIdxLo()+1;
  int nComp = comps.size();
  int chunk_size = points_on_line * nComp;

  int tot_num_to_send = 0;
  for (std::map<int,std::map<int,const MLloc*> >::const_iterator it=tosend_nodes.begin(), 
         End=tosend_nodes.end(); it!=End; ++it) 
  {
    sdispls[it->first] = tot_num_to_send;
    sendcnts[it->first] = it->second.size() * chunk_size;
    tot_num_to_send += sendcnts[it->first];
  }
    
  int num_recv_nodes = 0;
  int tot_num_to_recv = 0;
  for (std::map<int,std::map<int,const MLloc*> >::const_iterator it=remote_nodes.begin(), 
         End=remote_nodes.end(); it!=End ; ++it) 
  {
    rdispls[it->first] = tot_num_to_recv;
    num_recv_nodes += it->second.size();
    recvcnts[it->first] = it->second.size() * chunk_size;
    tot_num_to_recv += recvcnts[it->first];
  }
    
  // Load data to send, then communicate
  senddata.resize(tot_num_to_send);
  recvdata.resize(tot_num_to_recv);

  StreamData& Cstream = const_cast<StreamData&>(stream); // A bit ugly, but for the greater good...
    
  Vector<int> offsets = sdispls;
  for (std::map<int,std::map<int,const MLloc*> >::const_iterator it=tosend_nodes.begin(), 
         End=tosend_nodes.end(); it!=End; ++it) 
  { 
    int to_proc = it->first;
    const std::map<int,const MLloc*>& proc_nodes = it->second;
      
    int icnt = 0;
    for (std::map<int,const MLloc*>::const_iterator it1=proc_nodes.begin(), End1=proc_nodes.end(); 
         it1!=End1; ++it1, ++icnt) 
    {
      const MLloc* node = it1->second;
      int local_node_id = node->pt_idx;
      for (int n=0; n<comps.size(); ++n) {
        const FArrayBox& stream_fab = *(Cstream.getFab(node->amr_lev,node->box_idx,comps[n]));
        for (int j=stream.StreamIdxLo(); j<=stream.StreamIdxHi(); ++j) {
          IntVect iv(AMREX_D_DECL(local_node_id,j,0));
          senddata[offsets[to_proc]++] = stream_fab(iv,0);
        }
      }
    }
  }

#if BL_USE_MPI
  BL_MPI_REQUIRE( MPI_Alltoallv(tot_num_to_send == 0 ? 0 : senddata.dataPtr(),
                                sendcnts.dataPtr(),
                                sdispls.dataPtr(),
                                ParallelDescriptor::Mpi_typemap<Real>::type(),
                                tot_num_to_recv == 0 ? 0 : recvdata.dataPtr(),
                                recvcnts.dataPtr(),
                                rdispls.dataPtr(),
                                ParallelDescriptor::Mpi_typemap<Real>::type(),
                                ParallelDescriptor::Communicator()) );
#endif    

  if (num_recv_nodes>0) {
      
    Box d_box(IntVect(AMREX_D_DECL(0,stream.StreamIdxLo(),0)),
              IntVect(AMREX_D_DECL(num_recv_nodes-1,stream.StreamIdxHi(),0)));
      
    for (int n=0; n<comps.size(); ++n) {
        dest[comps[n]]->resize(d_box,1);
    }

    for (std::map<int,std::map<int,const MLloc*> >::const_iterator it=remote_nodes.begin(), 
           End=remote_nodes.end(); it!=End ; ++it)
    {
      int from_proc = it->first;

      if (rdispls[from_proc] + it->second.size() >= tot_num_to_recv) {
        Abort("Bad comm data");
      }

      Real* dest_loc = &recvdata[rdispls[from_proc]];
        
      for (std::map<int,const MLloc*>::const_iterator it1=it->second.begin(), 
             End1=it->second.end(); it1!=End1; ++it1)
      {
        int local_node_id = global_to_local_node_ids[it1->first];
        for (int n=0; n<comps.size(); ++n) {
          for (int j=d_box.smallEnd()[1]; j<=d_box.bigEnd()[1]; ++j) {
            IntVect idx(AMREX_D_DECL(local_node_id,j,0));
            (*dest[comps[n]])(idx,0) = *dest_loc++;
          }
        }
      }        
    }
  }
}



void
StreamData::PartitionElements (int nCompPass)
{
  Vector<int> comps(NComp());
  for (int n=0; n<comps.size(); ++n) {
    comps[n] = n;
  }
  PartitionElements(nCompPass, comps);
}

void
StreamData::PartitionElements (int nCompPass, const Vector<int>& comps)
{
  if (stream_verbose && ParallelDescriptor::IOProcessor() ) {
    cout << "Partitioning elements..." << std::endl;
  }

  for (int i=0; i<comps.size(); ++i) {
    if (comps[i]<0 || comps[i]>NComp()) {
        cout << "comp: " << comps[i] << endl;
      Error("desired component not on streamline");
    }
  }

  int num_comp = comps.size();
  
  // Build a structure that for each node points to where in the streamlines to get the data
  BuildGlobalNodeMap();
  
  Vector<int> sendcnts;
  Vector<int> recvcnts;
  Vector<int> sdispls;
  Vector<int> rdispls;
  Vector<Real> senddata;
  Vector<Real> recvdata;
  
  int num_tot = NComp();
  if (boundaryNodes.size() != num_tot) {
      boundaryNodes.resize(num_tot,nullptr);
      for (int i=0; i<num_tot; ++i) {
          boundaryNodes[i] = new FArrayBox();
      }
  }
  
  bool isio = ParallelDescriptor::IOProcessor();  
  int src_comp = 0;
  while (src_comp < num_comp) 
  {
    int nCompPassLoc = std::min(src_comp + nCompPass, num_comp) - src_comp;
    Vector<int> loc_comps(nCompPassLoc);
    for (int n=0; n<loc_comps.size(); ++n) {
      loc_comps[n] = comps[src_comp + n];
    }
    
    if (stream_verbose && isio) {
      std::cout << "Reading " << loc_comps.size() << " components of stream data ";
      if (1 || stream_verbose>1) {
        std::cout << ": ";
        for (int j=0; j<loc_comps.size(); ++j) {
          std::cout << PlotVarNames()[loc_comps[j]] << " ";
        }
      }
      std::cout << std::endl;
    }
    DoIt(loc_comps, sendcnts, recvcnts, sdispls, rdispls, senddata, recvdata,
         remote_nodes, tosend_nodes, boundaryNodes, *this, global_to_local_node_ids);
    
    src_comp += nCompPassLoc;
  }
  if (stream_verbose && ParallelDescriptor::IOProcessor()) {
    cout << "Finished partitioning" << endl;
  }
}

