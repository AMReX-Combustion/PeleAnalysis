#include <string>
#include <iostream>
#include <set>

#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include <AMReX_DataServices.H>
#include <AMReX_BCRec.H>
#include <AMReX_Interpolater.H>

#include <AMReX_BLFort.H>

using namespace amrex;

extern "C"
{
    void interpstream(const Real* loc, const int* loc_lo, const int* loc_hi, const int* nl,
                      const Real* fab, const int* fab_lo, const int* fab_hi, const int* np,
                      Real* strm, const int* strm_lo, const int* strm_hi,
                      const Real* dx, const Real* plo);

    void set_distance(const Real* loc, const int* loc_lo, const int* loc_hi,
                      Real* res, const int* res_lo, const int* res_hi);

}

void
read_ml_streamline_data(const std::string& outfile,
                        Vector<MultiFab*>& data,
                        vector<string>&    names,
                        Vector<int>&       faceData,
                        int&               nElts,
                        Vector<Vector<Vector<int> > >& inside_nodes);

void
write_ml_streamline_data(const std::string&       outfile,
                         const Vector<MultiFab*>& data,
                         const vector<string>&    names,
                         const Vector<int>&       faceData,
                         int                      nElts,
                         const Vector<Vector<Vector<int> > >& inside_nodes);

void
dump_ml_streamline_data(const std::string&       outfile,
                        const Vector<MultiFab*>& data,
                        int                      sComp,
                        const vector<string>&    names,
                        const Vector<int>&       faceData,
                        int                      nElts,
                        const Vector<Vector<Vector<int> > >& inside_nodes);

void
sample_pathlines(Vector<MultiFab*>&       sampledata,
                 int                      sComp,
                 int                      dComp,
                 int                      nComp,
                 const Vector<int>&       comps,
                 AmrData&                 amrData,
                 int                      nGrow,
                 const Vector<MultiFab*>& pathlines,
                 const Vector<Vector<Vector<int> > >& inside_nodes,
                 const Vector<Geometry*>& geoms);

void
set_sample_distance(Vector<MultiFab*>&       sampledata,
                    int                      dComp,
                    const Vector<MultiFab*>& pathlines);

void
set_sample_location(Vector<MultiFab*>&       sampledata,
                    int                      dComp,
                    const Vector<MultiFab*>& pathlines);

Real strt_time, strt_io, io_time;

int
main (int   argc,
      char* argv[])
{
    amrex::Initialize(argc,argv);
    {
    strt_time = ParallelDescriptor::second();
    io_time = 0;

    ParmParse pp;

    std::string plotfile; pp.get("plotfile",plotfile);
    DataServices::SetBatchMode();
    Amrvis::FileType fileType(Amrvis::NEWPLT);

    strt_io = ParallelDescriptor::second();
    DataServices dataServices(plotfile, fileType);
    if( ! dataServices.AmrDataOk()) {
        DataServices::Dispatch(DataServices::ExitRequest, NULL);
        // ^^^ this calls ParallelDescriptor::EndParallel() and exit()
    }
    AmrData& amrData = dataServices.AmrDataRef();
    io_time += ParallelDescriptor::second() - strt_io;

    int finestLevel = amrData.FinestLevel(); pp.query("finestLevel",finestLevel);
    int Nlev = finestLevel + 1;
    int nGrow = 4; pp.query("nGrow",nGrow); // Set this depending on how long the lines are going to be

    // Set components to read from plotfile
    Vector<int> comps;
    if (int nc = pp.countval("comps"))
    {
        comps.resize(nc);
        pp.getarr("comps",comps,0,nc);
    }
    else
    {
        int sCompPlt = 0; pp.query("sComp",sCompPlt);
        int nCompPlt = amrData.NComp(); pp.query("nComp",nCompPlt);
        BL_ASSERT(sCompPlt+nCompPlt <= amrData.NComp());
        comps.resize(nCompPlt);
        for (int i=0; i<nCompPlt; ++i)
            comps[i] = sCompPlt + i;
    }
    const int nCompPlt = comps.size();

    std::string pathFile; pp.get("pathFile",pathFile);
    vector<string> infileTokens = Tokenize(pathFile,".");

    
    string outfile = infileTokens[0];
    for (int i=1; i<infileTokens.size()-1; ++i)
        outfile += string(".") + infileTokens[i];
    outfile += "_sampled_from_" + plotfile;
    pp.query("outfile",outfile);

    // Read pathlines
    Vector<MultiFab*> pathlines(Nlev);
    Vector<string> strVarNames;
    Vector<int> faceData;
    int nElts;
    Vector<Vector<Vector<int> > > inside_nodes;
    Print() << "reading streamline data" << std::endl;
    strt_io = ParallelDescriptor::second();
    read_ml_streamline_data(pathFile,pathlines,strVarNames,faceData,nElts,inside_nodes);
    Print() << "done reading streamline data" << std::endl;
    io_time += ParallelDescriptor::second() - strt_io;

    // Make space for result
    int nCompSamp = nCompPlt + 1 + BL_SPACEDIM; // Add for distance function, and for locations
    Vector<MultiFab*> sampledata(Nlev);
    for (int lev=0; lev<Nlev; ++lev) {
        const BoxArray& pba = pathlines[lev]->boxArray();
        const DistributionMapping& dmpba = pathlines[lev]->DistributionMap();
        sampledata[lev] = new MultiFab(pba,dmpba,nCompSamp,0);
        sampledata[lev]->setVal(0.);
    }

    // Sample plotfile on paths (do all components at once if we've got lots of memory, nCompsPerPass=-1)
    int nCompsPerPass = -1; pp.query("nCompsPerPass", nCompsPerPass); 
    if (nCompsPerPass<0) 
    {
        nCompsPerPass = nCompPlt;
    }

    RealBox rb(&(amrData.ProbLo()[0]),
               &(amrData.ProbHi()[0]));
    int coord = 0;
    Vector<int> is_per(BL_SPACEDIM,1);
    pp.queryarr("is_per",is_per,0,BL_SPACEDIM);
    Print() << "Periodicity assumed for this case: ";
    for (int i=0; i<BL_SPACEDIM; ++i) {
        Print() << is_per[i] << " ";
    }
    Print() << std::endl;
    Vector<Geometry*> geoms(Nlev);
    for (int lev=0; lev<Nlev; ++lev) {
        geoms[lev] = new Geometry(amrData.ProbDomain()[lev],&rb,coord,&(is_per[0]));
    }

    int i=0;
    while (i<nCompPlt)
    {
        int nCompWork = std::min(nCompsPerPass, nCompPlt - i);
        
        // Sample the data in chunks
        sample_pathlines(sampledata,i,i+BL_SPACEDIM+1,nCompWork,comps,amrData,nGrow,pathlines,inside_nodes,geoms);
        
        // Move to next chunk
        i += nCompsPerPass;
    }

    Print() << "done sampling data" << std::endl;

    // Load locations at bottom
    set_sample_location(sampledata,0,pathlines);

    // Load distance from anchor point above locations
    set_sample_distance(sampledata,BL_SPACEDIM,pathlines);

    // Load locations onto paths
    // Write results
    Vector<string> inVarNames(nCompSamp);
    const Vector<string>& plotVarNames = amrData.PlotVarNames();
    inVarNames[0] = "X";
    inVarNames[1] = "Y";
    inVarNames[2] = "Z";
    inVarNames[BL_SPACEDIM] = "distance_from_seed";
    for (int n=0; n<nCompPlt; ++n)
        inVarNames[n+BL_SPACEDIM+1] = plotVarNames[comps[n]];


    if (pp.countval("streamSampleFile") > 0)
    {
      ParallelDescriptor::Barrier();
      if (ParallelDescriptor::IOProcessor())
        std::cerr << "Writing the streamline data " << std::endl;

      // Get output filename
      std::string streamFile; pp.get("streamSampleFile",streamFile);
      if (!streamFile.empty() && streamFile[streamFile.length()-1] != '/')
        streamFile += '/';
    
      strt_io = ParallelDescriptor::second();
      write_ml_streamline_data(streamFile,sampledata,inVarNames,faceData,nElts,inside_nodes);
      io_time += ParallelDescriptor::second() - strt_io;

      if (ParallelDescriptor::IOProcessor())
        std::cerr << "...done writing the streamline data " << std::endl;
    }
    else
    {
      std::string outFile;
      if (pp.countval("outFile")==0) Abort("Must specify streamSampleFile or outFile");
      pp.get("outFile",outFile);
      dump_ml_streamline_data(outFile,sampledata,0,inVarNames,faceData,nElts,inside_nodes);
    }
    


    
    //strt_io = ParallelDescriptor::second();
    //write_ml_streamline_data(outfile, sampledata, inVarNames, faceData, nElts, inside_nodes);
    //io_time += ParallelDescriptor::second() - strt_io;


    long nInsideNodes = 0;
    for (int k=0; k<inside_nodes.size(); ++k) {
        for (int j=0; j<inside_nodes[k].size(); ++j) {
            nInsideNodes += inside_nodes[k][j].size();
        }
    }

    amrex::Finalize();
    }
    return 0;
}

void
FillCFgrowCells(AmrData& amrd, const BoxArray& fine_ba,Vector<Geometry*>& geoms, int lev,
                MultiFab& crseSrc, int sComp, int nComp, int nGrow, MultiFab& growScr)
{
    int refRatio = amrd.RefRatio()[lev-1];
    //const BoxArray cfBoxesCRSE = CFBoxes(BoxArray(fine_ba).coarsen(refRatio),nGrow,*geoms[lev-1]);
    const BoxArray cfBoxesCRSE = BoxArray(GetBndryCells(BoxArray(fine_ba).coarsen(refRatio),nGrow));
    const DistributionMapping dmScrCRSE(cfBoxesCRSE);
    MultiFab growScrCRSE(cfBoxesCRSE,dmScrCRSE,nComp,0);

    crseSrc.FillBoundary(sComp,nComp,geoms[lev-1]->periodicity());
    
    BoxArray gcba = BoxArray(crseSrc.boxArray()).grow(nGrow);
    MultiFab scrCRSE(gcba,crseSrc.DistributionMap(),nComp,0); // Build growLess mf to get f-f and p-f grow data through parallel copy
    for (MFIter mfi(scrCRSE); mfi.isValid(); ++mfi)
        scrCRSE[mfi].copy(crseSrc[mfi],sComp,0,nComp);

    // Now we have good grow data in scrCRSE, copy to fill interp source data
    growScrCRSE.copy(scrCRSE,0,0,nComp);
    
    // Do interp from coarse data to fine
    const BoxArray cfBoxes = BoxArray(cfBoxesCRSE).refine(refRatio); // Get matching fine boxes for interp
    MultiFab growScrFC(cfBoxes,dmScrCRSE,nComp,0); // Container on matching fine boxes
    Vector<BCRec> bc(1); // unused for pc_interp...
    int idummy1=0, idummy2=0;
    for (MFIter mfi(growScrFC);mfi.isValid();++mfi) {
        pc_interp.interp(growScrCRSE[mfi],0,growScrFC[mfi],0,1,growScrFC[mfi].box(),
                         refRatio*IntVect::TheUnitVector(),*geoms[lev-1],*geoms[lev],bc,idummy1,idummy2,RunOn::Cpu);
    }    
    // Finally, build correct c-f boxes and copy-on-intersect from interp data
    const BoxArray cfBoxesFINE = BoxArray(GetBndryCells(fine_ba,nGrow));
    DistributionMapping dmBoxesFINE(cfBoxesFINE);
    growScr.define(cfBoxesFINE,dmBoxesFINE,nComp,0);
    growScr.copy(growScrFC); // Get c-f data into just-right size container
}

void
write_ml_streamline_data(const std::string&      outfile,
                         const Vector<MultiFab*>& data,
                         const vector<string>&   names,
                         const Vector<int>&       faceData,
                         int                     nElts,
                         const Vector<Vector<Vector<int> > >& inside_nodes)
{
    if (ParallelDescriptor::IOProcessor())
    {
        if (!UtilCreateDirectory(outfile, 0755))
            CreateDirectoryFailed(outfile);

        // Write Header
        const string oHeader = outfile + string("/Header");
        std::ofstream ofh;
        ofh.open(oHeader.c_str());
        ofh << "Oddball-multilevel-connected-data-format" << '\n';
        ofh << data.size() << '\n';
        ofh << names.size() << '\n';
        for (int i=0; i<names.size(); ++i)
            ofh << names[i] << '\n';
        ofh.close();

        // Write elements
        int nodesPerElt = faceData.size() / nElts;
        const string oElements = outfile + string("/Elements");
        std::ofstream ofe;
        ofe.open(oElements.c_str());
        ofe << nElts << '\n';
        ofe << nodesPerElt << '\n';
        for (int i=0; i<faceData.size(); ++i)
            ofe << faceData[i] << " ";
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
        sprintf(buf, "/Level_%d", lev);
        std::string DataFile = outfile + std::string(buf);
        
        if (ParallelDescriptor::IOProcessor())
            if (!UtilCreateDirectory(DataFile, 0755))
                CreateDirectoryFailed(DataFile);
        //
        // Force other processors to wait till directory is built.
        //
        ParallelDescriptor::Barrier();
        
        std::string FabDataFile = DataFile + "/Str";

        VisMF::Write(*data[lev],FabDataFile);
    }
}

void
dump_ml_streamline_data(const std::string&       outFile,
                        const Vector<MultiFab*>& data,
                        int                      sComp,
                        const vector<string>&    names,
                        const Vector<int>&       faceData,
                        int                      nElts,
                        const Vector<Vector<Vector<int> > >& inside_nodes)
{
  // Create a folder and have each processor write their own data, one file per streamline
  auto myProc = ParallelDescriptor::MyProc();
  auto nProcs = ParallelDescriptor::NProcs();
  int cnt = 0;
  
  if (!amrex::UtilCreateDirectory(outFile, 0755))
    amrex::CreateDirectoryFailed(outFile);
  ParallelDescriptor::Barrier();
  
  const Box nullBox(IntVect::TheZeroVector(),IntVect::TheZeroVector());
  for (int lev=0; lev<data.size(); ++lev)
  {
    const auto& ba = data[lev]->boxArray();
    const auto& dm = data[lev]->DistributionMap();
    int nComp = data[lev]->nComp();
    for (int j=0; j<ba.size(); ++j)
    {
      if (dm[j] == myProc)
      {
        const auto& b = ba[j];
        if (b!=nullBox)
        {
          std::string fileName = outFile + "/str_";
          fileName = Concatenate(fileName,myProc) + "_";
          fileName = Concatenate(fileName,cnt++);
          std::ofstream ofs(fileName.c_str());

          for (int i=0; i<names.size(); ++i) {
            ofs << names[i] << " ";
          }
          ofs << '\n';
          
          for (auto iv=b.smallEnd(); iv<=b.bigEnd(); b.next(iv))
          {
            for (int i=sComp; i<nComp; ++i)
            {
              ofs << (*data[lev])[j](iv,i) << " ";
            }
            ofs << '\n';
          }
          ofs.close();
        }
      }
      ParallelDescriptor::Barrier();
    }
  }
}

void
read_ml_streamline_data(const std::string& outfile,
                        Vector<MultiFab*>&  data,
                        vector<string>&    names,
                        Vector<int>&        faceData,
                        int&               nElts,
                        Vector<Vector<Vector<int> > >& inside_nodes)
{
    // Read header for path traces
    const string Header = outfile + string("/Header");
    std::ifstream ifh;
    ifh.open(Header.c_str());
    std::string typeLabel; ifh >> typeLabel;
    int NlevPath; ifh >> NlevPath;
    int nCompPath; ifh >> nCompPath;
    names.resize(nCompPath);
    for (int i=0; i<nCompPath; ++i)
        ifh >> names[i];
    ifh.close();

    // Read elements
    const string iElements = outfile + string("/Elements");
    std::ifstream ife;
    ife.open(iElements.c_str());
    ife >> nElts;
    int nodesPerElt; ife >> nodesPerElt;
    faceData.resize(nElts*nodesPerElt);
    for (int i=0; i<nElts*nodesPerElt; ++i)
        ife >> faceData[i];

    for (int lev=0; lev<NlevPath; ++lev)
    {
        char buf[64];
        sprintf(buf, "/Level_%d", lev);
        std::string ThisStreamFile = outfile + std::string(buf);
        ThisStreamFile += "/Str";
        data[lev] = new MultiFab();
        VisMF::Read(*data[lev],ThisStreamFile);
    }
    //
    // Force other processors to wait until parallel data has been read
    //
    ParallelDescriptor::Barrier();

    // Read element distribution (after we know Nlev from path traces)
    inside_nodes.resize(data.size());
    for (int lev=0; lev<inside_nodes.size(); ++lev)
    {
        inside_nodes[lev].resize(data[lev]->size());

        // Get number of nonempty boxes
        int num_non_zero;
        ife >> num_non_zero;

        for (int j=0; j<num_non_zero; ++j)
        {
            // Get index of non-empty box, and number of entries
            int box_id, num_ids;
            ife >> box_id >> num_ids;

            inside_nodes[lev][box_id].resize(num_ids);

            for (int k=0; k<num_ids; ++k)
                ife >> inside_nodes[lev][box_id][k];
        }
    }
    ife.close();
}

Box 
find_containing_box(const FArrayBox&   XYZ,
                    const Vector<Real>& dx,
                    const Vector<Real>& plo,
                    const Vector<int>&  idXYZ)
{
    // Find bounding box in real space
    IntVect se, be;
    IntVect strIV(IntVect::TheZeroVector());

    Real locL[BL_SPACEDIM];
    Real locH[BL_SPACEDIM];
    for (int d=0; d<BL_SPACEDIM; ++d)
    {
        locL[d] = XYZ(strIV,idXYZ[d]);
        locH[d] = locL[d];
    }

    for (int d=0; d<BL_SPACEDIM; ++d)
    {
        for (int j=1; j<XYZ.box().length(0); ++j)
        {
            strIV = IntVect(D_DECL(j,0,0));
            Real loc = XYZ(strIV,idXYZ[d]);
            locL[d] = std::min(locL[d], loc);
            locH[d] = std::max(locH[d], loc);
        }
         
        se[d] = int( (locL[d] - plo[d])/dx[d] );
        be[d] = int( (locH[d] - plo[d])/dx[d] );

        if (se[d] < 0) Print() << "se < 0: " << se[d] << std::endl;
    }
    return Box(se,be);
}

void
sample_pathlines(Vector<MultiFab*>&       sampledata,
                 int                      sComp,
                 int                      dComp,
                 int                      nComp,
                 const Vector<int>&       comps,
                 AmrData&                 amrData,
                 int                      nGrow,
                 const Vector<MultiFab*>& pathlines,
                 const Vector<Vector<Vector<int> > >& inside_nodes,
                 const Vector<Geometry*>& geoms)
{
    Vector<int> destFillComps(nComp);
    Vector<string> inVarNames(nComp);
    const Vector<string>& plotVarNames = amrData.PlotVarNames();
    for (int i=0; i<nComp; ++i)
    {
        destFillComps[i] = i;
        inVarNames[i] = plotVarNames[comps[sComp+i]];
    }
    const int Nlev = sampledata.size();
    Vector<MultiFab*> state(Nlev);

    // Assumes xyz in first three slots of pathlines
    Vector<int> idXYZ(BL_SPACEDIM);
    for (int d=0; d<BL_SPACEDIM; ++d)
        idXYZ[d] = d;

    const Vector<Real>& plo = amrData.ProbLo();
    for (int lev=0; lev<Nlev; ++lev)
    {
        const int NProcs = ParallelDescriptor::NProcs();        
        const int MyProc = ParallelDescriptor::MyProc();

        const MultiFab& paths = *pathlines[lev];
        const BoxArray& ba = paths.boxArray();

        Vector<int> bx_good(ba.size(),0);
        for (MFIter mfi(paths); mfi.isValid(); ++mfi) {
            if (inside_nodes[lev][mfi.index()].size()) {
                bx_good[mfi.index()] = 1;
            }
        }
        ParallelDescriptor::ReduceIntSum(bx_good.dataPtr(), bx_good.size());
	Vector<int> map_size(NProcs,0);
        const DistributionMapping& dm = paths.DistributionMap();
        for (int i=0; i<ba.size(); ++i) {
            if (bx_good[i]) {
	      map_size[ dm[i] ]++;
            }
        }

        // Build offsets into 1D int array in order to pass new boxes
        Vector<int> map_offset(NProcs,0);
        //Vector<int> bai_offset(NProcs*2*BL_SPACEDIM,0);
        Vector<int> bai_offset(NProcs,0);
        for (int i=1; i<NProcs; ++i) {
            map_offset[i] = map_offset[i-1] + map_size[i-1];
            bai_offset[i] = bai_offset[i-1] + map_size[i-1]*2*BL_SPACEDIM;
        }
        int total_map_size = map_offset[NProcs-1] + map_size[NProcs-1];
        if (total_map_size == 0) continue; // Go to next level...no streamlines on this level

        Vector<int> map_to_full(total_map_size,-1);
        Vector<int> dm_red(total_map_size,-1);
	Vector<int> pcnt(NProcs,0);
        for (int i=0; i<ba.size(); ++i) {
            if (bx_good[i]) {
	      int owner_proc = dm[i];
	      map_to_full[map_offset[owner_proc] + pcnt[owner_proc]] = i;
	      dm_red[map_offset[owner_proc] + pcnt[owner_proc]] = owner_proc;
	      pcnt[owner_proc]++;
            }
        }
	pcnt.clear(); // Just a temporary...

        // fill in my part of the ba_ints data structure
        Vector<int> ba_ints(total_map_size*2*BL_SPACEDIM,0);

        const Vector<Real>& dx = amrData.DxLevel()[lev];
        for (MFIter mfi(paths); mfi.isValid(); ++mfi) {
            if (bx_good[mfi.index()]) {
                Box minBox = find_containing_box(paths[mfi], dx, plo, idXYZ);
                minBox.grow(nGrow);
                for (int d=0; d<BL_SPACEDIM; ++d)
                {
                    ba_ints[bai_offset[MyProc] + d] = minBox.smallEnd(d);
                    ba_ints[bai_offset[MyProc] + BL_SPACEDIM + d] = minBox.bigEnd(d);
                }
                bai_offset[MyProc] += 2*BL_SPACEDIM; // increment offsets for next time
            }            
        }
        ParallelDescriptor::ReduceIntSum(ba_ints.dataPtr(), ba_ints.size());

        // Build the new ba (all procs know all data now)
        BoxArray dba(total_map_size);
        for (int i=0; i<dba.size(); ++i)
        {
            IntVect se( &ba_ints[2*BL_SPACEDIM*i] );
            IntVect be( &ba_ints[2*BL_SPACEDIM*i + BL_SPACEDIM] );
            dba.set(i,Box(se,be));
        }

        if (0 && ParallelDescriptor::IOProcessor()) {
	  std::multimap<int,int> orig;
	  const Vector<int>& dmpm = dm.ProcessorMap();
	  for (int i=0; i<ba.size(); ++i) {
	    orig.insert(std::pair<int,int>(dmpm[i],i));
	  }

	  Vector<int> cov(ba.size(),-1);
	  Box ZBOX(IntVect::TheZeroVector(),IntVect::TheZeroVector());
	  for (int i=0; i<ba.size(); ++i) {
	    if (ba[i] != ZBOX) {
	      cov[i] = dmpm[i];
	    }
	  }

          Print() << "reduced BA:" << std::endl;
	  for (int i=0; i<dba.size(); ++i)
          {
	      Print() << dba[i] << " " << dba[i].length(0)   << " "
                      << dba[i].length(1)   << " " << dba[i].length(2)
                      << "  i,dm_red[i]: " << i << ", " << dm_red[i] << " map_to_full = " << map_to_full[i] << ", grids on this proc: ";
	      typedef std::multimap<int,int>::iterator ii;
	      std::pair<ii,ii> res = orig.equal_range(dm_red[i]);
	      for (ii it=res.first; it!=res.second; ++it) {
                  Print() << it->second << " (" << cov[it->second] << ") ";
	      }
	      Print() << std::endl; 
            }
        }

        // build container to hold data for interp, and fill it with AmrData functions
        DistributionMapping dmRed(dm_red);
        MultiFab data(dba,dmRed,nComp,0);
        data.setVal(-20000);

        BoxArray vba = amrex::intersect(dba,amrData.ProbDomain()[lev]);
        for (int i=0; i<vba.size(); ++i)
        {
            if (!vba[i].ok()) {
                Print() << "bad vba box: " << vba[i] << ", " << i << std::endl;
	    }
        }
        MultiFab dataS(vba,DistributionMapping(vba),1,0);

        strt_io = ParallelDescriptor::second();
        for (int i = 0; i < inVarNames.size(); i++)
        {
            amrData.FillVar(dataS,lev,inVarNames[i],0);
            data.copy(dataS,0,destFillComps[i],1);
        }
        io_time += ParallelDescriptor::second() - strt_io;

        const Geometry& geom = *geoms[lev];
        const Box& pd = amrData.ProbDomain()[lev];
        if (geom.isAnyPeriodic())
        {
            Vector<IntVect> shifts;
            BoxList sbl;
            BoxList bl;
            for (int i=0; i<dba.size(); ++i)
            {
                geom.periodicShift(pd,dba[i],shifts);
                for (int j=0; j<shifts.size(); ++j)
                {
                    Box isect = (dba[i] + shifts[j]) & pd;
                    sbl.push_back(isect);
                    bl.push_back(isect - shifts[j]);
                }
            }
            if (sbl.size()>0)
            {
                BoxArray sba(sbl);
                BoxArray uba(bl);
                DistributionMapping dmp(uba);
                MultiFab sh_pieces(sba,dmp,nComp,0);
                MultiFab unsh_pieces(uba,dmp,nComp,0);
                sh_pieces.setVal(-400);
                // Fill shifted mf...guaranteed to lay completely inside pd
                strt_io = ParallelDescriptor::second();
                amrData.FillVar(sh_pieces,lev,inVarNames,destFillComps);
                io_time += ParallelDescriptor::second() - strt_io;

                // Now copy to mirror mf...where box-by-box, pieces are in periodic region outside pd
                for (MFIter mfi(sh_pieces); mfi.isValid(); ++mfi)
                {
                    const FArrayBox& src = sh_pieces[mfi];
                    FArrayBox& dst = unsh_pieces[mfi];
                    dst.copy(src,src.box(),0,dst.box(),0,nComp);
                }
                
                // Finally, do a mf-to-mf copy
                data.copy(unsh_pieces,0,0,nComp);
            }
        }

        Print() << "Sampling paths for level " << lev << std::endl;

        for (MFIter mfi(data); mfi.isValid(); ++mfi)
        {
            const FArrayBox& pth = (*pathlines[lev])[ map_to_full[mfi.index()] ];
            const FArrayBox& src = data[mfi];
            FArrayBox& sam = (*sampledata[lev])[  map_to_full[mfi.index()] ];
            int nCompPath = pth.nComp();
            
            interpstream(BL_TO_FORTRAN_ANYD(pth), &nCompPath,
                         BL_TO_FORTRAN_ANYD(src), &nComp,
                         BL_TO_FORTRAN_N_ANYD(sam,dComp),
                         dx.dataPtr(), plo.dataPtr());
        }
        
        Print() << "....paths sampled for level, comp, ncomp " << lev << ", " << sComp << ", " << nComp << std::endl;
    }

    for (int i=0; i<nComp; i++)
    {
        amrData.FlushGrids(comps[sComp+i]);
    }
}
    
void
set_sample_distance(Vector<MultiFab*>&       result,
                    int                      dComp,
                    const Vector<MultiFab*>& path)
{
    const int Nlev = result.size();
    for (int lev=0; lev<Nlev; ++lev)
    {
        for (MFIter mfi(*result[lev]); mfi.isValid(); ++mfi)
        {
            const FArrayBox& pth = (*path[lev])[mfi];
            FArrayBox& res = (*result[lev])[mfi];
            set_distance(BL_TO_FORTRAN_ANYD(pth),
                         BL_TO_FORTRAN_N_ANYD(res,dComp));
        }
    }
}

void
set_sample_location(Vector<MultiFab*>&       result,
                    int                      dComp,
                    const Vector<MultiFab*>& path)
{
    for (int lev=0; lev<result.size(); ++lev)
        result[lev]->copy(*path[lev],0,dComp,BL_SPACEDIM);
}

