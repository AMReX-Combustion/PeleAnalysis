#include <string>
#include <iostream>
#include <set>
#include <map>

#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include <AMReX_DataServices.H>

#include <AMReX_BLFort.H>

#include <StreamData.H>

using namespace amrex;

using std::vector;
using std::map;
using std::string;
using std::cerr;
using std::endl;
using std::cout;
using std::ofstream;
using std::ios;

static int NLINESMAX=32700;



typedef std::list<std::pair<int,int> > SelectionMapList;

void
downsampleStreamData(const StreamData&                stream,
                     int&                             num_lines,
                     Vector<BoxArray>&                 stream_data_layout,
                     Vector<map<int,int> >&            map_selected_to_full_gidx,
                     Vector<Vector<SelectionMapList> >& selectIdxPairs,
                     bool                             verbose)
{
    // This routine is parallel (all do the same thing), but only IOProc will write anything
    verbose = (ParallelDescriptor::IOProcessor() ? verbose : false);

    double rval = (double)num_lines / stream.NLines(); // fraction of total to write
    int nLev = stream.FinestLevel()+1;

    selectIdxPairs.clear();
    selectIdxPairs.resize(nLev);

    stream_data_layout.clear();
    stream_data_layout.resize(nLev);

    map_selected_to_full_gidx.clear();
    map_selected_to_full_gidx.resize(nLev);

    int slo = stream.StreamIdxLo();
    int shi = stream.StreamIdxHi();

    int cnt = 0;
    for (int lev=0; lev<nLev; ++lev)
    {
        const BoxArray& ba = stream.boxArray(lev);
        int nBoxes = ba.size();
        BoxList bldst;
        int last_added_srcboxid=-1;

        // Easiest to build this first on map based on full ba indexing
        map<int,SelectionMapList> map_baIdx_selectIdxPairs;
        map<int,int>& m_sel_to_full = map_selected_to_full_gidx[lev];

        for (int box_id=0; box_id<nBoxes; ++box_id)
        {
            const Vector<int>& in_lev = stream.InsideNodes(lev,box_id);
            for (int i=0; i<in_lev.size(); ++i)
            {
                double doit = Random(); // Note: MUST override default cpu-dependent seed!
                if (doit<rval)
                {
                    if (box_id!=last_added_srcboxid) {
                        last_added_srcboxid = box_id;
                        m_sel_to_full[m_sel_to_full.size()] = box_id;
                    }
                    int old_local_idx = i;
                    int new_global_idx = cnt;
                    map_baIdx_selectIdxPairs[box_id].push_back(
                        std::pair<int,int>(new_global_idx,old_local_idx));
                    cnt++;
                }
            }
            int num_this_bx = map_baIdx_selectIdxPairs[box_id].size();
            if (num_this_bx>0)
            {
                int first_idx = map_baIdx_selectIdxPairs[box_id].front().first;
                int last_idx = map_baIdx_selectIdxPairs[box_id].back().first;
                Box bx(IntVect(D_DECL(first_idx,slo,0)),
                       IntVect(D_DECL(last_idx,shi,0)));
                bldst.push_back(bx);
            }
        }
        stream_data_layout[lev] = BoxArray(bldst);

        // Build final selection map based on reduced ba
        int nSrcGrids = m_sel_to_full.size();
        selectIdxPairs[lev].resize(nSrcGrids);
        for (int k=0; k<nSrcGrids; ++k) {
            int fgidx = m_sel_to_full[k];
            selectIdxPairs[lev][k] = map_baIdx_selectIdxPairs[fgidx];
        }

        if (verbose)
        {
            for (int k=0; k<nSrcGrids; ++k) {
                const SelectionMapList& idx_map = selectIdxPairs[lev][k];
                int fgidx = m_sel_to_full[k];
                int nAvail = ba[fgidx].length(0);
                cout << k << " (" << idx_map.size() << "): " ;
                for (SelectionMapList::const_iterator it=idx_map.begin(); it!=idx_map.end(); ++it)
                {
                    cout << " ( " << it->first << " from " << it->second << " ) ";
                    BL_ASSERT(it->second < nAvail);
                }
                cout << " src data has " << nAvail << " particles" << endl;
            }
        }
    }
    num_lines = cnt;
    ParallelDescriptor::Barrier(); // Make sure they all get here so nLines is uniform across procs
}

#ifdef BL_WRITE_BINARY

#include "TECIO.h"
#define SIZET INTEGER4
void
write_tec_binary(const FArrayBox&     strm,
                 const string&        outfile,
                 const string&        file_label,
                 const Vector<string>& names,
                 const Vector<int>&    write_lines)
{
    ofstream osf;
    osf.open(outfile.c_str(),std::ios::out);

    int Ncomp = names.size();
    BL_ASSERT(Ncomp==strm.nComp());

    BL_ASSERT(write_lines.size()==strm.box().length(0));

    std::string vars = names[0];
    for (int j=1; j<Ncomp; ++j)
        vars += " " + names[j];

    bool verbose = true;

    INTEGER4 Debug = std::max(0,verbose-1);
    INTEGER4 VIsDouble = 1;
    TECINI100((char*)file_label.c_str(),
              (char*)vars.c_str(),
              (char*)binary_outfile.c_str(),
              (char*)".",
              &Debug,
              &VIsDouble);
    
    INTEGER4 ZoneType  = 0;  /* ORDERED */
    INTEGER4 ICellMax  = 0;  /* UNUSED (for FEBRICK ONLY) */
    INTEGER4 JCellMax  = 0;  /* UNUSED (for FEBRICK ONLY) */
    INTEGER4 KCellMax  = 0;  /* UNUSED (for FEBRICK ONLY) */
    INTEGER4 IsBlock   = 1;
    INTEGER4 NumFaceConnections = 0;
    INTEGER4 FaceNeighborMode   = 0;
    INTEGER4 ShareConnectivityFromZone = 0; /* UNUSED (for FE ZONES ONLY) */

    Vector<Real> loc(BL_SPACEDIM);
    const Box& box = strm.box();
    const IntVect ivst = box.smallEnd();
    const IntVect ivmid(D_DECL(ivst[0],0,ivst[2]));
    
    int nLines = box.length(0);
    int Npts = box.length(1);
    INTEGER4 IMax = Npts;  /* Max number of nodes in I index direction */
    INTEGER4 JMax = 1;     /* Max number of nodes in J index direction */
    INTEGER4 KMax = 1;     /* Max number of nodes in K index direction, used for ORDERED data */

    Vector<Real> dat(Npts*Ncomp);
    const IntVect ivd = box.smallEnd();
    for (int i=0; i<nLines; ++i)
    {
        if (write_lines[i] > 0)
        {
            int cnt=0;
            for (int n=0; n<Ncomp; ++n)
            {
                for (int j=0; j<Npts; ++j)
                {
                    dat[cnt++] = strm(ivst+IntVect(i,j,0),n);
                }
            }
                            
            char buf[72];
            sprintf(buf,"id%d",i);
            string label = std::string(buf);
            
            if (verbose>1)
                std::cerr << "Writing line: " << i << std::endl;
            
            TECZNE100((char*)label.c_str(),
                      &ZoneType,
                      &IMax,
                      &JMax,
                      &KMax,
                      &ICellMax,
                      &JCellMax,
                      &KCellMax,
                      &IsBlock,
                      &NumFaceConnections,
                      &FaceNeighborMode,
                      NULL,      /* ValueLocation */
                      NULL,      /* ShareVarFromZone */
                      &ShareConnectivityFromZone);        
            
            INTEGER4 III = dat.size();
            TECDAT100(&III,dat.dataPtr(),&VIsDouble);
        }
    }
    TECEND100(); 

#if 0        
        
        // Add distance
        int doffset=ivmid[1] - ivst[1];
        
        tmp(ivd + doffset*BASISV(0),Ncomp) = 0.;
        
        int iprev = 0;
        for (int i=-1; i>=box.smallEnd(1); --i)
        {
            IntVect ivsrc = ivmid + j*BASISV(0) + i*BASISV(1);
            IntVect ivsrcprev = ivmid + j*BASISV(0) + iprev*BASISV(1);
            IntVect ivdst = ivd + (doffset+i)*BASISV(0);
            IntVect ivdstprev = ivd + (doffset+iprev)*BASISV(0);
            
            tmp(ivdst,Ncomp) = tmp(ivdstprev,Ncomp) -
                std::sqrt( std::pow( strm(ivsrc,0)-strm(ivsrcprev,0),2)
                           + std::pow( strm(ivsrc,1)-strm(ivsrcprev,1),2)
                           + std::pow( strm(ivsrc,2)-strm(ivsrcprev,2),2) );
            
            iprev = i;
        }
        
        iprev = 0;
        for (int i=1; i<=box.bigEnd(1); ++i)
        {
            IntVect ivsrc = ivmid + j*BASISV(0) + i*BASISV(1);
            IntVect ivsrcprev = ivmid + j*BASISV(0) + iprev*BASISV(1);
            IntVect ivdst = ivd + (doffset+i)*BASISV(0);
            IntVect ivdstprev = ivd + (doffset+iprev)*BASISV(0);
            
            tmp(ivdst,Ncomp) = tmp(ivdstprev,Ncomp) +
                std::sqrt( std::pow( strm(ivsrc,0)-strm(ivsrcprev,0),2)
                           + std::pow( strm(ivsrc,1)-strm(ivsrcprev,1),2)
                           + std::pow( strm(ivsrc,2)-strm(ivsrcprev,2),2) );
            iprev = i;
        }
        
        INTEGER4 III = (Ncomp+1) * Npts;
#endif
}

#else
void
write_tec_ascii(const FArrayBox&     strm,
                const string&        outfile,
                const string&        file_label,
                const Vector<string>& names,
                const Vector<int>&    write_lines)
{
    ofstream osf;
    osf.open(outfile.c_str(),std::ios::out);

    int nComp = names.size();
    BL_ASSERT(nComp==strm.nComp());

    osf << "VARIABLES = ";
    for (int i=0; i<nComp; ++i)
        osf << names[i] << " ";
    osf << '\n';


    BL_ASSERT(write_lines.size()==strm.box().length(0));

    std::string vars = names[0];
    for (int j=1; j<nComp; ++j)
        vars += " " + names[j];
    
    const Box& box = strm.box();
    int nLines = box.length(0);
    int Npts = box.length(1);

    const IntVect ivst = box.smallEnd();
    for (int i=0; i<nLines; ++i)
    {
        if (write_lines[i] > 0)
        {
            char buf[72];
            sprintf(buf,"id%d",i);
            string label = std::string(buf);

            osf << "ZONE T=" << label << " I=" << Npts << " F=POINT" << endl;;                    
            for (int j=0; j<Npts; ++j)
            {
              IntVect ivt = ivst + IntVect(D_DECL(i,j,0));
              for (int n=0; n<nComp; ++n)
                osf << strm(ivt,n) << " ";
              osf << '\n';
            }
        }
    }

    osf.close();
}
#endif

static bool do_test(Real thisVal,
                    const std::string& sgn,
                    Real critVal);

int
main (int   argc,
      char* argv[])
{
    amrex::Initialize(argc,argv);
    //
    // Force everyone to have the same random sequence.
    //
    InitRandom(987654321);

    ParmParse pp;

    std::string infile; pp.get("infile",infile);
    std::string outfile; pp.get("outfile",outfile);

    int verbose=0; pp.query("verbose",verbose);

    StreamData stream(infile);
    int finestLevel=stream.FinestLevel(); pp.query("finestLevel",finestLevel);

    // How many lines in file
    if (ParallelDescriptor::IOProcessor() && verbose)
        std::cerr << "Found " << stream.NLines() << " streamlines in this file." << std::endl;
    
    const Vector<string>& names = stream.ComponentNames();
    Vector<int> comps;
    if (int nc = pp.countval("comps"))
    {
        comps.resize(nc);
        pp.getarr("comps",comps,0,nc);
    }
    else
    {
        int sComp = 0; pp.query("sComp",sComp);
        int nComp = names.size(); pp.query("nComp",nComp);
        BL_ASSERT(sComp+nComp <= names.size());
        comps.resize(nComp);
        for (int i=0; i<nComp; ++i)
            comps[i] = sComp + i;
    }

    Vector<string> selectedNames(comps.size());
    for (int i=0; i<comps.size(); ++i) {
        selectedNames[i] = names[comps[i]];
    }

    int distComp=-1; pp.query("distComp",distComp);
    Real distVal; 
    if (distComp>=0) {
        pp.get("distVal",distVal);
        BL_ASSERT(distComp<comps.size());
        char buf[64]; sprintf(buf,"%g",distVal);
        string distName = "distance_from_" + names[comps[distComp]] + "_eq_" + string(buf);
        selectedNames.resize(selectedNames.size()+1);
        selectedNames[selectedNames.size()-1] = distName;
    }

    bool no_filter = false;
    pp.query("no_filter",no_filter); // Explicitly shut off all filters


    int nLines = 0; pp.query("nLines",nLines);
    nLines = std::min(nLines,NLINESMAX);

    if (ParallelDescriptor::IOProcessor() && verbose)
        std::cerr << "  Looking to extract approximately " << nLines << " of them." << std::endl;

    Real RXY = -1; pp.query("RXY",RXY);
    string RXYsgn = ""; pp.query("RXYsgn",RXYsgn);

    Vector<Real> maxVals;
    Vector<string> maxSgns;
    Vector<int> maxComps;
    if (int nmax = pp.countval("maxComps"))
    {
        maxComps.resize(nmax);
        maxVals.resize(nmax);
        maxSgns.resize(nmax);
        pp.getarr("maxComps",maxComps,0,nmax);
        pp.getarr("maxVals",maxVals,0,nmax);
        pp.getarr("maxSgns",maxSgns,0,nmax);
    }

    Vector<Real> minVals;
    Vector<string> minSgns;
    Vector<int> minComps;
    if (int nmin = pp.countval("minComps"))
    {
        minComps.resize(nmin);
        minVals.resize(nmin);
        minSgns.resize(nmin);
        pp.getarr("minComps",minComps,0,nmin);
        pp.getarr("minVals",minVals,0,nmin);
        pp.getarr("minSgns",minSgns,0,nmin);
    }

    // Where fab(atComp)=valAt, fab(compAt) atSgn atVal
    Vector<Real> compAt;
    Vector<Real> valAt;
    Vector<Real> atVal;
    Vector<string> atSgns;
    Vector<int> atComps;
    if (int nat = pp.countval("atComps"))
    {
        compAt.resize(nat);
        atComps.resize(nat);
        valAt.resize(nat);
        atVal.resize(nat);
        atSgns.resize(nat);
        pp.getarr("compAt",compAt,0,nat);
        pp.getarr("atComps",atComps,0,nat);
        pp.getarr("valAt",valAt,0,nat);
        pp.getarr("atVal",atVal,0,nat);
        pp.getarr("atSgns",atSgns,0,nat);
    }

    string condStr = "";
    if (ParallelDescriptor::IOProcessor() && verbose)
    {
        cout << "Preparing to extract/write the following components (i,idx,name): " << endl;
        for (int i=0; i<comps.size(); ++i) {
            cout << "(" << i << "," << comps[i] << "," << names[comps[i]] << ") ";
        }
        cout << endl;

        if (selectedNames.size() > comps.size()) {
            cout << "....appending the following component(s): ";
            for (int i=comps.size(); i<selectedNames.size(); ++i) {
                cout << selectedNames[i] << endl;
            }
        }

        for (int i=0; i<maxComps.size(); ++i) {
            char buf[64]; sprintf(buf,"%g",maxVals[i]);
            condStr += "Max(" + names[comps[maxComps[i]]] + ") " + maxSgns[i] + " " + string(buf) + "; ";
        }

        for (int i=0; i<minComps.size(); ++i) {
            char buf[64]; sprintf(buf,"%g",minVals[i]);
            condStr += "Min(" + names[comps[minComps[i]]] + ") " + minSgns[i] + " " + string(buf) + "; ";
        }

        for (int i=0; i<atComps.size(); ++i) {
            char buf1[64]; sprintf(buf1,"%g",atVal[i]);
            char buf2[64]; sprintf(buf2,"%g",valAt[i]);

            condStr += names[comps[compAt[i]]] + "(" + names[comps[atComps[i]]]
                + " = " + string(buf1) + ") " + atSgns[i] + " " + string(buf2) + "; ";
        }

        if (RXY>0) {
            char buf[64]; sprintf(buf,"%g",RXY);
            condStr += "radiusXY " + RXYsgn + " " + string(buf) + " ";
        }

        if (condStr != "")
            std::cout << " ...subject to: " << condStr << std::endl;

    }

    Vector<BoxArray> stream_data_layout;
    Vector<Vector<SelectionMapList> > selectIdxPairs;
    Vector<map<int,int> > map_selected_to_full_gidx;

    downsampleStreamData(stream,nLines,stream_data_layout,
                         map_selected_to_full_gidx,selectIdxPairs,verbose);
    if (ParallelDescriptor::IOProcessor()) {
        cout << "Reduced dataset has " << nLines << " lines " << std::endl;
    }

    int slo = stream.StreamIdxLo();
    int shi = stream.StreamIdxHi();

    Box finalBox(IntVect(D_DECL(0,slo,0)),
                 IntVect(D_DECL(nLines-1,shi,0)));
    BoxArray finalBa(finalBox);
    DistributionMapping dmfinal(finalBa);
    MultiFab finalMF(finalBa,dmfinal,selectedNames.size(),0);

    for (int lev=0; lev<=finestLevel; ++lev)
    {
        const BoxArray& dstBa=stream_data_layout[lev];

        if (ParallelDescriptor::IOProcessor() && verbose)
        {
            cout << "Reduced  dataset has " << dstBa.size() << " boxes (original has "
                 << stream.boxArray(lev).size() << ") at level = " << lev << endl; 
        }
        if (dstBa.size()>0) 
        {
            MultiFab dst(dstBa,DistributionMapping(dstBa),comps.size(),0);
            for (MFIter mfi(dst); mfi.isValid(); ++mfi)
            {
                int gidx = mfi.index();
                FArrayBox& dstFab = dst[mfi];
                const SelectionMapList& line_map = selectIdxPairs[lev][gidx];
                int srcIdx = map_selected_to_full_gidx[lev][gidx];
                
                for (int n=0; n<comps.size(); ++n)
                {
                    const FArrayBox& srcFab = *(stream.getFab(lev,srcIdx,comps[n]));
                    for (SelectionMapList::const_iterator it = line_map.begin(); it!=line_map.end(); ++it)
                    {
                        int old_local_idx = it->second;                
                        Box srcBox(IntVect(D_DECL(old_local_idx,slo,0)),
                                   IntVect(D_DECL(old_local_idx,shi,0)));
                        
                        int new_global_idx = it->first;
                        BL_ASSERT(new_global_idx < nLines);
                        Box dstBox(IntVect(D_DECL(new_global_idx,slo,0)),
                                   IntVect(D_DECL(new_global_idx,shi,0)));
                        
                        BL_ASSERT(srcFab.box().contains(srcBox));
                        BL_ASSERT(dstFab.box().contains(dstBox));
                        
                        dstFab.copy(srcFab,srcBox,0,dstBox,n,1);
                    }
                }
            } // mfi
            
            finalMF.copy(dst);
        }
    }
    
    int procWithFab = finalMF.DistributionMap()[0];

    if (ParallelDescriptor::MyProc()==procWithFab)
    {
        FArrayBox& fab=finalMF[0];

        // Scan downselected streamlines for criteria
        Vector<int> write_lines(nLines,1);

        if (no_filter) {
            for (int i=0; i<write_lines.size(); ++i) {
                write_lines[i] = 1;
            }
        }
        else
        {
            IntVect se=finalBox.smallEnd();
            int nvals = shi-slo+1;

            // Check criteria on maximum values
            for (int n=0; n<maxComps.size(); ++n) 
            {
                int nc = maxComps[n];
                for (int i=0; i<nLines; ++i)
                {
                    Real maxVal=fab(se,nc);
                    for (int j=0; j<nvals; ++j) {
                        IntVect iv = se + i*BASISV(0) + j*BASISV(1);
                        maxVal = std::max(fab(iv,nc),maxVal);
                    }
                    if ( ! do_test(maxVal,maxSgns[n],maxVals[n]) ) {
                        write_lines[i] = 0;
                    }
                }
            }

            // Check criteria on minimum values
            for (int n=0; n<minComps.size(); ++n) 
            {
                int nc = minComps[n];
                for (int i=0; i<nLines; ++i)
                {
                    Real minVal=fab(se,nc);
                    for (int j=0; j<nvals; ++j) {
                        IntVect iv = se + i*BASISV(0) + j*BASISV(1);
                        minVal = std::min(fab(iv,nc),minVal);
                    }
                    if ( ! do_test(minVal,minSgns[n],minVals[n]) ) {
                        write_lines[i] = 0;
                    }
                }
            }

            // Check criteria on radius of seed point
            BL_ASSERT(comps.size()>2);
            if (RXY>0)
            {
                for (int i=0; i<nLines; ++i)
                {
                    IntVect ivs(D_DECL(i,0.5*(slo+shi),0));
                    Real radius = 0;
                    for (int n=0; n<2; ++n) {
                        Real val = fab(ivs,n);
                        radius += val*val;
                    }
                    radius = std::sqrt(radius);
                    if ( ! do_test(radius,RXYsgn,RXY) ) {
                        write_lines[i] = 0;
                    }
                }
            }

            // Check criteria on "atVals"
            for (int n=0; n<atComps.size(); ++n) 
            {
                int loc_comp = atComps[n];
                Real loc_val = atVal[n];
                int test_comp = compAt[n];
                Real test_val = valAt[n];

                for (int i=0; i<nLines; ++i)
                {
                    bool found_it = false;
                    IntVect ivl, ivh;
                    Real loc_val_lo, loc_val_hi;
                    for (int j=0; j<nvals-1 && !found_it; ++j) {
                        ivl = se + i*BASISV(0) + j*BASISV(1);
                        ivh = se + i*BASISV(0) + (j+1)*BASISV(1);
                        loc_val_lo = fab(ivl,loc_comp);
                        loc_val_hi = fab(ivh,loc_comp);
                        if ( ( (loc_val_lo > loc_val)  && (loc_val_hi < loc_val) ) ||
                             ( (loc_val_lo < loc_val)  && (loc_val_hi > loc_val) ) ) {
                            Real alpha = (loc_val - loc_val_lo)/(loc_val_hi-loc_val_lo);
                            BL_ASSERT(alpha>=0 && alpha<=1);
                            Real val_lo = fab(ivl,test_comp);
                            Real val_hi = fab(ivh,test_comp);
                            Real val_test = val_lo + alpha*(val_hi-val_lo);
                            write_lines[i] = do_test(val_test,atSgns[n],test_val);
                            found_it = true;
                        }
                    }
                }
            }

            // Add auxiliary quantities
            if (distComp>=0)
            {
                // Compute relative to the start of the streamline, then shift
                int nDist = comps.size();
                int loc_comp = distComp;
                Real loc_val = distVal;
                Real loc_offset = 0;

                IntVect ivl, ivh;
                BL_ASSERT(comps.size()>BL_SPACEDIM);
                for (int i=0; i<nLines; ++i)
                {
                    ivl = se + i*BASISV(0);
                    fab(ivl,nDist) = 0;
                    for (int j=1; j<nvals; ++j) {
                        ivl = se + i*BASISV(0) + (j-1)*BASISV(1);
                        ivh = se + i*BASISV(0) +   j  *BASISV(1);
                        Real ds = 0;
                        for (int n=0; n<BL_SPACEDIM; ++n) {
                            Real dx = fab(ivh,n) - fab(ivl,n);
                            ds += dx*dx;
                        }
                        ds = std::sqrt(ds);
                        fab(ivh,nDist) = fab(ivl,nDist) + ds;
                    }                    

                    bool found_it = false;
                    for (int j=0; j<nvals-1 && !found_it; ++j) {
                        ivl = se + i*BASISV(0) + j*BASISV(1);
                        ivh = se + i*BASISV(0) + (j+1)*BASISV(1);
                        Real loc_val_lo = fab(ivl,loc_comp);
                        Real loc_val_hi = fab(ivh,loc_comp);
                        if ( ( (loc_val_lo > loc_val)  && (loc_val_hi < loc_val) ) ||
                             ( (loc_val_lo < loc_val)  && (loc_val_hi > loc_val) ) ) {
                            Real alpha = (loc_val - loc_val_lo)/(loc_val_hi-loc_val_lo);
                            BL_ASSERT(alpha>=0 && alpha<=1);
                            Real val_lo = fab(ivl,nDist);
                            Real val_hi = fab(ivh,nDist);
                            loc_offset = val_lo + alpha*(val_hi-val_lo);
                            // Do the shift
                            for (int jj=0; jj<nvals; ++jj) {
                                IntVect iv = se + i*BASISV(0) + jj*BASISV(1);
                                fab(iv,nDist) = fab(iv,nDist) - loc_offset;
                            }                        
                            found_it = true;
                        }
                    }

                    if (!found_it) {
                        // set all distance values to something off the line
                        Real err_val=fab(IntVect(se + i*BASISV(0) + (nvals-1)*BASISV(1)),nDist) * 2;
                        for (int j=0; j<nvals; ++j) {
                            IntVect iv = se + i*BASISV(0) + j*BASISV(1);
                            fab(iv,nDist) = err_val;
                        }
                    }
                }


            }
        }

        char buf1[64]; sprintf(buf1,"%d",nLines);
        char buf2[64]; sprintf(buf2,"%d",stream.NLines());
        string file_label="streamlines from " + infile + " selected from a subset ("
            + string(buf1) + " of " + string(buf2) + ")";
        if (condStr != "")
            file_label += " and further conditioned such that: " + condStr;
#ifdef BL_WRITE_BINARY
        write_tec_binary(fab,outfile,file_label,selectedNames,write_lines);
#else
        write_tec_ascii(fab,outfile,file_label,selectedNames,write_lines);
#endif
    }
    amrex::Finalize();
    return 0;
}

static bool do_test(Real thisVal,
                    const std::string& sgn,
                    Real critVal)
{
    bool retVal = false;
    if (sgn=="ge") {
        retVal =  thisVal >= critVal;
    } else if (sgn=="gt") {
        retVal =  thisVal > critVal;
    } else if (sgn=="lt") {
        retVal =  thisVal < critVal;
    } else if (sgn=="le") {
        retVal =  thisVal <= critVal;
    } else if (sgn=="eq") {
        retVal =  thisVal == critVal;
    } else if (sgn=="ne") {
        retVal =  thisVal != critVal;
    }
    return retVal;
}

