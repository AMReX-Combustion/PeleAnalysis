#include <string>
#include <iostream>
#include <set>

#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include <AMReX_DataServices.H>
#include <AMReX_BCRec.H>
#include <AMReX_Interpolater.H>
#include <AppendToPlotFile.H>
#include <WritePlotFile.H>

#include <AMReX_BLFort.H>

using namespace amrex;
//typedef StateDescriptor::BndryFunc BndryFunc;

extern "C" {
    void pushvtog(const int* lo,  const int* hi,
                  const int* dlo, const int* dhi,
                  Real* U, const int* Ulo, const int* Uhi,
                  const int* nc);
    void normalize(const int* lo,  const int* hi,
                   const Real *U, const int* U_lo, const int* U_hi,
                   Real* S, const int* S_lo, const int* S_hi,
                   const Real* nmin, const Real* nmax, const Real* dx);
    void smooth(const int* lo,  const int* hi,
                const Real *Tin, const int* Tin_lo, const int* Tin_hi,
                Real *Tout, const int* Tout_lo, const int* Tout_hi);
    void mcurv(const int* lo,  const int* hi,
               const Real *T, const int* T_lo, const int* T_hi,
               Real *curv, const int* curv_lo, const int* curv_hi,
               Real *wrk, const int* wrk_lo, const int* wrk_hi,
               const Real* delta, const int* sym);
    void gcurv(const int* lo,  const int* hi,
               const Real *T, const int* T_lo, const int* T_hi,
               Real *curv, const int* curv_lo, const int* curv_hi,
               const Real* delta);
    void strainrate(const int* lo,  const int* hi,
                    const Real *U, const int* U_lo, const int* U_hi,
                    const Real *T, const int* T_lo, const int* T_hi,
                    Real *sr, const int* sr_lo, const int* sr_hi,
                    Real *wrk, const int* wrk_lo, const int* wrk_hi,
                    const Real* delta);
    void straintensor(const int* lo,  const int* hi,
                      const Real *U, const int* U_lo, const int* U_hi,
                      Real *sr, const int* sr_lo, const int* sr_hi,
                      const Real* delta);
    void progressgrad(const int* lo,  const int* hi,
                      const Real *c, const int* c_lo, const int* c_hi,
                      Real *g, const int* g_lo, const int* g_hi,
                      const Real* delta);
}

std::string
getFileRoot(const std::string& infile)
{
    vector<std::string> tokens = Tokenize(infile,std::string("/"));
    return tokens[tokens.size()-1];
}

void
FillCFgrowCells(int refRatio, const BoxArray& fine_ba,Vector<Geometry*>& geoms, int lev,
                MultiFab& crseSrc, int sComp, int nComp, int nGrow, MultiFab& growScr)
{
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

int
main (int   argc,
      char* argv[])
{
    amrex::Initialize(argc,argv);
    {
    // if (argc < 2)
    //    print_usage(argc,argv);
    
    ParmParse pp;
    
    // if (pp.contains("help"))
    //    print_usage(argc,argv);

    int verbose = 0; pp.query("verbose",verbose);
    if (verbose>1)
       AmrData::SetVerbose(true);
    
    std::string plotFileName; pp.get("infile",plotFileName);
    if (ParallelDescriptor::IOProcessor())
        std::cout << "infile = " << plotFileName << std::endl;
    std::string outfile(getFileRoot(plotFileName) + "_K"); pp.query("outfile",outfile);
    std::cout << "reading plt file = " << plotFileName << std::endl;
    
    DataServices::SetBatchMode();
    Amrvis::FileType fileType(Amrvis::NEWPLT);

    DataServices dataServices(plotFileName, fileType);
    if( ! dataServices.AmrDataOk()) {
        std::cout << "Somethin' in wrong with plt file = " << plotFileName << std::endl;
        DataServices::Dispatch(DataServices::ExitRequest, NULL);
    }
    AmrData& amrData = dataServices.AmrDataRef();

    std::cout << "Trying to get finestLevel" << std::endl;
    int finestLevel = amrData.FinestLevel();
    pp.query("finestLevel",finestLevel);
    int Nlev = finestLevel + 1;

    std::cout << "Trying to get data for progvar" << std::endl;
    std::string progressName = "temp"; pp.query("progressName",progressName);
    Real progMin = 1.0e20; pp.query("progMin",progMin);
    Real progMax = -1.0e20; pp.query("progMax",progMax);
    int floorIt = 0; pp.query("floorIt",floorIt);

    int useFileMinMax = 1; pp.query("useFileMinMax",useFileMinMax);
    Real progMinlvl = 1.0e20;
    Real progMaxlvl = -1.0e20;
    if (useFileMinMax || floorIt)
    {
        if (useFileMinMax) {
            for (int lev=0; lev<Nlev; ++lev) {
                amrData.MinMax(amrData.ProbDomain()[lev], progressName, lev, progMinlvl, progMaxlvl);
                progMin = std::min(progMin,progMinlvl);
                progMax = std::max(progMax,progMaxlvl);
                std::cout << "In there" << std::endl;
            }
            ParallelDescriptor::ReduceRealMin(progMin);
            ParallelDescriptor::ReduceRealMax(progMax);
        }

        std::cout << "progressName = " << progressName << " at index: " << amrData.StateNumber(progressName) << std::endl;
        std::cout << "useFileMinMax = " << useFileMinMax << std::endl;
        std::cout << "Min/Max = " << progMin << " / " << progMax << std::endl;

        ParallelDescriptor::Barrier();

        if (progMin >= progMax) {
            amrex::Abort("progMin must be less than progMax");
        }
    }
    std::cout << "Data for progvar read OK" << std::endl;

    int idC = -1;
    const Vector<std::string>& plotVarNames = amrData.PlotVarNames();
    for (int i=0; i<plotVarNames.size(); ++i)
    {
        if (plotVarNames[i] == progressName) idC = i;
    }
    if (idC<0)
        Print() << "Cannot find required data in pltfile" << std::endl;

    const int idCst = 0;
    const int idVst = idCst + 1;
    int nCompIn = idVst;

    Vector<std::string> inVarNames(nCompIn);
    inVarNames[idCst] = plotVarNames[idC];
    bool do_strain = false;
    pp.query("do_strain",do_strain);
    if (do_strain) {
      nCompIn += BL_SPACEDIM;
      inVarNames.resize(nCompIn);
      inVarNames[idVst+0] = "x_velocity";
      inVarNames[idVst+1] = "y_velocity";
#if BL_SPACEDIM>2
      inVarNames[idVst+2] = "z_velocity";
#endif
    }
    Vector<int> destFillComps(nCompIn);
    for (int i=0; i<nCompIn; ++i)
        destFillComps[i] = i;    

    Vector<int> auxComps;
    if (int nc = pp.countval("aux_comps"))
    {
        auxComps.resize(nc);
        pp.getarr("aux_comps",auxComps,0,nc);
    }
    else
    {
        int sComp = 0;
        pp.query("aux_sComp",sComp);
        int nComp = 0;
        pp.query("aux_nComp",nComp);
        BL_ASSERT(sComp+nComp <= amrData.NComp());
        auxComps.resize(nComp);
        for (int i=0; i<nComp; ++i)
            auxComps[i] = sComp + i;
    }
    if (auxComps.size()>0)
    {
        int nCompNEW = inVarNames.size()+auxComps.size();
        inVarNames.resize(nCompNEW);
        destFillComps.resize(nCompNEW);
        for (int i=0; i<auxComps.size(); ++i)
        {
            inVarNames[nCompIn+i] = plotVarNames[auxComps[i]];
            destFillComps[nCompIn+i] = nCompIn+i;
        }
        nCompIn += auxComps.size();
    }

    // Enforce symmetry in given coordinate direction (default is no)
    Vector<int> sym_dir(3,0);
    if (int ns = pp.countval("sym_dir"))
    {
        pp.getarr("sym_dir",sym_dir,0,ns);
    }

    const int idProg = nCompIn;
    const int idSmProg = idProg + 1;
    const int idKm = idSmProg + 1;
    int idSR = -1;
    int nCompOut = 0;
#if BL_SPACEDIM>2
    const int idKg = idKm + 1;
    if (do_strain) {
      idSR = idKg + 1;
      nCompOut = idSR + 1;
    }
    else {
      nCompOut = idKg + 1;
    }
#else
    if (do_strain) {
      idSR = idKm + 1;
      nCompOut = idSR + 1;
    }
    else {
      nCompOut = idKm + 1;
    }
#endif

    bool getStrainTensor = false;
    if (do_strain) pp.query("getStrainTensor",getStrainTensor);
    int idROST=-1;
    if (getStrainTensor)
    {
        idROST = nCompOut;
        nCompOut = idROST + 9; // Rate-of-strain, always return 3D set
    }

    bool getProgGrad = false; pp.query("getProgGrad",getProgGrad);
    int idProgGrad=-1;
    if (getProgGrad)
    {
        idProgGrad = nCompOut;
        nCompOut = idProgGrad + 3; // Always return 3D set
    }

    Vector<MultiFab*> state(Nlev);
    Vector<Geometry*> geoms(Nlev);
    const int nGrow = 2 ;

    // Smoothing
    int num_smooth_pre = 1; pp.query("num_smooth_pre",num_smooth_pre);
    int num_smooth_post = 0; pp.query("num_smooth_post",num_smooth_post);

    Print() << "Will read the following states: ";
    for (int i=0; i<nCompIn; ++i)
        Print() << " " << amrData.StateNumber(inVarNames[i]);
    Print() << '\n';
    Print() << "States out will be those plus: " << '\n';
    Print() << "   idProg:   " << idProg << '\n';
    Print() << "   idSmProg: " << idSmProg << '\n';
    Print() << "   idKm:     " << idKm << '\n';
#if BL_SPACEDIM>2
    Print() << "   idKg:     " << idKg << '\n';
#endif
    if (do_strain) {
      Print() << "   idSR:     " << idSR << '\n';
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

    FArrayBox tmp;
    for (int lev=0; lev<Nlev; ++lev)
    {
        const BoxArray ba = amrData.boxArray(lev);
        DistributionMapping dm(ba);
        state[lev] = new MultiFab(ba,dm,nCompOut,nGrow);
        const Vector<Real>& delta = amrData.DxLevel()[lev];
        Real dxn[3];
        for (int i=0; i<BL_SPACEDIM; ++i) dxn[i] = delta[i];

        // Get input state data onto intersection ba
        const int myNcompIsOne = 1; // gonna need this for fortran calls
        
        Print() << "Reading data for level " << lev << std::endl;

        amrData.FillVar(*state[lev],lev,inVarNames,destFillComps);

        Print() << "Data has been read for level " << lev << std::endl;

        // Build progress variable from state at idCst, put into idCsc
        for (MFIter mfi(*state[lev]); mfi.isValid(); ++mfi)
        {
            FArrayBox& fab = (*state[lev])[mfi];
            const Box& box = mfi.validbox();
            normalize(BL_TO_FORTRAN_BOX(box),
                      BL_TO_FORTRAN_N_ANYD(fab,idCst),
                      BL_TO_FORTRAN_N_ANYD(fab,idProg),
                      &progMin,&progMax,dxn);
        }

        Print() << "Progress variable computed for level " << lev << std::endl;

        geoms[lev] = new Geometry(amrData.ProbDomain()[lev],&rb,coord,&(is_per[0]));
        MultiFab bigMF(BoxArray(ba).grow(nGrow),dm,myNcompIsOne,0); // handy structure to simplify parallel copies

        // Fill unsmoothed progress variable grow cells using interp over c-f boundaries
        {
            MultiFab growScr; // Container in c-f grow cells, filled with interped coarse unsmoothed Progress value
            if (lev>0)
              FillCFgrowCells(amrData.RefRatio()[lev-1],ba,geoms,lev,*state[lev-1],idSmProg,myNcompIsOne,nGrow,growScr);

            MultiFab::Copy(*state[lev],*state[lev],idProg,idSmProg,myNcompIsOne,0); // initialize progress w/unsmoothed

            // Do extrap to fill grow cells.  If not base grid, overwrite c-f grow cells with coarse interp
            for (MFIter mfi(*state[lev]); mfi.isValid(); ++mfi)
            {
/*
                const Box& box = mfi.validbox();
                FORT_PUSHVTOG(box.loVect(),box.hiVect(),
                              *state[lev][mfi].dataPtr(idSmProg),
                              ARLIM(state[lev][mfi].loVect()),ARLIM(state[lev][mfi].hiVect()),
                              &myNcompIsOne,&floorIt);
*/
            }
            
            if (lev > 0)
            {
                // Get data into grow-free mf to allow parallel copy for grabbing c-f data
                for (MFIter mfi(*state[lev]); mfi.isValid(); ++mfi)
                    bigMF[mfi].copy((*state[lev])[mfi],idSmProg,0,myNcompIsOne); // get valid data, and extrap grow data
                
                bigMF.copy(growScr,0,0,myNcompIsOne); // Overwrite c-f with preferred c-f interp data
                
                for (MFIter mfi(bigMF); mfi.isValid(); ++mfi)
                    (*state[lev])[mfi].copy(bigMF[mfi],0,idSmProg,myNcompIsOne); // Put result back into idSmProg
                
            }

            Print() << "Smoothing progress variable on level " << lev << std::endl;

            for (int i=0; i<num_smooth_pre; ++i)
            {
                // Fix up fine-fine and periodic
                state[lev]->FillBoundary(idSmProg,myNcompIsOne,geoms[lev]->periodicity());
            
                // Smooth the data, use slot in state for idKm to hold temporary result
                BL_ASSERT(myNcompIsOne==1); // smooth only knows about a single component
                for (MFIter mfi(*state[lev]); mfi.isValid(); ++mfi)
                {
                    FArrayBox& fab = (*state[lev])[mfi];
                    const Box& box = mfi.validbox();
                    smooth(BL_TO_FORTRAN_BOX(box),
                           BL_TO_FORTRAN_N_ANYD(fab,idSmProg),
                           BL_TO_FORTRAN_N_ANYD(fab,idKm));
                }

                // Set result back into idSmProg
                MultiFab::Copy(*state[lev],*state[lev],idKm,idSmProg,myNcompIsOne,0);
            }
        }
              
        Print() << "Progress variable filled/smoothed on level " << lev << std::endl;
              
        // Fix up fine-fine and periodic for idSmProg
        state[lev]->FillBoundary(idSmProg,myNcompIsOne,geoms[lev]->periodicity());
              
        // Compute curvatures with smoothed progress variable.  Result in state, comp=idKm,idKg
        FArrayBox nWork;
        for (MFIter mfi(*state[lev]); mfi.isValid(); ++mfi)
        {
            FArrayBox& fab = (*state[lev])[mfi];
            const Box& box = mfi.validbox();
            const Box nodebox = amrex::surroundingNodes(box);
            nWork.resize(nodebox,BL_SPACEDIM);
            mcurv(BL_TO_FORTRAN_BOX(box),
                  BL_TO_FORTRAN_N_ANYD(fab,idSmProg),
                  BL_TO_FORTRAN_N_ANYD(fab,idKm),
                  BL_TO_FORTRAN_ANYD(nWork),
                  delta.dataPtr(), sym_dir.dataPtr());
#if BL_SPACEDIM>2
            gcurv(BL_TO_FORTRAN_BOX(box),
                  BL_TO_FORTRAN_N_ANYD(fab,idSmProg),
                  BL_TO_FORTRAN_N_ANYD(fab,idKg),
                  delta.dataPtr());
#endif
        }

        Print() << "Curvature has been computed on level " << lev << std::endl;
        VisMF::Write(*state[lev], "State_Lev"+std::to_string(lev));

//        amrex::Abort("Stop");    

        MultiFab growScr; // Container in c-f grow cells, to be filled with interped coarse Km
        MultiFab scr(ba,dm,myNcompIsOne,nGrow); // Scratch area for smoothing curvature
        MultiFab::Copy(scr,*state[lev],idKm,0,myNcompIsOne,0); // copy curvature into tmp slot for relaxations

        // Do extrap to fill grow cells of scr, if lev>0 improve this with interp from coarse data
        const Box& dbox = amrData.ProbDomain()[lev];
        for (MFIter mfi(scr); mfi.isValid(); ++mfi)
        {
            const Box& box = mfi.validbox();
            pushvtog(BL_TO_FORTRAN_BOX(box),
                     BL_TO_FORTRAN_BOX(dbox),
                     BL_TO_FORTRAN_ANYD(scr[mfi]),
                     &myNcompIsOne);
        }
        // Fix up fine-fine and periodic for scr
        scr.FillBoundary(0,myNcompIsOne,geoms[lev]->periodicity());

        if (lev > 0 && num_smooth_post>0)
        {
            FillCFgrowCells(amrData.RefRatio()[lev-1],ba,geoms,lev,*state[lev-1],idKm,myNcompIsOne,nGrow,growScr);

            // Get data into grow-free mf to allow parallel copy for grabbing c-f data
            for (MFIter mfi(scr); mfi.isValid(); ++mfi)
                bigMF[mfi].copy(scr[mfi],0,0,myNcompIsOne); // get valid data, and extrap grow data
            
            bigMF.copy(growScr,0,0,myNcompIsOne); // Overwrite c-f with preferred c-f interp data
            
            for (MFIter mfi(bigMF); mfi.isValid(); ++mfi)
                scr[mfi].copy(bigMF[mfi],0,0,myNcompIsOne); // Put result back into scr
        }

        // Now that all f-f grow cells filled, do smoothing
        Print() << "Smoothing curvature on level " << lev << std::endl;

        for (int i=0; i<num_smooth_post; ++i)
        {
            // Fix up fine-fine and periodic
            scr.FillBoundary(0,myNcompIsOne,geoms[lev]->periodicity());
            
            // Smooth the data, use idKm to hold result
            BL_ASSERT(myNcompIsOne==1); //  only knows about a single component
            for (MFIter mfi(*state[lev]); mfi.isValid(); ++mfi)
            {
                const Box& box = mfi.validbox();
                smooth(BL_TO_FORTRAN_BOX(box),
                       BL_TO_FORTRAN_N_ANYD(scr[mfi],0),
                       BL_TO_FORTRAN_N_ANYD((*state[lev])[mfi],idKm));
            }

            // Set result back into 0 component of scr for next iteration
            MultiFab::Copy(scr,*state[lev],idKm,0,myNcompIsOne,0);
        }

        if (do_strain)
        {
          Print() << "Filling velocity grow cells on level " << lev << std::endl;

          // Fill boundary cells for velocity
          // Do extrap to fill grow cells of scr, if lev>0 improve this with interp from coarse data
          const int nCompVEL=BL_SPACEDIM;
          const int nCompC=1;
          const int Ccomp = idSmProg;
          for (MFIter mfi(*state[lev]); mfi.isValid(); ++mfi)
          {
            const Box& box = mfi.validbox();
            pushvtog(BL_TO_FORTRAN_BOX(box),
                     BL_TO_FORTRAN_BOX(dbox),
                     BL_TO_FORTRAN_N_ANYD((*state[lev])[mfi],idVst),
                     &nCompVEL);
            pushvtog(BL_TO_FORTRAN_BOX(box),
                     BL_TO_FORTRAN_BOX(dbox),
                     BL_TO_FORTRAN_N_ANYD((*state[lev])[mfi],Ccomp),
                     &nCompC);
          }
          // Fix up fine-fine and periodic for scr
          state[lev]->FillBoundary(idVst,nCompVEL,geoms[lev]->periodicity());
          state[lev]->FillBoundary(Ccomp,nCompC,geoms[lev]->periodicity());

          Print() << "Computing tangential strain on level " << lev << std::endl;

          FArrayBox nWork2;
          for (MFIter mfi(*state[lev]); mfi.isValid(); ++mfi)
          {
            const Box& box = mfi.validbox();
            const Box nodebox = amrex::surroundingNodes(box);
            nWork2.resize(nodebox,BL_SPACEDIM);
            FArrayBox& s = (*state[lev])[mfi];
            strainrate(BL_TO_FORTRAN_BOX(box),
                       BL_TO_FORTRAN_N_ANYD(s,idVst),
                       BL_TO_FORTRAN_N_ANYD(s,Ccomp),
                       BL_TO_FORTRAN_N_ANYD(s,idSR),
                       BL_TO_FORTRAN_ANYD(nWork2),
                       delta.dataPtr());
          }


          if (getStrainTensor)
          {
            Print() << "Filling rate of strain on level " << lev << std::endl;

            const int nCompV=BL_SPACEDIM;

            // Do extrap to fill grow cells of scr, if lev>0 improve this with interp from coarse data
            for (MFIter mfi(*state[lev]); mfi.isValid(); ++mfi)
            {
                FArrayBox& fab = (*state[lev])[mfi];
                const Box& box = mfi.validbox();
                pushvtog(BL_TO_FORTRAN_BOX(box),
                         BL_TO_FORTRAN_BOX(dbox),
                         BL_TO_FORTRAN_N_ANYD(fab,idVst),
                         &nCompV);
            }

            // Fix up fine-fine and periodic for scr
            state[lev]->FillBoundary(idVst,nCompV,geoms[lev]->periodicity());

            if (lev > 0)
            {
                FillCFgrowCells(amrData.RefRatio()[lev-1],ba,geoms,lev,*state[lev-1],idVst,nCompV,nGrow,growScr);

                // Get data into grow-free mf to allow parallel copy for grabbing c-f data
                MultiFab bigMF1(BoxArray(ba).grow(nGrow),dm,nCompV,0);

                for (MFIter mfi(scr); mfi.isValid(); ++mfi)
                    bigMF1[mfi].copy((*state[lev])[mfi],idVst,0,nCompV); // get valid data, and extrap grow data
                
                bigMF1.copy(growScr,0,0,nCompV); // Overwrite c-f with preferred c-f interp data

                for (MFIter mfi(bigMF1); mfi.isValid(); ++mfi)
                    (*state[lev])[mfi].copy(bigMF1[mfi],0,idVst,nCompV); // Put result back into scr
            }

            for (MFIter mfi(*state[lev]); mfi.isValid(); ++mfi)
            {
                const Box& box = mfi.validbox();
                FArrayBox& s = (*state[lev])[mfi];
                straintensor(BL_TO_FORTRAN_BOX(box),
                             BL_TO_FORTRAN_N_ANYD(s,idVst),
                             BL_TO_FORTRAN_N_ANYD(s,idROST),
                             delta.dataPtr());
            }
          }
        }

        //if (getProgGrad)
        //{
        //    Print() << "Computing gradient of progress on level " << lev << std::endl;
        //    
        //    const int nCompCG=1;

        //    // Do extrap to fill grow cells of scr, if lev>0 improve this with interp from coarse data
        //    for (MFIter mfi(*state[lev]); mfi.isValid(); ++mfi)
        //    {
        //        FArrayBox& fab = (*state[lev])[mfi];
        //        const Box& box = mfi.validbox();
        //        pushvtog(BL_TO_FORTRAN_BOX(box),
        //                 BL_TO_FORTRAN_BOX(dbox),
        //                 BL_TO_FORTRAN_N_ANYD(fab,idCst),
        //                 &nCompCG);
        //    }

        //    // Fix up fine-fine and periodic for scr
        //    state[lev]->FillBoundary(idCst,nCompCG,geoms[lev]->periodicity());

        //    if (lev > 0)
        //    {
        //        FillCFgrowCells(amrData.RefRatio()[lev-1],ba,geoms,lev,*state[lev-1],idCst,nCompCG,nGrow,growScr);

        //        // Get data into grow-free mf to allow parallel copy for grabbing c-f data
        //        MultiFab bigMFG(BoxArray(ba).grow(nGrow),dm,nCompCG,0);
        //        for (MFIter mfi(scr); mfi.isValid(); ++mfi)
        //            bigMFG[mfi].copy((*state[lev])[mfi],idCst,0,nCompCG); // get valid data, and extrap grow data

        //        const int nCompC=1;
        //        bigMFG.copy(growScr,0,0,nCompC); // Overwrite c-f with preferred c-f interp data
        //    
        //        for (MFIter mfi(bigMFG); mfi.isValid(); ++mfi)
        //            (*state[lev])[mfi].copy(bigMFG[mfi],0,idCst,nCompCG); // Put result back into scr
        //    }

        //    for (MFIter mfi(*state[lev]); mfi.isValid(); ++mfi)
        //    {
        //        const Box& box = mfi.validbox();
        //        FArrayBox& s = (*state[lev])[mfi];
        //        progressgrad(BL_TO_FORTRAN_BOX(box),
        //                     BL_TO_FORTRAN_N_ANYD(s,idCst),
        //                     BL_TO_FORTRAN_N_ANYD(s,idProgGrad),
        //                     delta.dataPtr());
        //    }
        //}

        //Print() << "Derive finished for level " << lev << std::endl;
    }

    Vector<std::string> nnames(nCompOut);
    for (int i=0; i<nCompIn; ++i)
        nnames[i] = inVarNames[i];
    nnames[idProg] = "Progress";
    nnames[idSmProg] = "SmoothedProgress";
    nnames[idKm] = "MeanCurvature_" + progressName;
#if BL_SPACEDIM>2
    nnames[idKg] = "GaussianCurvature_" + progressName;
#endif
    if (do_strain) nnames[idSR] = "StrainRate_" + progressName;

    if (getStrainTensor)
    {
        for (int i=0; i<9; ++i)
        {
            int n=i%3;
            int m=i/3;
            char buf[40];
            sprintf(buf,"ROST_%1d%1d",n+1,m+1);
            nnames[idROST+i] = string(buf);
        }
    }
    
    if (getProgGrad)
    {
        nnames[idProgGrad+0] = progressName + "_g1";
        nnames[idProgGrad+1] = progressName + "_g2";
        nnames[idProgGrad+2] = progressName + "_g3";
    }

    bool verb=true;
    bool appendPlotFile=false; pp.query("appendPlotFile",appendPlotFile);
    if (appendPlotFile)
    {

        int nStateOut = nCompOut - nCompIn;
        Vector<MultiFab*> ostate(Nlev);
        for (int lev=0; lev<Nlev; ++lev)
        {
            const BoxArray ba = state[lev]->boxArray();
            ostate[lev] = new MultiFab(ba,DistributionMapping(ba),nStateOut,0);
            MultiFab::Copy(*ostate[lev],*state[lev],nCompIn,0,nStateOut,0);
        }
        Vector<std::string> namesOut(nStateOut);
        for (int i=0; i<nStateOut; ++i)
            namesOut[i] = nnames[nCompIn+i];

        std::string newMFBaseName = "NEWDAT"; pp.query("newMFBaseName",newMFBaseName);
        std::string newHeaderName = "NewHeader"; pp.query("newHeaderName",newHeaderName);
        AppendToPlotFile(amrData,ostate,plotFileName,namesOut,newMFBaseName,newHeaderName,verb);
        Print() << "...finished.  Note: to see new data, you must rename NewHeader in the" << std::endl;
        Print() << "              pltfile to Header (probably want to save the original first)" << std::endl;
    }
    else
    {
        // Remove GC from outstate
        Vector<MultiFab*> ostate(Nlev);
        for (int lev=0; lev<Nlev; ++lev)
        {
            const BoxArray ba = state[lev]->boxArray();
            ostate[lev] = new MultiFab(ba,DistributionMapping(ba),nCompOut,0);
            MultiFab::Copy(*ostate[lev],*state[lev],0,0,nCompOut,0);
        }
        Print() << "Writing new data to " << outfile << std::endl;
        WritePlotFile(ostate,amrData,outfile,verb,nnames);
    }
    }
    amrex::Finalize();
    return 0;
}
