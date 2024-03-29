#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_DataServices.H>
#include <AMReX_WritePlotFile.H>
#include <AMReX_PlotFileUtil.H>

using namespace amrex;

int
main (int   argc,
      char* argv[])
{
  Initialize(argc,argv);
  {
    ParmParse pp;

    // Open first plotfile header and create an amrData object pointing into it
    int nf = pp.countval("infiles");
    AMREX_ALWAYS_ASSERT(nf>0);
    Vector<std::string> plotFileNames; pp.getarr("infiles",plotFileNames,0,nf);

    std::string outfile("JUNK");
    pp.query("outfile",outfile);

    DataServices::SetBatchMode();
    Amrvis::FileType fileType(Amrvis::NEWPLT);
    DataServices dataServices0(plotFileNames[0], fileType);
    if( ! dataServices0.AmrDataOk()) {
      DataServices::Dispatch(DataServices::ExitRequest, NULL);
    }
    AmrData& amrData0 = dataServices0.AmrDataRef();

    Vector<int> comps;
    if (int nc = pp.countval("comps"))
    {
      comps.resize(nc);
      pp.getarr("comps",comps,0,nc);
    }
    else
    {
      int sComp = 0;
      pp.query("sComp",sComp);
      int nComp = amrData0.NComp();
      pp.query("nComp",nComp);
      AMREX_ASSERT(sComp+nComp <= amrData0.NComp());
      comps.resize(nComp);
      for (int i=0; i<nComp; ++i)
        comps[i] = sComp + i;
    }
    Vector<std::string> pltnames(comps.size());
    for (int i=0; i<comps.size(); ++i) {
      pltnames[i] = amrData0.PlotVarNames()[comps[i]];
    }
    int finestLevel = amrData0.FinestLevel(); pp.query("finestLevel",finestLevel);
    AMREX_ALWAYS_ASSERT(finestLevel >= 0 && finestLevel<=amrData0.FinestLevel());

    Box subbox = amrData0.ProbDomain()[finestLevel];
    /* Deactivate this for now. Lead to small differences to the original
     * pltfile box due to arithmetic operations
    if (int nx=pp.countval("box"))
    {
      Vector<int> inBox;
      pp.getarr("box",inBox,0,nx);
      int d=AMREX_SPACEDIM;
      AMREX_ASSERT(inBox.size()==2*d);
      subbox=Box(IntVect(D_DECL(inBox[0],inBox[1],inBox[2])),
                 IntVect(D_DECL(inBox[d],inBox[d+1],inBox[d+2])),
                 IndexType::TheCellType());
    }

    Vector<Real> plo(AMREX_SPACEDIM), phi(AMREX_SPACEDIM);
    const IntVect ilo = subbox.smallEnd();
    const IntVect ihi = subbox.bigEnd();

    for (int i =0 ; i< AMREX_SPACEDIM; i++) {
       plo[i] = amrData0.ProbLo()[i]+(ilo[i])*amrData0.DxLevel()[finestLevel][i];
       phi[i] = amrData0.ProbLo()[i]+(ihi[i]+1)*amrData0.DxLevel()[finestLevel][i];
    }
    */

    const Vector<int>& ratio = amrData0.RefRatio();

    Vector<Box> domain(finestLevel+1);
    Vector<int> ratioTot(finestLevel+1);
    domain[finestLevel] = subbox;
    for (int lev=finestLevel-1; lev>=0; --lev) {
      domain[lev] = coarsen(domain[lev+1],ratio[lev]);
    }
    for (int lev=1; lev<=finestLevel; ++lev) {
      domain[lev] = refine(domain[lev-1],ratio[lev-1]);
    }
    for (int lev=0; lev<=finestLevel; ++lev) {
      ratioTot[lev] = (lev == 0 ? 1 : ratioTot[lev-1] * ratio[lev-1]);
    }

    AMREX_ALWAYS_ASSERT(domain[0].numPts() > 0);
    Vector<Geometry> geom(finestLevel+1);
    RealBox rb(AMREX_D_DECL(amrData0.ProbLo()[0],amrData0.ProbLo()[1],amrData0.ProbLo()[2]),
               AMREX_D_DECL(amrData0.ProbHi()[0],amrData0.ProbHi()[1],amrData0.ProbHi()[2]));
    Array<int,AMREX_SPACEDIM> is_per = {0};
    for (int lev=0; lev<=finestLevel; ++lev) {
      geom[lev].define(domain[lev],rb,amrData0.CoordSys(),is_per);
    }
    BoxArray ba_res(domain[0]);
    int max_grid_size = 32; pp.query("max_grid_size",max_grid_size);
    ba_res.maxSize(max_grid_size);

    bool all_same_boxes = false; pp.query("all_same_boxes",all_same_boxes);
    Print() << "All same boxes: " << all_same_boxes << std::endl;

    if (all_same_boxes) {
      Vector<MultiFab*> mf_avg(finestLevel+1);
      for (int lev=0; lev<=finestLevel; ++lev) {
        const BoxArray& ba_lev = amrData0.boxArray(lev);
        mf_avg[lev] = new MultiFab(ba_lev,DistributionMapping(ba_lev),comps.size(),0);
        mf_avg[lev]->setVal(0);
      }

      for (int iFile=0; iFile<plotFileNames.size(); ++iFile) {
        DataServices dataServices(plotFileNames[iFile], fileType);
        if( ! dataServices.AmrDataOk())
        {
          DataServices::Dispatch(DataServices::ExitRequest, NULL);
        }
        AmrData& amrData = dataServices.AmrDataRef();

        Print() << "Reading data in " << plotFileNames[iFile] << std::endl;

        for (int lev=0; lev<=finestLevel; ++lev)
        {
          AMREX_ALWAYS_ASSERT_WITH_MESSAGE(amrData.boxArray(lev) == mf_avg[lev]->boxArray(),"BoxArrays not equal");

          for (int i=0; i<comps.size(); ++i)
          {
            MultiFab mfComp(mf_avg[lev]->boxArray(),mf_avg[lev]->DistributionMap(),1,0);
            mfComp.setVal(0);
            mfComp.ParallelCopy(amrData.GetGrids(lev,comps[i]));
            MultiFab::Add(*mf_avg[lev],mfComp,0,i,1,0);
          }
        }
      }
      Print() << "Writing output to " << outfile << std::endl;
      for (int lev=0; lev<=finestLevel; ++lev) {
        mf_avg[lev]->mult(1.0/nf);
      }
      Vector<int> isteps(finestLevel+1, 0);
      Vector<IntVect> refRatios(finestLevel,{AMREX_D_DECL(2, 2, 2)});
      amrex::WriteMultiLevelPlotfile(outfile, finestLevel+1, GetVecOfConstPtrs(mf_avg), pltnames,
                                     geom, 0.0, isteps, refRatios);
    }
    else
    {
      MultiFab mf_avg(ba_res,DistributionMapping(ba_res),comps.size(),0);
      mf_avg.setVal(0);

      MultiFab mf_avgDn(ba_res,DistributionMapping(ba_res),comps.size(),0);
      for (int iFile=0; iFile<plotFileNames.size(); ++iFile) {
        DataServices dataServices(plotFileNames[iFile], fileType);
        if( ! dataServices.AmrDataOk()) {
          DataServices::Dispatch(DataServices::ExitRequest, NULL);
        }
        AmrData& amrData = dataServices.AmrDataRef();

        Print() << "Reading data in " << plotFileNames[iFile] << std::endl;

        for (int lev=0; lev<=finestLevel; ++lev) {
          for (int i=0; i<comps.size(); ++i)
          {
            MultiFab alias(mf_avgDn,make_alias,i,1);
            if (lev==0) {
              alias.ParallelCopy(amrData.GetGrids(0,comps[i],domain[0]));
            } else {
              average_down(amrData.GetGrids(lev,comps[i],domain[lev]),alias,geom[lev],geom[0],0,1,ratioTot[lev]);
            }
            amrData.FlushGrids(comps[i]);
          }
        }

        MultiFab::Add(mf_avg,mf_avgDn,0,0,comps.size(),0);
      }

      mf_avg.mult(1.0/nf);
      writePlotFile(outfile.c_str(),mf_avg,geom[0],IntVect(AMREX_D_DECL(2,2,2)),0.0,pltnames);
    }
  }
  Finalize();
  return 0;
}

