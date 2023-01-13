#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_DataServices.H>
#include <AMReX_WritePlotFile.H>
#include <AMReX_PlotFileUtil.H>
#include <PltFileManager.H>

using namespace amrex;

struct pltfileMesh {
   void init (int a_nlevs) {
      nLevels = a_nlevs;
      BAVec.resize(nLevels);
   };
   int nLevels{0};
   Vector<BoxArray> BAVec;
};

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
    for (auto &file : plotFileNames) {
       if (!amrex::FileSystem::Exists(file)) {
          Abort(file+" does not exists");
       }
    }

    DataServices::SetBatchMode();
    Amrvis::FileType fileType(Amrvis::NEWPLT);

    // Build a Vector of BoxList where each entry is the union of the pltfile BoxList
    // at each level
    int finestPltLevel = 0;
    Vector<BoxList> avgBL(1);
    Vector<pltfileMesh> meshes(nf);
    Box fixedCoarseDomain;
    Vector<Geometry> geoms;
    for (int f{0}; f<nf; ++f) {
        DataServices dataServices(plotFileNames[f], fileType);
        if( ! dataServices.AmrDataOk()) {
            DataServices::Dispatch(DataServices::ExitRequest, NULL);
        }
        AmrData& amrData = dataServices.AmrDataRef();
        if (f==0) {
          fixedCoarseDomain = amrData.ProbDomain()[0];
        } else {
          AMREX_ALWAYS_ASSERT(fixedCoarseDomain == amrData.ProbDomain()[0]);
        }
        // Stash away the mesh data of this pltfile
        if (amrData.FinestLevel() > finestPltLevel) {
            finestPltLevel = amrData.FinestLevel();
        }
        meshes[f].init(amrData.FinestLevel()+1);
        for (int lev{0}; lev<meshes[f].nLevels; ++lev) {
            meshes[f].BAVec[lev] = amrData.boxArray(lev);
        }
        // Append to the union of BoxList of each level
        if (f==0) {
            avgBL.resize(meshes[f].nLevels);
            for (int lev{0}; lev<meshes[f].nLevels; ++lev) {
                avgBL[lev] = meshes[f].BAVec[lev].boxList();
            }
            geoms.resize(meshes[f].nLevels);
            RealBox rb(AMREX_D_DECL(amrData.ProbLo()[0],amrData.ProbLo()[1],amrData.ProbLo()[2]),
                       AMREX_D_DECL(amrData.ProbHi()[0],amrData.ProbHi()[1],amrData.ProbHi()[2]));
            Array<int,AMREX_SPACEDIM> is_per = {0};
            for (int lev{0}; lev<meshes[f].nLevels; ++lev) {
                geoms[lev].define(amrData.ProbDomain()[lev],rb,amrData.CoordSys(),is_per);
            }
        } else {
            for (int lev{0}; lev<meshes[f].nLevels; ++lev) {
               if (lev >= avgBL.size()) {
                   avgBL.push_back(meshes[f].BAVec[lev].boxList());
               } else {
                   auto currAvgBA = BoxArray(avgBL[lev]);
                   auto baCompDom = complementIn(amrData.ProbDomain()[lev], currAvgBA);
                   auto baAdd = intersect(baCompDom,meshes[f].BAVec[lev]);
                   for (int bx{0}; bx<baAdd.size(); ++bx) {
                      avgBL[lev].push_back(baAdd[bx]);
                   }
               }
            }
        }
    }

    /*
    Vector<MultiFab> MFtest(finestPltLevel+1);
    for (int lev{0}; lev<finestPltLevel+1; ++lev) {
        auto finalAvgBA = BoxArray(avgBL[lev]);
        auto dmap = DistributionMapping(finalAvgBA);
        MFtest[lev].define(finalAvgBA,dmap,1,0);
        MFtest[lev].setVal(lev+1);
    }
    Vector<int> isteps(finestPltLevel+1, 0);
    Vector<IntVect> refRatios(finestPltLevel,{AMREX_D_DECL(2, 2, 2)});
    amrex::WriteMultiLevelPlotfile("Test", finestPltLevel+1, GetVecOfConstPtrs(MFtest), {"test"},
                                   geoms, 0.0, isteps, refRatios);
    */

    std::string outfile("PltAvg");
    pp.query("outfile",outfile);

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
      for (int i=0; i<nComp; ++i) {
        comps[i] = sComp + i;
      }
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
      Vector<MultiFab> MFweight(finestPltLevel+1);
      Vector<MultiFab> MFout(finestPltLevel+1);
      for (int lev{0}; lev<finestPltLevel+1; ++lev) {
          auto finalAvgBA = BoxArray(avgBL[lev]);
          auto dmap = DistributionMapping(finalAvgBA);
          MFweight[lev].define(finalAvgBA,dmap,1,0);
          MFweight[lev].setVal(0.0);
          MFout[lev].define(finalAvgBA,dmap,comps.size(),0);
          MFout[lev].setVal(0.0);
      }

      for (int f{0}; f<plotFileNames.size(); ++f) {
        Print() << " Treating " << plotFileNames[f] << "\n";
        pele::physics::pltfilemanager::PltFileManager pltData(plotFileNames[f]);
        Vector<std::string> plt_vars = pltData.getVariableList();
        int nvars = plt_vars.size();

        for (int lev=0; lev<pltData.getNlev(); ++lev) {
            MultiFab temp(MFout[lev].boxArray(),
                          MFout[lev].DistributionMap(),
                          MFout[lev].nComp(),0);
            pltData.fillPatchFromPlt(0, geoms[lev], 0, 0, nvars, temp);
            MultiFab::Add(MFout[lev],temp,0,0,nvars,0);
        }
      }

      for (int lev{0}; lev<finestPltLevel+1; ++lev) {
          MFout[lev].mult(1.0/nf);
      }
      Vector<int> isteps(finestPltLevel+1, 0);
      Vector<IntVect> refRatios(finestPltLevel,{AMREX_D_DECL(2, 2, 2)});
      amrex::WriteMultiLevelPlotfile(outfile, finestPltLevel+1, GetVecOfConstPtrs(MFout), pltnames,
                                     geoms, 0.0, isteps, refRatios);
      }
  }
  Finalize();
  return 0;
}

