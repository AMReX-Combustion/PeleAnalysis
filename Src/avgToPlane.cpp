#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include <AMReX_DataServices.H>
#include <WritePlotFile.H>

using namespace amrex;

static
void 
print_usage (int,
             char* argv[])
{
  std::cerr << "usage:\n";
  std::cerr << argv[0] << " infile=f1 [options] \n\tOptions:\n";
  exit(1);
}

std::string
getFileRoot(const std::string& infile)
{
  vector<std::string> tokens = Tokenize(infile,std::string("/"));
  return tokens[tokens.size()-1];
}

Box ProjectBox(const Box& srcBox,
               int        dir,
               int        loc)
{
  Box dstBox = srcBox;
  dstBox.setSmall(dir,loc);
  dstBox.setBig(dir,loc);
  return dstBox;
}

int
main (int   argc,
      char* argv[])
{
  Initialize(argc,argv);
  {
    if (argc < 2)
      print_usage(argc,argv);

    ParmParse pp;

    if (pp.contains("help"))
      print_usage(argc,argv);

    // Open plotfile header and create an amrData object pointing into it
    std::string plotFileName; pp.get("infile",plotFileName);
    DataServices::SetBatchMode();
    Amrvis::FileType fileType(Amrvis::NEWPLT);
    DataServices dataServices(plotFileName, fileType);
    if( ! dataServices.AmrDataOk()) {
      DataServices::Dispatch(DataServices::ExitRequest, NULL);
    }
    AmrData& amrData = dataServices.AmrDataRef();

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
      int nComp = amrData.NComp();
      pp.query("nComp",nComp);
      BL_ASSERT(sComp+nComp <= amrData.NComp());
      comps.resize(nComp);
      for (int i=0; i<nComp; ++i)
        comps[i] = sComp + i;
    }

    int finestLevel = amrData.FinestLevel(); pp.query("finestLevel",finestLevel);
    AMREX_ALWAYS_ASSERT(finestLevel >= 0 && finestLevel<=amrData.FinestLevel());
    Box subbox = amrData.ProbDomain()[finestLevel];

    if (int nx=pp.countval("box"))
    {
      Vector<int> inBox;
      pp.getarr("box",inBox,0,nx);
      int d=BL_SPACEDIM;
      BL_ASSERT(inBox.size()==2*d);
      subbox=Box(IntVect(D_DECL(inBox[0],inBox[1],inBox[2])),
                 IntVect(D_DECL(inBox[d],inBox[d+1],inBox[d+2])),
                 IndexType::TheCellType());
    }

    Vector<std::string> names(comps.size());
   
    Vector<Real> plo(BL_SPACEDIM), phi(BL_SPACEDIM);
    Vector<Box> psize(finestLevel+1);
    const IntVect ilo = subbox.smallEnd();
    const IntVect ihi = subbox.bigEnd();

    for (int i =0 ; i< BL_SPACEDIM; i++) {   
       plo[i] = amrData.ProbLo()[i]+(ilo[i])*amrData.DxLevel()[finestLevel][i];
       phi[i] = amrData.ProbLo()[i]+(ihi[i]+1)*amrData.DxLevel()[finestLevel][i];
    }

    int dir = 0; pp.query("dir",dir); AMREX_ALWAYS_ASSERT(dir>=0 && dir<AMREX_SPACEDIM);
    int loc = 0;
    
    BoxArray ba_sub(subbox);
    int max_grid_size = 32; pp.query("max_grid_size",max_grid_size);
    ba_sub.maxSize(max_grid_size);

    if (ba_sub.size() > 0)
    {
      BoxList bl_flat;
      for (int i=0; i<ba_sub.size(); ++i) {
        bl_flat.push_back(ProjectBox(ba_sub[i],dir,loc));
      }
      BoxArray ba_flat(bl_flat);
      DistributionMapping dmap_sub(ba_sub);

      MultiFab mf_full(ba_sub,dmap_sub,1,0);
      MultiFab mf_flat(ba_flat,dmap_sub,1,0); mf_flat.setVal(0);

      Box res_box = ProjectBox(subbox,dir,loc);
      const Vector<int> dm_vec({0});
      DistributionMapping res_dm(dm_vec);
      MultiFab res_mf(BoxArray(res_box),res_dm,comps.size(),0); res_mf.setVal(0);

      for (int i=0; i<comps.size(); ++i)
      {
        names[i] = amrData.PlotVarNames()[comps[i]];
        if (ParallelDescriptor::IOProcessor()) 
          std::cout << "Filling " << names[i] << " on level " << finestLevel << std::endl;

        mf_full.ParallelCopy(amrData.GetGrids(finestLevel,comps[i],subbox),0,0,1);
        amrData.FlushGrids(comps[i]);                

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(mf_full, false); mfi.isValid(); ++mfi) // do not tile, to avoid race on += below
        {
          auto const& full_fab = mf_full[mfi];
          Array4<Real const> const& full_arr = full_fab.const_array();

          const Box& tile_box  = mfi.tilebox();
          Box flat_box = ProjectBox(tile_box,dir,loc);
          auto& flat_fab = mf_flat[mfi];
          Array4<Real> const& flat_arr = flat_fab.array();

          AMREX_LAUNCH_HOST_DEVICE_FUSIBLE_LAMBDA ( tile_box, thread_box,
          {
            const auto lo = lbound(thread_box);
            const auto hi = ubound(thread_box);

            if (dir==0) {
              int ilo = flat_box.smallEnd(dir);
              for     (int k = lo.z; k <= hi.z; ++k) {
                for   (int j = lo.y; j <= hi.y; ++j) {
                    for (int i = lo.x; i <= hi.x; ++i) {
                      flat_arr(ilo,j,k) += full_arr(i,j,k);
                    }
                }
              }
            }
            else if (dir==1) {
              int jlo = flat_box.smallEnd(dir);
              for     (int k = lo.z; k <= hi.z; ++k) {
                for   (int j = lo.y; j <= hi.y; ++j) {
                    for (int i = lo.x; i <= hi.x; ++i) {
                      flat_arr(i,jlo,k) += full_arr(i,j,k);
                    }
                }
              }
            }
            else if (dir==2) {
              int klo = flat_box.smallEnd(dir);
              for     (int k = lo.z; k <= hi.z; ++k) {
                for   (int j = lo.y; j <= hi.y; ++j) {
                    for (int i = lo.x; i <= hi.x; ++i) {
                      flat_arr(i,j,klo) += full_arr(i,j,k);
                    }
                }
              }
            }
          });
        }
        res_mf.ParallelAdd(mf_flat,0,i,1);
      }

      std::string outfile=getFileRoot(plotFileName) + "_avg.fab";
      pp.query("outfile",outfile);
      if (ParallelDescriptor::IOProcessor()) {
        std::ofstream ofs(outfile.c_str());
        res_mf[0].writeOn(ofs);
        ofs.close();
      }
    }
  }
  Finalize();
  return 0;
}
