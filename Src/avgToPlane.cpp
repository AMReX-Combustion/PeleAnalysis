#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include <AMReX_DataServices.H>
#include <AMReX_WritePlotFile.H>

using namespace amrex;

void STORE_PPM_STR (const std::string& file,
                    int width, int height, int iimage[],
                    const int r[], const int g[], const int b[]);
void STORE_PGM_STR (const std::string& file,
                    int width, int height, int iimage[]);
void LOAD_PALETTE_STR (const std::string& file, int r[], int g[], int b[], int a[]);
void pixelizeData(const FArrayBox& data, int slicedir, int sliceloc,
                  BaseFab<int>& image, Real data_min, Real data_max, int nVals);

static std::string AND("_");
static std::string PPM(".ppm");
static std::string PGM(".pgm");
static std::string FABF(".fab");

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
  std::vector<std::string> tokens = Tokenize(infile,std::string("/"));
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
      subbox=Box(IntVect(AMREX_D_DECL(inBox[0],inBox[1],inBox[2])),
                 IntVect(AMREX_D_DECL(inBox[d],inBox[d+1],inBox[d+2])),
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

      if (ParallelDescriptor::IOProcessor()) {
        FArrayBox& data = res_mf[0];
        Real data_min = data.min();
        Real data_max = data.max();
        pp.query("min",data_min);
        pp.query("max",data_max);
      
        BaseFab<int> image;
        const int nVals = 256;
        pixelizeData(data,dir,loc,image,data_min,data_max,nVals);
      
        const int width = image.box().length(0);
        const int height= image.box().length(1);

        std::string palfile;
        pp.get("palette",palfile);
        int r[256], g[256], b[256], t[256];
        LOAD_PALETTE_STR(palfile,r,g,b,t);
            
        std::string outfile = getFileRoot(plotFileName) + PPM;
        STORE_PPM_STR(outfile, width, height, image.dataPtr(), r, g, b);
      }
      
    }
  }
  Finalize();
  return 0;
}

void pixelizeData(const FArrayBox& data, int slicedir, int sliceloc,
                  BaseFab<int>& image, Real data_min, Real data_max, int nVals)
{
  const Box& box = data.box();
  const Real del = data_max - data_min;
  BL_ASSERT(std::abs(del)>0.0);
  const int nvm1 = nVals-1;
  int cnt=0;
  int d[2];
  for (int dir=0; dir<BL_SPACEDIM; ++dir)
    if (dir != slicedir)
      d[cnt++] = dir;

  const IntVect se = box.smallEnd();
  const IntVect be = box.bigEnd();

  Print() << Box(se,be) << std::endl;

  IntVect img(AMREX_D_DECL(box.length(d[0]) - 1,box.length(d[1]) - 1,0));
  image.resize(Box(IntVect::TheZeroVector(),img),1);

  IntVect div;
  if (slicedir>=0 && slicedir<BL_SPACEDIM)
      div[slicedir] = sliceloc;

  for (int i=se[d[0]]; i<=be[d[0]]; ++i) {
    for (int j=se[d[1]]; j<=be[d[1]]; ++j) {
      div[d[0]] = i;
      div[d[1]] = j;
      image(IntVect(AMREX_D_DECL(i - se[d[0]],j - se[d[1]],0)),0) =
        std::max(0,(int)(nvm1*std::min( (data(div,0) - data_min)/del,1.0))); 
    }
  }
  //BaseFab<int> imageRev(image.box(), 1);
  //imageRev.copy(image);
  //int rMult(1);
  //image.copyRev(image.box(), imageRev, image.box(), 1, &rMult);
}

#include <cstdlib>
#include <cstdio>

static const char PGM_MAGIC1 = 'P';
static const char RPGM_MAGIC2 = '5';
static const char RPGM_MAGIC3 = '6';

void
STORE_PPM_STR (const std::string& file, int width, int height, int iimage[],
	       const int r[], const int g[], const int b[])
{
  FILE *image_fp;	/* file descriptor for image (output) */
  
  /* create image file */
  if ((image_fp = fopen(file.c_str(), "w")) == NULL)
  {
      fprintf(stderr, "cannot open output file %s\n", file.c_str());
      exit(1);
  }
  
  /* translate Fortran image data to chars */
  unsigned char* image = new unsigned char[3*width*height];
  for (int i = 0; i < width*height; i++ )
  {
      const int j = iimage[i];
      if ( j < 0 || j > 255 ) 
      {
	  fprintf(stderr,"out of bounds on image[%d] = %d\n", i, j);
	  exit(1);
      }
      image[3*i+0] = r[j];
      image[3*i+1] = g[j];
      image[3*i+2] = b[j];
  }
  fprintf(image_fp, "%c%c\n%d %d\n%d\n", PGM_MAGIC1, RPGM_MAGIC3,
	  width, height, 255);
  fwrite(image, 1, 3*width*height, image_fp);
  delete [] image;
  fclose(image_fp);
}

void
STORE_PGM_STR (const std::string& file, int width, int height, int iimage[])
{
  FILE *image_fp;	/* file descriptor for image (output) */

  /* create image file */
  if ((image_fp = fopen(file.c_str(), "w")) == NULL)
    {
      fprintf(stderr, "cannot open output file %s\n", file.c_str());
      exit(1);
    }

  /* translate Fortran image data to chars */
  unsigned char* image = new unsigned char[width*height];
  for (int i = 0; i < width*height; ++i )
  {
      image[i] = iimage[i];
  }
  fprintf(image_fp, "%c%c\n%d %d\n%d\n", PGM_MAGIC1, RPGM_MAGIC2,
	  width, height, 255);
  fwrite(image, 1, width*height, image_fp);
  delete [] image;
  fclose(image_fp);
}

void
LOAD_PALETTE_STR (const std::string& file, int r[], int g[], int b[], int a[])
{
  FILE *pal_fp;	/* file descriptor for image (output) */
  unsigned char c[256];
  int i;

  if ((pal_fp = fopen(file.c_str(), "rb")) == NULL)
    {
      fprintf(stderr, "cannot open palette file %s\n", file.c_str());
      exit(1);
    }
  fread(c, 1, 256, pal_fp);
  for ( i = 0; i < 256; ++i )
    {
      r[i] = c[i];
    }
  fread(c, 1, 256, pal_fp);
  for ( i = 0; i < 256; ++i )
    {
      g[i] = c[i];
    }
  fread(c, 1, 256, pal_fp);
  for ( i = 0; i < 256; ++i )
    {
      b[i] = c[i];
    }
  i = fread(c, 1, 256, pal_fp);
  if ( i == 256 )
    {
      for ( i = 0; i < 256; ++i )
	{
	  a[i] = c[i];
	}
    }
  else
    {
      for ( i = 0; i < 256; ++i )
	{
	  a[i] = 0;
	}
    }
  fclose(pal_fp);
}

