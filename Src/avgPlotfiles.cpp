#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MultiFabUtil.H>
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
      BL_ASSERT(sComp+nComp <= amrData0.NComp());
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

    Vector<Real> plo(BL_SPACEDIM), phi(BL_SPACEDIM);
    Vector<Box> psize(finestLevel+1);
    const IntVect ilo = subbox.smallEnd();
    const IntVect ihi = subbox.bigEnd();

    for (int i =0 ; i< BL_SPACEDIM; i++) {   
       plo[i] = amrData0.ProbLo()[i]+(ilo[i])*amrData0.DxLevel()[finestLevel][i];
       phi[i] = amrData0.ProbLo()[i]+(ihi[i]+1)*amrData0.DxLevel()[finestLevel][i];
    }

    const Vector<int>& ratio = amrData0.RefRatio();

    Vector<Box> domain(finestLevel+1);
    Vector<int> ratioTot(finestLevel,1);
    domain[finestLevel] = subbox;
    for (int lev=finestLevel-1; lev>=0; --lev) {
      domain[lev] = coarsen(domain[lev+1],ratio[lev]);
    }
    for (int lev=1; lev<=finestLevel; ++lev) {
      domain[lev] = refine(domain[lev-1],ratio[lev-1]);
      ratioTot[lev-1] = (lev == 1 ? ratio[lev-1] : ratioTot[lev-2] * ratio[lev-1]);
    }
    AMREX_ALWAYS_ASSERT(domain[0].numPts() > 0);
    Vector<Geometry> geom(finestLevel+1);
    RealBox rb(AMREX_D_DECL(plo[0],plo[1],plo[2]), AMREX_D_DECL(phi[0],phi[1],phi[2]));
    Array<int,AMREX_SPACEDIM> is_per = {0};
    for (int lev=0; lev<=finestLevel; ++lev) {
      geom[lev].define(domain[lev],rb,amrData0.CoordSys(),is_per);
    }
    
    BoxArray ba_res(domain[0]);
    int max_grid_size = 32; pp.query("max_grid_size",max_grid_size);
    ba_res.maxSize(max_grid_size);

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
            average_down(amrData.GetGrids(lev+1,comps[i],domain[lev+1]),alias,geom[lev+1],geom[lev],0,1,ratioTot[lev]);
          }
          amrData.FlushGrids(comps[i]);
        }
      }

      MultiFab::Add(mf_avg,mf_avgDn,0,0,comps.size(),0);
    }

    mf_avg.mult(1.0/nf);
    writePlotFile(outfile.c_str(),mf_avg,geom[0],IntVect(AMREX_D_DECL(2,2,2)),0.0,pltnames);
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

  IntVect img(D_DECL(box.length(d[0]) - 1,box.length(d[1]) - 1,0));
  image.resize(Box(IntVect::TheZeroVector(),img),1);

  IntVect div;
  if (slicedir>=0 && slicedir<BL_SPACEDIM)
      div[slicedir] = sliceloc;

  for (int i=se[d[0]]; i<=be[d[0]]; ++i) {
    for (int j=se[d[1]]; j<=be[d[1]]; ++j) {
      div[d[0]] = i;
      div[d[1]] = j;
      image(IntVect(D_DECL(i - se[d[0]],j - se[d[1]],0)),0) =
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

