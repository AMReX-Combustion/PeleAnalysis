
#include <AMReX.H>
#include <AMReX_Print.H>
#include <AMReX_ParmParse.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_DataServices.H>

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

int main(int argc, char* argv[])
{
  amrex::Initialize(argc,argv);
  {
    ParmParse pp;
    std::string file;
    pp.get("file",file);
    DataServices::SetBatchMode();
    Amrvis::FileType fileType(Amrvis::NEWPLT);
    DataServices dataServices(file, fileType);
    if( ! dataServices.AmrDataOk()) {
      DataServices::Dispatch(DataServices::ExitRequest, NULL);
    }
    AmrData& amrData = dataServices.AmrDataRef();

    int finestLevel = amrData.FinestLevel();
    pp.query("finestLevel",finestLevel);

    const std::vector<std::string>& pieces = Tokenize(file,std::string("/"));

    int slicedir, sliceloc;
    pp.get("slicedir",slicedir);
    pp.get("sliceloc",sliceloc);
    std::string varname;
    pp.get("varname",varname);
    AMREX_ALWAYS_ASSERT(amrData.CanDerive(varname));

    Box domain = amrData.ProbDomain()[finestLevel];
    domain.setSmall(slicedir,sliceloc);
    domain.setBig(slicedir,sliceloc);
    FArrayBox data(domain,1);
    amrData.FillVar(&data,data.box(),finestLevel,varname,0);
    Print() << "min,max: " << data.min() << ", " << data.max() << std::endl;

    std::string outtype("image"); pp.query("outtype",outtype);
    if (outtype=="image" || outtype=="gray")
    {
      Real data_min = data.min();
      Real data_max = data.max();
      pp.query("min",data_min);
      pp.query("max",data_max);
      
      BaseFab<int> image;
      const int nVals = 256;
      pixelizeData(data,slicedir,sliceloc,image,data_min,data_max,nVals);
      
      const int width = image.box().length(0);
      const int height= image.box().length(1);

      if (outtype=="image")
      {
        std::string palfile;
        pp.get("palette",palfile);
        int r[256], g[256], b[256], t[256];
        LOAD_PALETTE_STR(palfile,r,g,b,t);
            
        std::string outfile = pieces[pieces.size()-1] + PPM;
        pp.query("outfile",outfile);
        STORE_PPM_STR(outfile, width, height, image.dataPtr(), r, g, b);
      }
      else
      {
        std::string outfile = pieces[pieces.size()-1] + PGM;
        pp.query("outfile",outfile);
        STORE_PGM_STR(outfile, width, height, image.dataPtr());
      }
    }
  }
  amrex::Finalize();
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

