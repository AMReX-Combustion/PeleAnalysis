#include <string>
#include <iostream>
#include <set>

#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include <AMReX_DataServices.H>
#include <AMReX_PlotFileUtil.H>

using namespace amrex;

static
void 
print_usage (int,
             char* argv[])
{
  std::cerr << "usage:\n";
  std::cerr << argv[0] << " infile infile=f1 [options] \n\tOptions:\n";
  exit(1);
}

std::string
getFileRoot(const std::string& infile)
{
  std::vector<std::string> tokens = Tokenize(infile,std::string("/"));
  return tokens[tokens.size()-1];
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

    if (pp.contains("verbose"))
      AmrData::SetVerbose(true);

    std::string plotFileName; pp.get("infile",plotFileName);
    DataServices::SetBatchMode();
    Amrvis::FileType fileType(Amrvis::NEWPLT);

    DataServices dataServices(plotFileName, fileType);
    if( ! dataServices.AmrDataOk()) {
      DataServices::Dispatch(DataServices::ExitRequest, NULL);
      // ^^^ this calls ParallelDescriptor::EndParallel() and exit()
    }
    AmrData& amrData = dataServices.AmrDataRef();

    int finestLevel = amrData.FinestLevel();
    pp.query("finestLevel",finestLevel);
    int Nlev = finestLevel + 1;

    Vector<int> is_per(AMREX_SPACEDIM,1);
    pp.queryarr("is_per",is_per,0,AMREX_SPACEDIM);
    Print() << "Periodicity assumed for this case: ";
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        Print() << is_per[idim] << " ";
    }

    RealBox rb(&(amrData.ProbLo()[0]),&(amrData.ProbHi()[0]));

    Vector<MultiFab> outdata(Nlev);
    Vector<Geometry> geoms(Nlev);

    int nGrow = 0;
    int nComp = amrData.NComp();

    Vector<std::string> inNames = amrData.PlotVarNames();
    Vector<std::string> outNames = amrData.PlotVarNames();
    
    Vector<int> destFillComps(nComp);
    for (int i=0; i<nComp; ++i) destFillComps[i] = i;

    for (int lev = 0; lev < Nlev; ++lev) {
      
      const BoxArray ba = amrData.boxArray(lev);
      const DistributionMapping dm(ba);

      outdata[lev] = MultiFab(ba,dm,nComp,nGrow);
      MultiFab indata(ba,dm,nComp,nGrow);

      int coord = 0;
      geoms[lev] = Geometry(amrData.ProbDomain()[lev],&rb, coord, &(is_per[0]));      
      
      Print() << "Reading data for level " << lev << std::endl;
      amrData.FillVar(indata,finestLevel,inNames,destFillComps);
      Print() << "Data has been read for level " << lev << std::endl;
      
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (MFIter mfi(indata,TilingIfNotGPU()); mfi.isValid(); ++mfi)
	{
	  const Box& bx = mfi.tilebox();
	  auto const& out_a = outdata[lev].array(mfi);
	  auto const& in_a = indata.array(mfi);
	  amrex::ParallelFor(bx, [=]
			     AMREX_GPU_DEVICE (int i, int j, int k) noexcept
				 {
				   // Do something funky
				   for (int n = 0; n < nComp; n++) {
				     out_a(i,j,k,n) = in_a(i,j,k,n); 
				   }
				 });
	}

    }
    
    std::string outfile(getFileRoot(plotFileName) + "_temp");
    Print() << "Writing new data to " << outfile << std::endl;
    Vector<int> isteps(Nlev, 0);
    Vector<IntVect> refRatios(Nlev-1,{AMREX_D_DECL(2, 2, 2)});
    amrex::WriteMultiLevelPlotfile(outfile, Nlev, GetVecOfConstPtrs(outdata), outNames,
                                   geoms, 0.0, isteps, refRatios);

  }
  Finalize();
  return 0;
}
