#include <AMReX_ParmParse.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_FillPatchUtil.H>
#include <AMReX_Extrapolater.H>
#include <StreamPC.H>

using namespace amrex;

int
main (int   argc,
      char* argv[])
{
  Initialize(argc,argv);
  {
    ParmParse pp;

    std::string infile; pp.get("infile",infile);
    Vector<std::string> inVarNames = {D_DECL("x_velocity", "y_velocity", "z_velocity")};

    PlotFileData pf(infile);
    int finestLevel = pf.finestLevel();
    Vector<Geometry> geoms(finestLevel+1);
    Vector<BoxArray> grids(finestLevel+1);
    Vector<DistributionMapping> dms(finestLevel+1);
    Vector<int> ratios(finestLevel);

    Array<int,AMREX_SPACEDIM> is_per = {D_DECL(0, 0, 0)};
    RealBox rb(pf.probLo(),pf.probHi());

    int Nlev = finestLevel + 1;
    Vector<Vector<MultiFab>> pfdata(Nlev);
    for (int lev=0; lev<Nlev; ++lev) {
      geoms[lev].define(pf.probDomain(lev),rb,pf.coordSys(),is_per);
      grids[lev] = pf.boxArray(lev);
      dms[lev] = pf.DistributionMap(lev);
      if (lev < finestLevel) ratios[lev] = pf.refRatio(lev);

      pfdata[lev].resize(AMREX_SPACEDIM);
      for (int d=0; d<AMREX_SPACEDIM; ++d) {
        pfdata[lev][d] = pf.get(lev,inVarNames[d]);
      }
    }

    Real time=0;
    PhysBCFunctNoOp f;
    PCInterp cbi;
    BCRec bc;
    int nGrow = 1;
    int nComp = inVarNames.size();
    Vector<MultiFab> vectorField(Nlev);
    for (int lev=0; lev<Nlev; ++lev) {
      vectorField[lev].define(grids[lev],dms[lev],nComp,nGrow);
      for (int d=0; d<AMREX_SPACEDIM; ++d) {
        if (lev==0) {
          FillPatchSingleLevel(vectorField[lev],time,{&pfdata[lev][d]},{time},0,d,1,geoms[0],f,0);
        }
        else
        {
          FillPatchTwoLevels(vectorField[lev],time,{&pfdata[lev-1][d]},{time},{&pfdata[lev][d]},{time},0,d,1,
                             geoms[lev-1],geoms[lev],f,0,f,0,ratios[lev-1]*IntVect::Unit,&cbi,{bc},0);
        }
      }
      vectorField[lev].FillBoundary(geoms[lev].periodicity());
      Extrapolater::FirstOrderExtrap(vectorField[lev],geoms[lev],0,AMREX_SPACEDIM);
    }

    StreamParticleContainer spc(geoms,dms,grids,ratios);
    spc.InitParticles();

    Real hRK = 0.3; pp.query("hRK",hRK);
    int redist_int = (int)(nGrow / hRK);
    int Nsteps = RealData::nPointOnStream-1;

    Real dt = 1.e-6;
    for (int step=0; step<Nsteps; ++step)
    {
      Print() << "Step " << step << std::endl;
      
      if (step % redist_int == 0) spc.Redistribute();

      spc.ComputeNextLocation(step,dt,vectorField);
    }

    spc.WritePlotFile("junkPlt", "particles");
    //spc.WriteAsciiFile ("part");
  }
  Finalize();
  return 0;
}
