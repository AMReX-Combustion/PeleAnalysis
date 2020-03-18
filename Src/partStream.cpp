#include <AMReX_ParmParse.H>
#include <AMReX_PlotFileUtil.H>
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
    Vector<Vector<MultiFab>> pfdata(finestLevel+1);
    Vector<Geometry> geoms(finestLevel+1);
    Vector<BoxArray> grids(finestLevel+1);
    Vector<DistributionMapping> dms(finestLevel+1);
    Vector<int> ratios(finestLevel);

    Array<int,AMREX_SPACEDIM> is_per = {D_DECL(0, 0, 0)};
    RealBox rb(pf.probLo(),pf.probHi());
    
    int nComp = inVarNames.size();
    for (int lev=0; lev<=finestLevel; ++lev) {
      pfdata[lev].resize(nComp);
      for (int n=0; n<nComp; ++n) {
        pfdata[lev][n] = pf.get(lev,inVarNames[n]);
      }

      geoms[lev].define(pf.probDomain(lev),rb,pf.coordSys(),is_per);
      grids[lev] = pf.boxArray(lev);
      dms[lev] = pf.DistributionMap(lev);
      if (lev < finestLevel) ratios[lev] = pf.refRatio(lev);
    }

    StreamParticleContainer spc(geoms,dms,grids,ratios);
    spc.InitParticles();

    spc.WritePlotFile("junkPlt", "particles");

  }
  Finalize();
  return 0;
}
