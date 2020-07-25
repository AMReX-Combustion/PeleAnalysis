#include <AMReX_ParmParse.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_FillPatchUtil.H>
#include <StreamPC.H>

using namespace amrex;

static
Vector<Vector<Real>>
GetSeedLocations (const StreamParticleContainer& spc)
{
  Vector<Vector<Real>> locs;

  ParmParse pp;
  int nc=pp.countval("oneSeedPerCell");
  int ni=pp.countval("isoFile");
  int ns=pp.countval("seedLoc");
  int nrL=pp.countval("seedRakeL");
  int nrR=pp.countval("seedRakeR");
  AMREX_ALWAYS_ASSERT((nc>0) ^ ((ni>0) ^ ((ns>0) ^ ((nrL>0) && nrR>0))));
  if (nc>0)
  {
    int finestLevel = spc.numLevels() - 1;
    std::vector< std::pair<int,Box> > isects;
    FArrayBox mask;
    for (int lev=0; lev<=finestLevel; ++lev)
    {
      const auto& geom = spc.Geom(lev);
      const auto& dx = geom.CellSize();
      const auto& plo = geom.ProbLo();

      BoxArray baf;
      if (lev < finestLevel) {
        baf = BoxArray(spc.ParticleBoxArray(lev+1)).coarsen(spc.GetParGDB()->refRatio(lev));
      }
      for (MFIter mfi = spc.MakeMFIter(lev); mfi.isValid(); ++mfi)
      {
        const Box& tile_box  = mfi.tilebox();
        if (BL_SPACEDIM<3 || tile_box.contains(IntVect(D_DECL(0,50,107)))) {

          mask.resize(tile_box,1);
          mask.setVal(1);
          if (lev < finestLevel) {
            isects = baf.intersections(tile_box);
            for (const auto& p : isects) {
              mask.setVal(0,p.second,0,1);
            }
          }

          for (IntVect iv = tile_box.smallEnd(); iv <= tile_box.bigEnd(); tile_box.next(iv))
          {
            if (mask(iv,0) > 0)
            {
              locs.push_back({AMREX_D_DECL(plo[0] + (iv[0] + 0.5)*dx[0],
                                           plo[1] + (iv[1] + 0.5)*dx[1],
                                           plo[2] + (iv[2] + 0.5)*dx[2])});
            }
          }
        }
      }
    }
  }
  else if (ni>0)
  {
    // Read in isosurface
    AMREX_ALWAYS_ASSERT(AMREX_SPACEDIM==3);
    std::string isoFile; pp.get("isoFile",isoFile);
    if (ParallelDescriptor::IOProcessor())
      std::cerr << "Reading isoFile... " << isoFile << std::endl;

    std::ifstream ifs;
    ifs.open(isoFile.c_str());
    std::string line;
    std::getline(ifs,line);
    auto surfNames = Tokenize(line,std::string(", "));
    int nCompSeedNodes = surfNames.size();
    int nElts, nodesPerElt;
    ifs >> nElts;
    ifs >> nodesPerElt;

    FArrayBox tnodes;
    tnodes.readFrom(ifs);
    int nSeedNodes = tnodes.box().numPts();

    Real* ndat = tnodes.dataPtr();
    for (int i=0; i<nSeedNodes; ++i)
    {
      int o=i*nCompSeedNodes;
      locs.push_back({AMREX_D_DECL(ndat[o+0], ndat[o+1], ndat[o+2])});
    }
    tnodes.clear();

    Vector<int> faceData(nElts*nodesPerElt);
    ifs.read((char*)faceData.dataPtr(),sizeof(int)*faceData.size());
    ifs.close();
  }
  else if (pp.countval("seedLoc")>0)
  {
    Vector<Real> loc(BL_SPACEDIM);
    pp.getarr("seedLoc",loc,0,BL_SPACEDIM);
    locs.push_back({AMREX_D_DECL(loc[0], loc[1], loc[2])});
  }
  else
  {
    int seedRakeNum;
    pp.get("seedRakeNum",seedRakeNum);
    AMREX_ALWAYS_ASSERT(seedRakeNum >= 2);
    Vector<Real> locL(BL_SPACEDIM), locR(BL_SPACEDIM);
    pp.getarr("seedRakeL",locL,0,BL_SPACEDIM);
    pp.getarr("seedRakeR",locR,0,BL_SPACEDIM);

    for (int i=0; i<seedRakeNum; ++i) {
      locs.push_back({AMREX_D_DECL(locL[0] + (i/double(seedRakeNum-1))*(locR[0] - locL[0]),
                                   locL[1] + (i/double(seedRakeNum-1))*(locR[1] - locL[1]),
                                   locL[2] + (i/double(seedRakeNum-1))*(locR[2] - locL[2]))});
    }
  }
  return locs;
}

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
    int nGrow = 3;
    pp.query("nGrow",nGrow);
    AMREX_ALWAYS_ASSERT(nGrow>=1);
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
    }

    int Nsteps = 50;
    pp.query("Nsteps",Nsteps);
    StreamParticleContainer spc(Nsteps,geoms,dms,grids,ratios);

    auto locs = GetSeedLocations(spc);

    spc.InitParticles(locs);

    Real hRK = 0.1; pp.query("hRK",hRK);
    AMREX_ALWAYS_ASSERT(hRK>=0 && hRK<=0.5);
    Real dt = hRK * geoms[finestLevel].CellSize()[0];
    for (int step=0; step<Nsteps-1; ++step)
    {
      //Print() << "step " << step << std::endl;
      spc.ComputeNextLocation(step,dt,vectorField);
    }

    std::string outfile = "junkPlt";
    Print() << "Writing paticles to " << outfile << std::endl;
    spc.WritePlotFile(outfile, "particles");

    std::string tecfile = "tec.dat";
    Print() << "Writing streamlines in Tecplot ascii format to " << tecfile << std::endl;
    spc.WriteStreamAsTecplot(tecfile);
  }
  Finalize();
  return 0;
}
