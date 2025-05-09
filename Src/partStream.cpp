#include <AMReX_ParmParse.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_FillPatchUtil.H>
#include <StreamPC.H>

using namespace amrex;
static
Vector<Vector<Real>>
GetSeedLocations (const StreamParticleContainer& spc, Vector<int>& faceData)
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
        if (AMREX_SPACEDIM<3 || tile_box.contains(IntVect(AMREX_D_DECL(0,50,107)))) {

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
    std::string isoFile; pp.get("isoFile",isoFile);
    if (ParallelDescriptor::IOProcessor())
      std::cerr << "Reading isoFile... " << isoFile << std::endl;

    std::ifstream ifs;
    ifs.open(isoFile.c_str());
    // AJA added dummy line read; sometimes time, sometimes `decimated'
    std::string topline;
    std::getline(ifs,topline);
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

    faceData.resize(nElts*nodesPerElt);
    ifs.read((char*)faceData.dataPtr(),sizeof(int)*faceData.size());
    ifs.close();
  }
  else if (pp.countval("seedLoc")>0)
  {
    Vector<Real> loc(AMREX_SPACEDIM);
    pp.getarr("seedLoc",loc,0,AMREX_SPACEDIM);
    locs.push_back({AMREX_D_DECL(loc[0], loc[1], loc[2])});
  }
  else
  {
    int seedRakeNum;
    pp.get("seedRakeNum",seedRakeNum);
    AMREX_ALWAYS_ASSERT(seedRakeNum >= 2);
    Vector<Real> locL(AMREX_SPACEDIM), locR(AMREX_SPACEDIM);
    pp.getarr("seedRakeL",locL,0,AMREX_SPACEDIM);
    pp.getarr("seedRakeR",locR,0,AMREX_SPACEDIM);

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
    Vector<std::string> vectorVarNames(AMREX_SPACEDIM);
    pp.getarr("vectorField",vectorVarNames);
    int nVars= pp.countval("vars");
    Vector<std::string> inVarNames(nVars);
    pp.getarr("vars",inVarNames);
    PlotFileData pf(infile);
    int finestLevel = pf.finestLevel();
    Vector<Geometry> geoms(finestLevel+1);
    Vector<BoxArray> grids(finestLevel+1);
    Vector<DistributionMapping> dms(finestLevel+1);
    Vector<int> ratios(finestLevel);
    int nComp = AMREX_SPACEDIM + nVars;

    
    
    // get base for output files
    std::string outfile = infile;
    pp.query("outfile",outfile);

    int writeParticles(0);
    pp.query("writeParticles",writeParticles);
    std::string particlefile = outfile+"_particles";
    pp.query("particlefile",particlefile);
    int writeStreams(0);
    pp.query("writeStreams",writeStreams);
    std::string streamfile = outfile+"_stream";
    pp.query("streamfile",streamfile);
    //
    int writeStreamBin(1);
    pp.query("writeStreamBin",writeStreamBin);
    std::string streamBinfile = outfile+"_streamBin";
    pp.query("streamBinfile",streamBinfile);
    
    Vector<std::string> outVarNames = {AMREX_D_DECL("X","Y","Z")};
    for (int n = 0; n < nVars; n++) {
      outVarNames.push_back(inVarNames[n]);
    }

    IntVect pp_is_per;
    pp.getarr("is_per",pp_is_per);
    Array<int,AMREX_SPACEDIM> is_per = {AMREX_D_DECL(pp_is_per[0],pp_is_per[1],pp_is_per[2])};
    RealBox rb(pf.probLo(),pf.probHi());

    int Nlev = finestLevel + 1;
    Vector<Vector<MultiFab>> pfdata(Nlev);
    for (int lev=0; lev<Nlev; ++lev) {
      geoms[lev].define(pf.probDomain(lev),rb,pf.coordSys(),is_per);
      grids[lev] = pf.boxArray(lev);
      dms[lev] = pf.DistributionMap(lev);
      if (lev < finestLevel) ratios[lev] = pf.refRatio(lev);

      pfdata[lev].resize(nComp);
      Print() << "Loading data on level " << lev << std::endl;
      for (int d=0; d<AMREX_SPACEDIM; ++d) {
        pfdata[lev][d] = pf.get(lev,vectorVarNames[d]);
      }
      for (int n = 0; n < nVars; n++) {
	pfdata[lev][n+AMREX_SPACEDIM] = pf.get(lev,inVarNames[n]);
      }
    }

    int nGrow = 3;
    pp.query("nGrow",nGrow);
      
    Real time=0;
    PhysBCFunctNoOp f;
    PCInterp cbi;
    BCRec bc;
    AMREX_ALWAYS_ASSERT(nGrow>=1);
    Vector<MultiFab> vectorField(Nlev);
    for (int lev=0; lev<Nlev; ++lev) {
      vectorField[lev].define(grids[lev],dms[lev],nComp,nGrow);
      for (int d=0; d<nComp; ++d) {
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
    StreamParticleContainer spc(Nsteps,geoms,dms,grids,ratios,nComp,pp_is_per,outVarNames);
    
    Vector<int> faceData;
    auto locs = GetSeedLocations(spc,faceData);
    if (writeStreamBin == 1 && faceData.empty()) {
      Abort("Writing stream binary without surface definition!");
    }
    int nStreamPairs = locs.size();
    // Initialise particles
    Print() << "Initialising particles..." << std::endl;
    spc.InitParticles(locs);
    
    //Check if particles initialised fine
    if (spc.OK()) {
      Print() << "SPC is happy with initialisation :-)" << std::endl;
    }
    
#if AMREX_DEBUG //AJA 
    Print() << "Checking if initialised properly..." << std::endl;
    spc.InspectParticles(nStreamPairs);
#endif
   
    // Interpolate at start
    Print() << "Interpolation at the seed points..." << std::endl;
    spc.InterpDataAtLocation(0,vectorField);

    
    Print() << "Computing streams and interpolating..." << std::endl;
    int cSpace=0; pp.query("cSpace",cSpace);
    Real hRK = 0.1; pp.query("hRK",hRK);
    AMREX_ALWAYS_ASSERT((cSpace == 1) || (hRK>=0 && hRK<=0.5));
    //Real dt = hRK;// * geoms[finestLevel].CellSize()[0];
    for (int step=0; step<Nsteps-1; ++step)
    {
      // find next location
      spc.ComputeNextLocation(step,hRK,vectorField,cSpace);

      // interpolate all data
      spc.InterpDataAtLocation(step+1,vectorField);

    }

    // check in again
    if (spc.OK()) {
      Print() << "SPC is happy afterwards :-)" << std::endl;
    }
    
#if AMREX_DEBUG //AJA checks
    Print() << "Checking if we broke things afterwards..." << std::endl;
    spc.InspectParticles(nStreamPairs);
#endif
    //check we didn't break them again
    spc.OK();
    
    if (writeParticles) {
      Print() << "Writing particles in plotfile to " << particlefile << std::endl;
      spc.WritePlotFile(particlefile, "particles");
    }
    if (writeStreams) {
      Print() << "Writing streamlines in Tecplot ascii format to " << streamfile << std::endl;
      spc.WriteStreamAsTecplot(streamfile);
    }
        //
    // Write streamBin
    //
    if (writeStreamBin) {
      Print() << "Writing streamlines as binary " << streamBinfile << std::endl;
      spc.WriteStreamAsBinary(streamBinfile,faceData);
    }
    
  }
  Finalize();
  return 0;
}
