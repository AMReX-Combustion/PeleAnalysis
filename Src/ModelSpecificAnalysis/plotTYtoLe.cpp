#include <string>
#include <iostream>
#include <set>

#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include <AMReX_DataServices.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_BCRec.H>
#include <AMReX_Interpolater.H>

#include <AMReX_BLFort.H>
#include <mechanism.H>
#include <PelePhysics.H>

using namespace amrex;

pele::physics::PeleParams<pele::physics::transport::TransParm<
  pele::physics::PhysicsType::eos_type,
  pele::physics::PhysicsType::transport_type>>
  trans_parms;

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

    trans_parms.initialize();

    int finestLevel = amrData.FinestLevel();
    pp.query("finestLevel",finestLevel);
    int Nlev = finestLevel + 1;

    int idYin = -1;
    int idTin = -1;
    int idRin = -1;
    Vector<std::string> spec_names;
    pele::physics::eos::speciesNames<pele::physics::PhysicsType::eos_type>(spec_names);
    const Vector<std::string>& plotVarNames = amrData.PlotVarNames();
    const std::string spName= "Y(" + spec_names[0] + ")";
    const std::string TName = "temp";
    const std::string RName = "density";
    for (int i=0; i<plotVarNames.size(); ++i)
    {
      if (plotVarNames[i] == spName) idYin = i;
      if (plotVarNames[i] == TName)  idTin = i;
      if (plotVarNames[i] == RName)  idRin = i;
    }
    if (idYin<0 || idTin<0 || idRin<0)
      Print() << "Cannot find required data in pltfile" << std::endl;

    const int nCompIn  = NUM_SPECIES + 2;
    const int idLeout   = 0;
    const int nCompOut = idLeout + NUM_SPECIES;

    Vector<std::string> outNames(nCompOut);
    Vector<std::string> inNames(nCompIn);
    Vector<int> destFillComps(nCompIn);
    const int idYlocal = 0; // Xs start here
    const int idTlocal = NUM_SPECIES;   // T starts here
    const int idRlocal = NUM_SPECIES+1; // R starts here
    for (int i=0; i<NUM_SPECIES; ++i)
    {
      destFillComps[i] = idYlocal + i;
      inNames[i] =  "Y(" + spec_names[i] + ")";
      outNames[idLeout + i] = "Le(" + spec_names[i] + ")";
    }
    destFillComps[idTlocal] = idTlocal;
    destFillComps[idRlocal] = idRlocal;
    inNames[idTlocal] = TName;
    inNames[idRlocal] = RName;

    Vector<std::unique_ptr<MultiFab>> outdata(Nlev);
    Vector<std::unique_ptr<MultiFab>> tempdata(Nlev);
    Vector<Geometry> geoms(Nlev);
    amrex::RealBox real_box({AMREX_D_DECL(amrData.ProbLo()[0],
                                          amrData.ProbLo()[1],
                                          amrData.ProbLo()[2])},
                            {AMREX_D_DECL(amrData.ProbHi()[0],
                                          amrData.ProbHi()[1],
                                          amrData.ProbHi()[2])});
    amrex::Array<int, AMREX_SPACEDIM> is_periodic{AMREX_D_DECL(1, 1, 1)};
    geoms[0] = amrex::Geometry((amrData.ProbDomain())[0],
                               real_box,
                               amrData.CoordSys(),
                               is_periodic);
    const int nGrow = 0;

    for (int lev=0; lev<Nlev; ++lev)
    {
      const BoxArray ba = amrData.boxArray(lev);
      const DistributionMapping dm(ba);
      outdata[lev].reset(new MultiFab(ba,dm,nCompOut,nGrow));
      tempdata[lev].reset(new MultiFab(ba,dm,2*NUM_SPECIES+3,nGrow));
      if ( lev > 0 ) {
         geoms[lev] = amrex::refine(geoms[lev - 1], 2);
      }
      MultiFab indata(ba,dm,nCompIn,nGrow);

      Print() << "Reading data for level " << lev << std::endl;
      amrData.FillVar(indata,lev,inNames,destFillComps);
      Print() << "Data has been read for level " << lev << std::endl;

      // Get the transport data pointer
      auto const* ltransparm = trans_parms.device_parm();

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (amrex::MFIter mfi(indata, amrex::TilingIfNotGPU()); mfi.isValid();
           ++mfi) {

        const Box& bx = mfi.tilebox();
        Array4<Real> const& Y_a = indata.array(mfi,idYlocal);
        Array4<Real> const& T_a = indata.array(mfi,idTlocal);
        Array4<Real> const& rho_a = indata.array(mfi,idRlocal);
        Array4<Real> const& D_a = tempdata[lev]->array(mfi,0);
	Array4<Real> const& chi_a = tempdata[lev]->array(mfi,NUM_SPECIES);
        Array4<Real> const& mu_a = tempdata[lev]->array(mfi,2*NUM_SPECIES);
        Array4<Real> const& xi_a = tempdata[lev]->array(mfi,2*NUM_SPECIES+1);
        Array4<Real> const& lam_a = tempdata[lev]->array(mfi,2*NUM_SPECIES+2);
        Array4<Real> const& Le_a = outdata[lev]->array(mfi,idLeout);

        amrex::launch(bx, [=] AMREX_GPU_DEVICE(amrex::Box const& tbx) {
          auto trans = pele::physics::PhysicsType::transport();
          trans.get_transport_coeffs(
				     tbx, Y_a, T_a, rho_a, D_a, chi_a, mu_a, xi_a, lam_a, ltransparm);
        });
        amrex::ParallelFor(bx, [=]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
          Real Cpmix = 0.0;
          Real Yloc[NUM_SPECIES] = {0.0};
          for (int n=0; n<NUM_SPECIES; ++n) {
              Yloc[n] = Y_a(i,j,k,n);
          }
          CKCPBS(T_a(i,j,k),Yloc,Cpmix);
          for (int n=0; n<NUM_SPECIES; ++n) {
            Le_a(i,j,k,n) = D_a(i,j,k,n) / (lam_a(i,j,k) / Cpmix);
          }
        });
      }

      Print() << "Derive finished for level " << lev << std::endl;
    }

    std::string outfile(getFileRoot(plotFileName) + "_Le");
    Print() << "Writing new data to " << outfile << std::endl;
    Vector<int> isteps(Nlev, 0);
    Vector<IntVect> refRatios(Nlev-1,{AMREX_D_DECL(2, 2, 2)});
    amrex::WriteMultiLevelPlotfile(outfile, Nlev, GetVecOfConstPtrs(outdata), outNames,
                                   geoms, 0.0, isteps, refRatios);
  }
  Finalize();
  return 0;
}
