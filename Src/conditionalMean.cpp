#include <string>
#include <iostream>
using std::cout;
using std::cerr;
using std::endl;

#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include <AMReX_DataServices.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_MLMG.H>
#include <AMReX_MLPoisson.H>
#include <AMReX_MLABecLaplacian.H>
#include <mixfracBilger.H>
#include <mechanism.H>
#include <PelePhysics.H>

using namespace amrex;

static
void 
print_usage (int,
             char* argv[])
{
    std::cerr << "usage:\n";
    std::cerr << argv[0] << " infile=<filename>\n"
         << "   binComp=i            : Variable id number to condition on\n"
         << "   AvgComps=j k l       : Variable id numbers to average\n"
         << "   min=m;  max=m        : min/max values for bins%i\n"
              << "   nBins=n              : Number of bins in PDF (default=64)\n"
         << "   finestLevel=n        : Finest level at which to evaluate PDF\n"
         << "   outSuffix=str        : Suffix to add to the pltfile name as an alt dir for results (default="")\n"
         << "   aja=true/false       : Put the header in a separate file for gnuplot/matlab (default=false)\n"
              << "   infile= plt1 plt2    : List of plot files" << std::endl;
    exit(1);
}

typedef BaseFab<int> IntFab;
typedef FabArray<IntFab> IntMF;

void
parseComposition(Vector<std::string> compositionIn,
                 std::string         compositionType,
                 Real               *massFrac)
{
   Real compoIn[NUM_SPECIES] = {0.0};

   // Get species names
   Vector<std::string> specNames;
   pele::physics::eos::speciesNames<pele::physics::PhysicsType::eos_type>(specNames);

   // For each entry in the user-provided composition, parse name and value
   std::string delimiter = ":";
   int specCountIn = compositionIn.size();
   for (int i = 0; i < specCountIn; i++ ) {
      long unsigned sep = compositionIn[i].find(delimiter);
      if ( sep == std::string::npos ) {
         Abort("Error parsing '"+compositionIn[i]+"' --> unable to find delimiter :");
      }
      std::string specNameIn = compositionIn[i].substr(0, sep);
      Real value = std::stod(compositionIn[i].substr(sep+1,compositionIn[i].length()));
      int foundIt = 0;
      for (int k = 0; k < NUM_SPECIES; k++ ) {
         if ( specNameIn == specNames[k] ) {
            compoIn[k] = value;
            foundIt = 1;
         }
      }
      if ( !foundIt ) {
         Abort("Error parsing '"+compositionIn[i]+"' --> unable to match to any species name");
      }
   }

   // Ensure that it sums to 1.0:
   Real sum = 0.0;
   for (int k = 0; k < NUM_SPECIES; k++ ) {
      sum += compoIn[k];
   }
   for (int k = 0; k < NUM_SPECIES; k++ ) {
      compoIn[k] /= sum;
   }

   // Fill the massFrac array, convert from mole fraction if necessary
   if ( compositionType == "mass" ) {                // mass
      for (int i = 0; i < NUM_SPECIES; i++ ) {
         massFrac[i] = compoIn[i];
      }
   } else if ( compositionType == "mole" ) {         // mole
      auto eos = pele::physics::PhysicsType::eos();
      eos.X2Y(compoIn,massFrac);
   } else {
      Abort("Unknown mixtureFraction.type ! Should be 'mass' or 'mole'");
   }
}

void
initMixtureFraction(mixFracBilgerData &a_mixFrac)
{
    // Get default fuel and oxy tank composition: pure fuel vs air
    Vector<std::string> specNames;
    pele::physics::eos::speciesNames<pele::physics::PhysicsType::eos_type>(specNames);
    amrex::Real YF[NUM_SPECIES], YO[NUM_SPECIES];
    for (int i=0; i<NUM_SPECIES; ++i) {
        YF[i] = 0.0;
        YO[i] = 0.0;
        if (!specNames[i].compare("O2"))  YO[i] = 0.233;
        if (!specNames[i].compare("N2"))  YO[i] = 0.767;
    }

    auto eos = pele::physics::PhysicsType::eos();
    // Overwrite with user-defined value if provided in input file
    ParmParse pp("mixtureFraction");
    std::string MFformat;
    int hasUserMF = pp.contains("format");
    if ( hasUserMF ) {
        pp.query("format", MFformat);
        if ( !MFformat.compare("Cantera")) {             // use a Cantera-like format with <SpeciesName>:<Value>, default in 0.0
            std::string MFCompoType;
            pp.query("type", MFCompoType);
            Vector<std::string> compositionIn;
            int entryCount = pp.countval("oxidTank");
            compositionIn.resize(entryCount);
            pp.getarr("oxidTank",compositionIn,0,entryCount);
            parseComposition(compositionIn, MFCompoType, YO);
            entryCount = pp.countval("fuelTank");
            compositionIn.resize(entryCount);
            pp.getarr("fuelTank",compositionIn,0,entryCount);
            parseComposition(compositionIn, MFCompoType, YF);
        } else if ( !MFformat.compare("RealList")) {     // use a list of Real. MUST contains an entry for each species in the mixture
            std::string MFCompoType;
            pp.query("type", MFCompoType);
            if ( !MFCompoType.compare("mass") ) {
               int entryCount = pp.countval("oxidTank");
               AMREX_ALWAYS_ASSERT(entryCount==NUM_SPECIES);
               Vector<amrex::Real> compositionIn(NUM_SPECIES);
               pp.getarr("oxidTank",compositionIn,0,NUM_SPECIES);
               for (int i=0; i<NUM_SPECIES; ++i) {
                  YO[i] = compositionIn[i];
               }
               entryCount = pp.countval("fuelTank");
               AMREX_ALWAYS_ASSERT(entryCount==NUM_SPECIES);
               pp.getarr("fuelTank",compositionIn,0,NUM_SPECIES);
               for (int i=0; i<NUM_SPECIES; ++i) {
                  YF[i] = compositionIn[i];
               }
            } else if ( !MFCompoType.compare("mole") ) {
               amrex::Real XF[NUM_SPECIES], XO[NUM_SPECIES];
               int entryCount = pp.countval("oxidTank");
               AMREX_ALWAYS_ASSERT(entryCount==NUM_SPECIES);
               Vector<amrex::Real> compositionIn(NUM_SPECIES);
               pp.getarr("oxidTank",compositionIn,0,NUM_SPECIES);
               for (int i=0; i<NUM_SPECIES; ++i) {
                  XO[i] = compositionIn[i];
               }
               entryCount = pp.countval("fuelTank");
               AMREX_ALWAYS_ASSERT(entryCount==NUM_SPECIES);
               pp.getarr("fuelTank",compositionIn,0,NUM_SPECIES);
               for (int i=0; i<NUM_SPECIES; ++i) {
                  XF[i] = compositionIn[i];
               }

               eos.X2Y(XO,YO);
               eos.X2Y(XF,YF);
            } else {
               Abort("Unknown mixtureFraction.type ! Should be 'mass' or 'mole'");
            }
        } else {
            Abort("Unknown mixtureFraction.format ! Should be 'Cantera' or 'RealList'");
        }
    }

    // Only interested in CHON -in that order. Compute Bilger weights
    amrex::Real atwCHON[4] = {0.0};
    pele::physics::eos::atomic_weightsCHON<pele::physics::PhysicsType::eos_type>(atwCHON);
    a_mixFrac.Beta_mix[0] = ( atwCHON[0] != 0.0 ) ? 2.0/atwCHON[0] : 0.0;
    a_mixFrac.Beta_mix[1] = ( atwCHON[1] != 0.0 ) ? 1.0/(2.0*atwCHON[1]) : 0.0;
    a_mixFrac.Beta_mix[2] = ( atwCHON[2] != 0.0 ) ? -1.0/atwCHON[2] : 0.0;
    a_mixFrac.Beta_mix[3] = 0.0;

    // Compute each species weight for the Bilger formulation based on elemental compo
    // Only interested in CHON -in that order.
    int ecompCHON[NUM_SPECIES*4];
    pele::physics::eos::element_compositionCHON<pele::physics::PhysicsType::eos_type>(ecompCHON);
    amrex::Real mwt[NUM_SPECIES];
    eos.molecular_weight(mwt);
    a_mixFrac.Zfu = 0.0;
    a_mixFrac.Zox = 0.0;
    for (int i=0; i<NUM_SPECIES; ++i) {
        a_mixFrac.spec_Bilger_fact[i] = 0.0;
        for (int k = 0; k < 4; k++) {
            a_mixFrac.spec_Bilger_fact[i] += a_mixFrac.Beta_mix[k] * (ecompCHON[i*4 + k]*atwCHON[k]/mwt[i]);
        }
        a_mixFrac.Zfu += a_mixFrac.spec_Bilger_fact[i]*YF[i];
        a_mixFrac.Zox += a_mixFrac.spec_Bilger_fact[i]*YO[i];
    }
}


int
main (int   argc,
      char* argv[])
{
    amrex::Initialize(argc,argv);

    if (argc < 2)
        print_usage(argc,argv);

    ParmParse pp;

    if (pp.contains("help"))
        print_usage(argc,argv);

    bool isioproc = ParallelDescriptor::IOProcessor();
    int ioproc = ParallelDescriptor::IOProcessorNumber();
    int verbose=0; pp.query("verbose",verbose);
    if ( !(isioproc) ) {
        verbose=0;
    }

    if (verbose>1)
        AmrData::SetVerbose(true);

    // Get plot file
    int nPlotFiles(pp.countval("infile"));
    if(nPlotFiles <= 0) {
        std::cerr << "Bad nPlotFiles:  " << nPlotFiles << std::endl;
        std::cerr << "Exiting." << std::endl;
        DataServices::Dispatch(DataServices::ExitRequest, NULL);
    }
    if (verbose)
        std::cout << "Processing " << nPlotFiles << " plotfiles..." << std::endl;

    std::string outSuffix = ""; pp.query("outSuffix",outSuffix);

    mixFracBilgerData mFrac_d;
    initMixtureFraction(mFrac_d);

    // Transport pointer
    static pele::physics::transport::TransportParams<
              pele::physics::PhysicsType::transport_type> trans_parms;
    trans_parms.allocate();

    // Make an array of srings containing paths of input plot files
    Vector<std::string> plotFileNames(nPlotFiles);
    for(int iPlot = 0; iPlot < nPlotFiles; ++iPlot) {
        pp.get("infile", plotFileNames[iPlot], iPlot);
        if (verbose)
            std::cout << "   " << plotFileNames[iPlot] << std::endl;
    }

    // Get finest level argument
    int finestLevel(-1);
    pp.query("finestLevel",finestLevel);

    // Number of bins
    int nBins(64); pp.query("nBins",nBins);

    // Binning on Z-coor
    bool doZcoorBins = false; pp.query("doZcoorBin",doZcoorBins);
    int nBinsZ = 1; pp.query("nBinsZ",nBinsZ);

    // Variables to bin, and to average
    int binComp = -1; pp.get("binComp",binComp);
    bool useMixFrac = false; pp.query("useMixFrac",useMixFrac);
    bool writePltMixFrac = false; pp.query("plotMixFrac",writePltMixFrac);
    bool isPeleC = false; pp.query("isPeleC",isPeleC);

    // Get species names
    Vector<std::string> spec_names;
    pele::physics::eos::speciesNames<pele::physics::PhysicsType::eos_type>(spec_names);
    // Species + temp + volFrac
    int readInnComp = NUM_SPECIES + 8;
    Vector<std::string> readInComps(readInnComp,"");
    readInComps[0] = "volFrac";
    for (int comp{1}; comp < NUM_SPECIES+1; ++comp) {
       readInComps[comp] = "Y("+spec_names[comp-1]+")";
    }
    readInComps[NUM_SPECIES+1] = "temp";
    readInComps[NUM_SPECIES+2] = "density";
    readInComps[NUM_SPECIES+3] = "HeatRelease";
    readInComps[NUM_SPECIES+4] = "I_R(OH)";
    readInComps[NUM_SPECIES+5] = "I_R(H2O2)";
    readInComps[NUM_SPECIES+6] = "I_R(CH2O)";
    readInComps[NUM_SPECIES+7] = "I_R(OC12H23OOH)";

    // List of field to average and map into mfVec multifab components
    // Y(OH), Y(H2O2), Y(CH2O), Y(OC12H23OOH), temp, HR, I_R(OH), I_R(H2O2), I_R(CH2O), I_R(OC12H23OOH), ScalDiss
    const int nAvgComps = 11;
    Array<int,nAvgComps> compIdx = {4,8,12,34,NUM_SPECIES+1,NUM_SPECIES+3,NUM_SPECIES+4,NUM_SPECIES+5,
                                    NUM_SPECIES+6,NUM_SPECIES+7,NUM_SPECIES+8};
    
    Real binMin=0; pp.get("binMin",binMin);
    Real binMax=1; pp.get("binMax",binMax);
    if (binMax <= binMin) {
        amrex::Abort("Bad bin min,max");
    }

    bool floor=false; pp.query("floor",floor);
    bool ceiling=false; pp.query("ceiling",ceiling);

    Vector<int> destFillComps(readInnComp);
    for (int i=0; i<destFillComps.size(); ++i) {
        destFillComps[i] = i;
    }
    
    // Adding a scaldiss and mask
    int nComp = readInnComp + 2;
    Vector<Real> dataIV(nComp);

    // List of outgoing data fields
    int nOutComp = nAvgComps;
    Vector<std::string> outComps{"Y(OH)","Y(H2O2)","Y(CH2O)","Y(OC12H23OOH)",
                                 "temp","HeatRelease","I_R(OH)","I_R(H2O2)",
                                 "I_R(CH2O)","I_R(OC12H23OOH)","ScalDiss"};

    bool aja(false);
    pp.query("aja",aja);
    if (aja && isioproc) {
      std::cout << "Output for aja" << std::endl;
    }

    // Loop over files
    for (int iPlot=0; iPlot<nPlotFiles; iPlot++) {

        const Real plt_time_io = ParallelDescriptor::second();
        
        // Open file and get an amrData pointer
        const std::string& infile = plotFileNames[iPlot];
        if (verbose) {
            std::cout << "\nOpening " << infile << "..." << std::endl;
        }
        DataServices::SetBatchMode();
        Amrvis::FileType fileType(Amrvis::NEWPLT);
        DataServices dataServices(infile, fileType);
        if (!dataServices.AmrDataOk()) {
            DataServices::Dispatch(DataServices::ExitRequest, NULL);
        }

        AmrData& amrData = dataServices.AmrDataRef();
        int ngrow = 0;

        Box domain;
        if (iPlot==0)
        {
            domain = amrData.ProbDomain()[0];

            if (finestLevel<0) {
                finestLevel = amrData.FinestLevel();
            }
        }

        int thisFinestLevel = std::min(finestLevel,amrData.FinestLevel());
        int Nlev = thisFinestLevel+1;

        // Geometry data
        RealBox rb(&(amrData.ProbLo()[0]), 
                   &(amrData.ProbHi()[0]));
        int coord = 0;
        Vector<int> is_per(AMREX_SPACEDIM,0);

        // Solution containers
        Vector<MultiFab> mixfracVec(Nlev);
        Vector<MultiFab> mfVec(Nlev);
        Vector<DistributionMapping> dmaps(Nlev);
        Vector<BoxArray> grids(Nlev);
        Vector<Geometry> geoms(Nlev);
        for (int iLevel=0; iLevel<Nlev; ++iLevel) {
            grids[iLevel] = amrData.boxArray(iLevel);
            dmaps[iLevel] = DistributionMapping(grids[iLevel]);
            mfVec[iLevel].define(grids[iLevel],dmaps[iLevel],nComp, ngrow);
            mixfracVec[iLevel].define(grids[iLevel],dmaps[iLevel],AMREX_SPACEDIM+2, 1); // Contains MixFrac, Diffus, MixFrac_gradients
            geoms[iLevel] = Geometry(amrData.ProbDomain()[iLevel],&rb,coord,&(is_per[0]));
        }
   
        Real binZMin=amrData.ProbLo()[2]; pp.query("binZMin",binZMin);
        Real binZMax=amrData.ProbHi()[2]; pp.query("binZMax",binZMax);

        // Conditionals containers
        Vector<Real>  binVals(nBins*nAvgComps*nBinsZ,0.0);
        Vector<Real>  binValsSq(nBins*nAvgComps*nBinsZ,0.0);
        Vector<Real>  binHits(nBins*nBinsZ,0.0);

        for (int iLevel=0; iLevel<Nlev; ++iLevel) {
            if (verbose) {
               Print() << " Reading data at level " << iLevel << "\n";
            }

            amrData.FillVar(mfVec[iLevel], iLevel, readInComps, destFillComps);

            // Prepare mask, zeroing out fine-covered and EB-covered data
            mfVec[iLevel].setVal(1.0,nComp-1,1); // initialize mask to value
            if (iLevel<thisFinestLevel) {
                int ratio = amrData.RefRatio()[iLevel];
                BoxArray baf = BoxArray(amrData.boxArray(iLevel+1)).coarsen(ratio);
                for (MFIter mfi(mfVec[iLevel]); mfi.isValid(); ++mfi) {
                    const Box& box = mfi.validbox();
                    FArrayBox& fab = mfVec[iLevel][mfi];
                    std::vector< std::pair<int,Box> > isects = baf.intersections(box);
                    for (int ii = 0; ii < isects.size(); ii++)
                        fab.setVal(0.0,isects[ii].second,nComp-1,1);
                }
            }

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            for (MFIter mfi(mfVec[iLevel]); mfi.isValid(); ++mfi) {
                const Box& bx = mfi.validbox();
                auto        vf_a   = mfVec[iLevel].const_array(mfi,0);
                auto const& mask_a = mfVec[iLevel].array(mfi,nComp-1);
                amrex::ParallelFor(bx, [vf_a,mask_a] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {   
                    if (vf_a(i,j,k) < 0.999) {
                       mask_a(i,j,k) = 0.0;
                    }
                });
            }
   
            // Compute mixture fraction and thermal diff
            if (verbose) {
               Print() << " Mixture fraction at level " << iLevel << "\n";
            }
            auto const* ltransparm = trans_parms.device_trans_parm();
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            for (MFIter mfi(mfVec[iLevel]); mfi.isValid(); ++mfi) {
                const Box& bx = mfi.validbox();
                auto const& Ys  = mfVec[iLevel].const_array(mfi,1);
                auto const& T   = mfVec[iLevel].const_array(mfi,NUM_SPECIES+1);
                auto const& rho = mfVec[iLevel].const_array(mfi,NUM_SPECIES+2);
                auto       mixt_frac = mixfracVec[iLevel].array(mfi,0);
                auto       diffus    = mixfracVec[iLevel].array(mfi,1);

                amrex::Real denom_inv = 1.0 / (mFrac_d.Zfu - mFrac_d.Zox);

                amrex::ParallelFor(bx, [Ys, T, rho, mixt_frac, denom_inv, mFrac_d, diffus, ltransparm, isPeleC]
                AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {   
                    mixt_frac(i,j,k) = 0.0_rt;
                    for (int n = 0; n<NUM_SPECIES; ++n) {
                        mixt_frac(i,j,k) += Ys(i,j,k,n) * mFrac_d.spec_Bilger_fact[n];
                    }   
                    mixt_frac(i,j,k) = ( mixt_frac(i,j,k) - mFrac_d.Zox ) * denom_inv;

                    // Thermal diffusivity
                    auto eos = pele::physics::PhysicsType::eos();
                    Real rho_cgs = rho(i,j,k);
                    if (!isPeleC) {
                        rho_cgs *= 1.0e-3;
                    }
                    amrex::Real temp = T(i,j,k);
                    amrex::Real Ysp[NUM_SPECIES] = {0.0};
                    for (int n = 0; n < NUM_SPECIES; n++) {
                        Ysp[n] = Ys(i,j,k,n);
                    }
                    amrex::Real dummy_rhoDi[NUM_SPECIES] = {0.0};
                    amrex::Real lambda_cgs = 0.0;
                    amrex::Real dummy_mu = 0.0;
                    amrex::Real dummy_xi = 0.0;

                    bool get_xi = false;
                    bool get_mu = false;
                    bool get_lam = true;
                    bool get_Ddiag = false;
                    auto trans = pele::physics::PhysicsType::transport();
                    trans.transport(get_xi, get_mu, get_lam, get_Ddiag, temp,
                                    rho_cgs, Ysp, dummy_rhoDi, dummy_mu, dummy_xi, lambda_cgs, ltransparm);

                    amrex::Real cpmix = 0.0;
                    eos.TY2Cp(temp, Ysp, cpmix);
                    if (!isPeleC) {
                        cpmix *= 1.0e-4;
                        lambda_cgs *= 1.0e-5; // it's actually MKS now
                    }
                    diffus(i,j,k) = lambda_cgs / rho(i,j,k) / cpmix;
                });
            }

            mixfracVec[iLevel].FillBoundary(0,1,geoms[iLevel].periodicity());
        }

        if (verbose) {
           Print() << " Getting scalar dissipation rate using MLMG \n";
        }

        // Get face-centered gradients from MLMG
        LPInfo info;
        info.setAgglomeration(1);
        info.setConsolidation(1);
        info.setMetricTerm(false);
        info.setMaxCoarseningLevel(0);
        MLPoisson poisson({geoms}, {grids}, {dmaps}, info);
        poisson.setMaxOrder(4);
        std::array<LinOpBCType, AMREX_SPACEDIM> lo_bc;
        std::array<LinOpBCType, AMREX_SPACEDIM> hi_bc;
        for (int idim = 0; idim< AMREX_SPACEDIM; idim++){
           if (is_per[idim] == 1) {
              lo_bc[idim] = hi_bc[idim] = LinOpBCType::Periodic;
           } else {
              lo_bc[idim] = hi_bc[idim] = LinOpBCType::Neumann;
           }
        }
        poisson.setDomainBC(lo_bc, hi_bc);

        // Need to apply the operator to ensure CF consistency with composite solve
        int nGrowGrad = 0;                   // No need for ghost face on gradient
        Vector<Array<MultiFab,AMREX_SPACEDIM> > grad(Nlev);
        Vector<std::unique_ptr<MultiFab>> phi;
        Vector<MultiFab> laps;
        for (int lev = 0; lev < Nlev; ++lev) {
          for (int idim = 0; idim <AMREX_SPACEDIM; idim++) {
             const auto& ba = grids[lev];
             grad[lev][idim].define(amrex::convert(ba,IntVect::TheDimensionVector(idim)),
                                    dmaps[lev], 1, nGrowGrad);
          }    
          phi.push_back(std::make_unique<MultiFab> (mixfracVec[lev],amrex::make_alias,0,1));
          poisson.setLevelBC(lev, phi[lev].get());
          laps.emplace_back(grids[lev], dmaps[lev], 1, 1);
        }

        MLMG mlmg(poisson);
        mlmg.apply(GetVecOfPtrs(laps), GetVecOfPtrs(phi));
        mlmg.getFluxes(GetVecOfArrOfPtrs(grad), GetVecOfPtrs(phi), MLMG::Location::FaceCenter);

        for (int lev = 0; lev < Nlev; ++lev) {
            // Get the scalar dissipation rate D \nabla mixFrac \cdot \nabla mixFrac
            MultiFab gradAlias(mixfracVec[lev], amrex::make_alias, 2, AMREX_SPACEDIM);
            average_face_to_cellcenter(gradAlias, 0, GetArrOfConstPtrs(grad[lev]));
            gradAlias.mult(-1.0);
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            for (MFIter mfi(mfVec[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {    
               const Box& bx = mfi.tilebox();
               auto const& grad_a   = gradAlias.const_array(mfi);
               auto const& diffus   = mixfracVec[lev].array(mfi,1);
               auto const& scaldiss = mfVec[lev].array(mfi,nComp-2);
               amrex::ParallelFor(bx, [=]
               AMREX_GPU_DEVICE (int i, int j, int k) noexcept
               {    
                  scaldiss(i,j,k) = diffus(i,j,k) * (AMREX_D_TERM(  grad_a(i,j,k,0) * grad_a(i,j,k,0),
                                                                  + grad_a(i,j,k,1) * grad_a(i,j,k,1),
                                                                  + grad_a(i,j,k,2) * grad_a(i,j,k,2)));
               });  
            } 
        }

        if (writePltMixFrac) {
            Vector<std::unique_ptr<MultiFab>> scalDiss;
            for (int lev = 0; lev < Nlev; ++lev) {
              scalDiss.push_back(std::make_unique<MultiFab> (mfVec[lev],amrex::make_alias,nComp-2,1));
            }
            Vector<int> isteps(Nlev, 0);
            Vector<IntVect> refRatios(Nlev-1,{AMREX_D_DECL(2, 2, 2)});
            amrex::WriteMultiLevelPlotfile("ScalDiss", Nlev, GetVecOfConstPtrs(scalDiss),{"ScalDiss"},
                                           geoms, 0.0, isteps, refRatios);
            amrex::WriteMultiLevelPlotfile("MixFrac", Nlev, GetVecOfConstPtrs(phi),{"MixtureFraction"},
                                           geoms, 0.0, isteps, refRatios);
        }

        for (int iLevel=0; iLevel<Nlev; ++iLevel) {

            MultiFab vol(grids[iLevel], dmaps[iLevel], 1, ngrow);
            geoms[iLevel].GetVolume(vol);

            const auto geomdata = geoms[iLevel].data();

            if (verbose) {
               Print() << " CM at level " << iLevel << "\n";
            }
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            for(MFIter mfi(mfVec[iLevel]); mfi.isValid(); ++mfi) {
                const Box& box = mfi.validbox();
                FArrayBox& fab = mfVec[iLevel][mfi];
                FArrayBox& vfab = vol[mfi];
                FArrayBox& mixfracfab = mixfracVec[iLevel][mfi];
                for (IntVect iv=box.smallEnd(); iv<=box.bigEnd(); box.next(iv)) {
                    int myZBin = 0;
                    if (doZcoorBins) {
                        Real zcoor = geomdata.ProbLo(2) + (static_cast<Real>(iv[2])+0.5)*geomdata.CellSize(2);
                        if (zcoor < binZMin || zcoor > binZMax ) {
                           continue;
                        }
                        Real binZwidth_inv = static_cast<Real>(nBinsZ)/(binZMax-binZMin);
                        myZBin = std::floor((zcoor - binZMin) * binZwidth_inv);
                    }
                    Vector<Real> voldata(1,0.0);
                    Vector<Real> mixfracdata(AMREX_SPACEDIM+2,0.0);
                    vfab.getVal(voldata.dataPtr(),iv);
                    mixfracfab.getVal(mixfracdata.dataPtr(),iv);
                    fab.getVal(dataIV.dataPtr(),iv);
                    if (dataIV[nComp-1] > 0.0) {   // Mask fine/EB covered
                        Real binVal = mixfracdata[0];
                        if (binVal>=binMin && binVal<binMax) {
                            int myVBin =   (int)(nBins*(binVal-binMin)/(binMax-binMin));
                            if (myVBin<0 || myVBin>=nBins) {
                                amrex::Abort("Bad bin");
                            }
                            Real myWeight = voldata[0];
                            for (int f{0}; f<nAvgComps; ++f) {
                                Real val = dataIV[compIdx[f]];
                                int binIdx = myZBin*nBins*nAvgComps + myVBin*nAvgComps + f;
                                binVals[binIdx]    += myWeight * val;
                                binValsSq[binIdx]  += myWeight * val * val;
                            }
                            int hbinIdx = myZBin*nBins + myVBin;
                            binHits[hbinIdx] += myWeight;
                        }
                    }
                }
            } // MFI
        } // Level

        ParallelDescriptor::ReduceRealSum(binHits.dataPtr(),binHits.size(),ioproc);
        ParallelDescriptor::ReduceRealSum(binVals.dataPtr(),binVals.size(),ioproc);
        ParallelDescriptor::ReduceRealSum(binValsSq.dataPtr(),binValsSq.size(),ioproc);
        if (isioproc) {
            for (int zB{0}; zB<nBinsZ; ++zB) {
                std::string filename = "CondMean_Time" + std::to_string(amrData.Time());
                if (nBinsZ>1) {
                    const Real dz = (binZMax - binZMin)/static_cast<Real>(nBinsZ);
                    Real zCent = binZMin + (0.5+static_cast<Real>(zB))*dz;
                    filename += "_Zcent"+std::to_string(zCent);
                }
                filename += ".dat";
                std::ofstream ofs(filename.c_str());
                std::string variables = "VARIABLES = mixFrac ";
                ofs << variables;
                for (int f{0}; f<nAvgComps; ++f) {
                    const std::string varName = outComps[f];
                    ofs << varName+"_avg " << varName+"_stddev " << varName+"_Int ";
                }
                ofs << "\n";
                const Real dv = (binMax - binMin)/nBins;
                for (int n=0; n<nBins; n++) {
                    const Real v = binMin + dv*(0.5+(Real)n);
                    ofs << v << " ";
                    for (int f{0}; f<nAvgComps; ++f) {
                        int hbinIdx = zB*nBins + n;
                        if (binHits[hbinIdx] > 0.0) {
                            int binIdx = zB*nBins*nAvgComps+n*nAvgComps+f;
                            Real bh = binHits[hbinIdx];
                            ofs << binVals[binIdx]/bh << " ";
                            ofs << std::sqrt( (binValsSq[binIdx]/bh) - (binVals[binIdx]/bh)*(binVals[binIdx]/bh) ) << " ";
                            ofs << binVals[binIdx] << " ";
                        } else {
                            ofs << " 0.0 0.0 0.0 ";
                        }
                    }
                    ofs << "\n";
                }
                ofs.close();
            }
        }

        const Real end_plt_time_io = ParallelDescriptor::second();
        Print() << "Treating plot time: " << end_plt_time_io - plt_time_io << '\n';
    } // iPlot
        

//    // Output result
//    if (isioproc)
//    {            
//        std::string filename;
//        // Output data for tecplot
//        //filename = infile + outSuffix + "/CM_" + compNames[0] + ".dat";
//        if (aja) {
//            filename = plotFileNames[0] + "/CM_" + compNames[0] + ".key";
//            std::cout << "Opening file " << filename << std::endl;
//        } else {
//            // Default
//            filename = "CM_" + compNames[0] + ".dat";
//            std::cout << "Opening file " << filename << std::endl;
//        }
//        std::ofstream ofs(filename.c_str());
//        std::string variables = "VARIABLES = " + compNames[0];
//        for (int i=1; i<compNames.size(); ++i)
//            variables += " " + compNames[i] + "_sum";
//        for (int i=1; i<compNames.size(); ++i)
//            variables += " " + compNames[i] + "_sumSq";
//        for (int i=1; i<compNames.size(); ++i)
//            variables += " " + compNames[i] + "_avg";
//        for (int i=1; i<compNames.size(); ++i)
//            variables += " " + compNames[i] + "_std";
//        if (writeBinMinMax)
//        {
//            for (int i=1; i<compNames.size(); ++i)
//                variables += " " + compNames[i] + "_min";
//            for (int i=1; i<compNames.size(); ++i)
//                variables += " " + compNames[i] + "_max";
//        }
//        variables += " N ";
//        variables += " p ";
//        variables += '\n';
//        ofs << variables.c_str();
//        ofs << "ZONE I=" << nBins << " DATAPACKING=POINT\n";
//        if (aja) {
//            ofs.close();
//            filename = plotFileNames[0] + "/CM_" + compNames[0] + ".dat";
//            std::cout << "Opening file " << filename << std::endl;
//            ofs.open(filename.c_str());
//        }
//        const Real dv = (binMax - binMin)/nBins;
//        int ntot = 0;
//        for (int i=0; i<nBins; i++)
//        {
//            ntot += binHits[i];
//        }
//        for (int i=0; i<nBins; i++) {
//            const Real v = binMin + dv*(0.5+(Real)i);
//            ofs << v << " ";
//            // Sum
//            for (int j=0; j<nAvgComps; j++) {
//                ofs << binVals[i*nAvgComps + j] << " ";
//            }
//            // SumSq
//            for (int j=0; j<nAvgComps; j++) {
//                ofs << binValsSq[i*nAvgComps + j] << " ";
//            }
//       if (binHits[i]>0) {
//      // Avg
//      for (int j=0; j<nAvgComps; j++) {
//          ofs << binVals[i*nAvgComps + j]/(Real)binHits[i] << " ";
//      }
//      // Std Dev
//      for (int j=0; j<nAvgComps; j++) {
//          int idx = i*nAvgComps + j;
//          Real bh = (Real)binHits[i];
//                    ofs << std::sqrt( (binValsSq[idx]/bh) - (binVals[idx]/bh)*(binVals[idx]/bh) ) << " ";
//      }
//       } else {
//      // Avg & Std Dev
//      for (int j=0; j<nAvgComps*2; j++) {
//                    ofs << "0.0 ";
//      }
//       }
//       if (writeBinMinMax)
//            {
//                for (int j=0; j<nAvgComps; j++) {
//                    ofs << binMinVals[i*nAvgComps + j] << " ";
//                }
//                for (int j=0; j<nAvgComps; j++) {
//                    ofs << binMaxVals[i*nAvgComps + j] << " ";
//                }
//            }
//       ofs << Real(binHits[i]) << " ";
//            ofs << Real(binHits[i])/ntot << '\n';
//        }
//        cout << "total bins: " << ntot << endl;
//        ofs.close();
//        
//    } // IOProcessor
    
    amrex::Finalize();
    return 0;
}

