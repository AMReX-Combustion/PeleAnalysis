#include <string>
#include <iostream>
using std::cerr;
using std::cout;
using std::endl;

#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include <AMReX_DataServices.H>

using namespace amrex;

static std::string
ProtectSlashes(const std::string str)
{
  std::string s = str;

  std::vector<int> where;

  const char* cstr = s.c_str();

  for (int i = 0; i < s.size(); i++)
    if (cstr[i] == '/')
      where.push_back(i);

  for (int i = where.size() - 1; i >= 0; i--)
    s.replace(where[i], 1, "_");

  return s;
}

static void
print_usage(int, char* argv[])
{
  std::cerr << "usage:\n";
  std::cerr
    << argv[0] << " infile=<filename>\n"
    << "   output_gnuplot=[0/1]    : Output JPDFs for gnuplot             "
       "(default=0)\n"
    << "   output_matlab =[0/1]    : Output JPDFs for matlab              "
       "(default=0)\n"
    << "   output_tecplot=[0/1]    : Output JPDFs for tecplot             "
       "(default=0)\n"
    << "   output_scatter=[0/1]    : Output scatter plot of non-zero bins "
       "(default=0)\n"
    << "   output_fab=[0/1]        : Output JPDFs in fab format           "
       "(default=0)\n"
    << "   output_plotfile=[0/1]   : Output JPDFs in plotfile format      "
       "(default=1)\n"
    << "   vars= var1 var2 var3    : Variable list; JPDF will be evaluated for "
       "each pair\n"
    << "   useminmax%i= min max    : Use min/max for var %i\n"
    << "   condMean_vars = var     : Variable that will be averaged "
       "conditioned on vars\n"
    << "   nBins=n                 : Number of bins in each direction for JPDF "
       "(default=64)\n"
    << "   finestLevel=n           : Finest level at which to evaluate JPDFs\n"
    << "   outSuffix=str           : Suffix to add to the plfile name as an "
       "alt dir for results (default="
       ")\n"
    << "   do_conditioning=[0/1/2] : 0 -> no conditioning (default)\n"
    << "                             1 -> conditioning on a loaded variable \n"
    << "                             2 ->  Evaluation with c(1-c)\n"
    << "     cVar=[0...nVars]      : Index in 'vars' of var for conditioning "
       "(default=0)\n"
    << "     cMin = min            : Min value for conditioning (default=0.0)\n"
    << "     cMax = max            : Max value for conditioning (default=1.0)\n"
    << "     norm_cVal=[0/1]       : Normalize cVar (default=0, forced for "
       "do_conditioning=2)\n"
    << "        cNormMin = min     : Min value for normalization "
       "(default=0.0)\n"
    << "        cNormMax = max     : Max value for normalization "
       "(default=1.0)\n"
    << "   do_average=[0/1]        : Performs an average over all provided "
       "plot files (default=0)\n"
    << "   infile= plt1 plt2       : List of plot files" << std::endl;
  exit(1);
}

int
main(int argc, char* argv[])
{
  Initialize(argc, argv);
  {
    if (argc < 2)
      print_usage(argc, argv);

    ParmParse pp;

    if (pp.contains("help"))
      print_usage(argc, argv);

    if (pp.contains("verbose"))
      AmrData::SetVerbose(true);

    bool verbose = false;
    if (ParallelDescriptor::IOProcessor())
      verbose = true;

    // Get output file types
    if (verbose)
      std::cout << "Output types:" << std::endl;
    int output_gnuplot(0);
    pp.query("output_gnuplot", output_gnuplot);
    if (output_gnuplot && verbose)
      std::cout << "   + gnuplot" << std::endl;
    int output_matlab(0);
    pp.query("output_matlab", output_matlab);
    if (output_matlab && verbose)
      std::cout << "   + matlab" << std::endl;
    int output_tecplot(0);
    pp.query("output_tecplot", output_tecplot);
    if (output_tecplot && verbose)
      std::cout << "   + tecplot" << std::endl;
    int output_fab(0);
    pp.query("output_fab", output_fab);
    if (output_fab && verbose)
      std::cout << "   + fab" << std::endl;
    int output_plotfile(1);
    pp.query("output_plotfile", output_plotfile);
    if (output_plotfile && verbose)
      std::cout << "   + plotfile" << std::endl;
    int output_scatter(0);
    pp.query("output_scatter", output_scatter);
    if (output_scatter && verbose)
      std::cout << "   + scatter" << std::endl;

    // Do conditioning
    int do_conditioning(0);
    int cVar(0);
    int norm_cVal(0);
    Real cMin(0.0);
    Real cMax(1.0);
    Real cNormMin(0.0);
    Real cNormMax(1.0);

    pp.query("do_conditioning", do_conditioning);
    if (verbose)
      std::cout << "do_conditioning = " << do_conditioning << std::endl;

    if (do_conditioning > 0) {
      pp.query("cVar", cVar);
      if (verbose)
        std::cout << "cVar = " << cVar << std::endl;

      pp.query("norm_cVal", norm_cVal);
      if (verbose)
        std::cout << "norm_cVal = " << norm_cVal << std::endl;

      if (do_conditioning == 2)
        norm_cVal = 1; // Force normalisation for c(1-c) evaluation

      if (norm_cVal == 1) { // Get min max for c normalisation
        pp.query("cNormMin", cNormMin);
        if (verbose)
          std::cout << "cNormMin = " << cNormMin << std::endl;
        pp.query("cNormMax", cNormMax);
        if (verbose)
          std::cout << "cNormMax = " << cNormMax << std::endl;
      }

      pp.query("cMin", cMin);
      if (verbose)
        std::cout << "cMin = " << cMin << std::endl;

      pp.query("cMax", cMax);
      if (verbose)
        std::cout << "cMax = " << cMax << std::endl;
    }

    int do_average(0);
    pp.query("do_average", do_average);
    if (do_average && verbose)
      std::cout << "   + forming average" << std::endl;

    // Get plot file
    int nPlotFiles(pp.countval("infile"));
    if (nPlotFiles <= 0) {
      std::cerr << "Bad nPlotFiles:  " << nPlotFiles << std::endl;
      std::cerr << "Exiting." << std::endl;
      DataServices::Dispatch(DataServices::ExitRequest, NULL);
    }
    if (verbose)
      std::cout << "Processing " << nPlotFiles << " plotfiles..." << std::endl;

    std::string outSuffix = "";
    pp.query("outSuffix", outSuffix);

    // Make an array of srings containing paths of input plot files
    Vector<std::string> plotFileNames(nPlotFiles);
    for (int iPlot = 0; iPlot < nPlotFiles; ++iPlot) {
      pp.get("infile", plotFileNames[iPlot], iPlot);
      if (verbose)
        std::cout << "   " << plotFileNames[iPlot] << std::endl;
    }

    // Get finest level argument
    int inFinestLevel(-1);
    pp.query("finestLevel", inFinestLevel);

    // Number of bins in each direction
    int nBins(64);
    pp.query("nBins", nBins);

    // Reading variables ...
    // Variables for jPDF; a joint pdf of each pair will be evaluated
    int nVars_jpdf(pp.countval("vars"));
    if (nVars_jpdf < 2)
      Error("Need to specify at least two variables.");

    Vector<std::string> whichVar_jPDF(nVars_jpdf);
    if (verbose)
      std::cout << "Variable list:" << std::endl;
    // Read in variable list
    for (int v = 0; v < nVars_jpdf; v++) {
      pp.get("vars", whichVar_jPDF[v], v);
      if (verbose)
        std::cout << "   " << whichVar_jPDF[v] << std::endl;
    }

    // Variables for conditional mean
    int nVars_condMean(pp.countval("condMean_vars"));
    int nVars_tot(nVars_jpdf); // Total num of varables is jPDF vars plus unique
                               // condMean vars
    Vector<std::string> whichVar_condMean(nVars_condMean);
    if (verbose)
      std::cout << "Variable list for conditional mean:" << std::endl;
    // Read in variable list
    for (int v = 0; v < nVars_condMean; v++) {
      pp.get("condMean_vars", whichVar_condMean[v], v);
      if (verbose)
        std::cout << "   " << whichVar_condMean[v] << std::endl;
      // Check if variable is contained in whichVar_jPDF
      int is_unique(1);
      for (int k = 0; k < nVars_jpdf; k++) {
        if (whichVar_condMean[v] == whichVar_jPDF[k]) {
          is_unique = 0;
          break;
        }
      }
      nVars_tot += is_unique;
    }

    // Create a set of jPDF vars plus unique condMean vars; save the indices of
    // the condMeanVars is ix_condMean
    Vector<std::string> whichVar(nVars_tot);
    Vector<int> ix_condMean(nVars_condMean);
    if (verbose)
      std::cout << "Combined variable list:" << std::endl;
    for (int v = 0; v < nVars_jpdf; v++) {
      whichVar[v] = whichVar_jPDF[v];
      if (verbose)
        std::cout << "   " << v << " : " << whichVar_jPDF[v] << std::endl;
    }
    int ix_fill(nVars_jpdf);
    for (int v = 0; v < nVars_condMean; v++) {
      int dupplicat_ix(-1);
      for (int k = 0; k < nVars_jpdf; k++) {
        if (whichVar_condMean[v] == whichVar_jPDF[k]) {
          dupplicat_ix = k;
          break;
        }
      }
      if (dupplicat_ix == -1) {
        whichVar[ix_fill] = whichVar_condMean[v];
        ix_condMean[v] = ix_fill;
        if (verbose)
          std::cout << "   " << ix_condMean[v] << "=" << ix_fill << " : "
                    << whichVar[ix_fill] << std::endl;
        ix_fill++;
      } else {
        ix_condMean[v] = dupplicat_ix;
        if (verbose)
          std::cout << "   -- condMeanVar '" << whichVar_condMean[v]
                    << "' already on position " << ix_condMean[v] << std::endl;
      }
    }

    // Intersect variable number (nVars_tot does not include intersect variable)
    int ix_isVar = nVars_tot;

    // Copy the names of the variable for output filenames, replacing dodgy
    // characters
    Vector<std::string> whichVarOut(nVars_tot);
    for (int iVar = 0; iVar < nVars_tot; iVar++) {
      whichVarOut[iVar] = ProtectSlashes(whichVar[iVar]);
      if (verbose)
        std::cout << whichVar[iVar] << " -> " << whichVarOut[iVar] << std::endl;
    }

    // Initialise some stuff here for average later
    int nPairs(nVars_jpdf * (nVars_jpdf - 1) / 2);
    Real domainVol;
    Vector<Vector<Real>> binAvArray(nPairs);
    Vector<Vector<Real>> binAvX1Array(nPairs);
    Vector<Vector<Real>> binAvX2Array(nPairs);
    Vector<Vector<Vector<Real>>> binAvCondMeanArray(
      nPairs,
      Vector<Vector<Real>>(nVars_condMean, Vector<Real>(nBins * nBins, 0)));
    Vector<Real> vMin(nVars_jpdf);
    Vector<Real> vMax(nVars_jpdf);

    // Loop over files
    for (int iPlot = 0; iPlot < nPlotFiles; iPlot++) {

      // Open file and get an amrData pointer
      std::string infile = plotFileNames[iPlot];
      if (verbose)
        std::cout << "\nOpening " << infile << "..." << std::endl;
      DataServices::SetBatchMode();
      Amrvis::FileType fileType(Amrvis::NEWPLT);
      DataServices dataServices(infile, fileType);
      if (!dataServices.AmrDataOk())
        DataServices::Dispatch(DataServices::ExitRequest, NULL);
      AmrData& amrData = dataServices.AmrDataRef();
      if (verbose)
        std::cout << "   ...done." << std::endl;

      // Check the names of the variables are present in the plotfile
      Vector<int> destFills(nVars_tot);
      for (int v = 0; v < nVars_tot; v++) {
        destFills[v] = v;
        if (amrData.StateNumber(whichVar[v]) < 0) {
          std::string message = "Bad variable name (" + whichVar[v] + ")";
          Error(message.c_str());
        }
      }

      // Get finest level
      int finestLevel = amrData.FinestLevel();
      if (inFinestLevel > -1 && inFinestLevel < finestLevel) {
        finestLevel = inFinestLevel;
        if (verbose)
          std::cout << "Finest level: " << finestLevel << std::endl;
      }
      int nLevels = finestLevel + 1;

      // Get domain size
      Vector<Real> probLo = amrData.ProbLo();
      Vector<Real> probHi = amrData.ProbHi();

      // Get probDomain
      Vector<Box> probDomain = amrData.ProbDomain();

      // Get min/max for each component
      for (int iVar = 0; iVar < nVars_jpdf;
           iVar++) { // The variable to average is not included here
        vMin[iVar] = 1e100;
        vMax[iVar] = -vMin[iVar];
        for (int iLevel = 0; iLevel < nLevels; iLevel++) {
          Real min, max;
          amrData.MinMax(probDomain[iLevel], whichVar[iVar], iLevel, min, max);
          if (vMin[iVar] > min)
            vMin[iVar] = min;
          if (vMax[iVar] < max)
            vMax[iVar] = max;
        }
      }

      for (int iVar = 0; iVar < nVars_jpdf;
           iVar++) { // The variable to average is not included here
        char argName[12];
        sprintf(argName, "useminmax%i", iVar + 1);
        int nMinMax = pp.countval(argName);
        if (nMinMax > 0) {
          if (nMinMax != 2) {
            Abort("Need to specify 2 values for useMinMax");
          } else {
            pp.get(argName, vMin[iVar], 0);
            pp.get(argName, vMax[iVar], 1);
            if (verbose)
              std::cout << "Var" << iVar + 1 << " (" << whichVar[iVar]
                        << ") using min/max: " << vMin[iVar] << " / "
                        << vMax[iVar] << std::endl;
          }
        }
      }

      // Declare space for a joint pdf for each variable pair
      Vector<Vector<Real>> binArray(nPairs);
      Vector<Vector<Real>> binX1Array(nPairs);
      Vector<Vector<Real>> binX2Array(nPairs);

      for (int iPair = 0; iPair < nPairs; iPair++) {
        binArray[iPair].resize(nBins * nBins, 0);
        binX1Array[iPair].resize(nBins * nBins, 0);
        binX2Array[iPair].resize(nBins * nBins, 0);
      }

      // Vector for the conditional means
      // for each pair, nVars_condMean vectors with nBins^2 bins each are
      // prepared
      Vector<Vector<Vector<Real>>> binCondMeanArray(
        nPairs,
        Vector<Vector<Real>>(nVars_condMean, Vector<Real>(nBins * nBins, 0)));

      if (do_average && iPlot == 0) {
        for (int iPair = 0; iPair < nPairs; iPair++) {
          binAvArray[iPair].resize(nBins * nBins, 0);
          binAvX1Array[iPair].resize(nBins * nBins, 0);
          binAvX2Array[iPair].resize(nBins * nBins, 0);
        }
      }

      // Load the data
      if (verbose)
        std::cout << "Loading data..." << std::endl;
      Vector<std::unique_ptr<MultiFab>> mf(nLevels);
      for (int iLevel = 0; iLevel < nLevels; iLevel++) {
        if (verbose)
          std::cout << "   Level " << iLevel << "..." << std::endl;
        int ngrow(0);

        // Make space for each component and one for the intersect flag
        DistributionMapping dm(amrData.boxArray(iLevel));
        mf[iLevel].reset(
          new MultiFab(amrData.boxArray(iLevel), dm, nVars_tot + 1, ngrow));

        // Put the value 1 in the intersect flag
        mf[iLevel]->setVal(1., ix_isVar, 1);

        // Load the data (make a copy of whichVar to make sure it's the same
        // size as destFills)
        Vector<std::string> loadWhichVar(nVars_tot);
        for (int v = 0; v < nVars_tot; v++)
          loadWhichVar[v] = whichVar[v];
        amrData.FillVar(*mf[iLevel], iLevel, loadWhichVar, destFills);
      }
      if (verbose)
        std::cout << "      ...done." << std::endl;

      // Set the intersect flag to zero if necessary
      for (int iLevel = 0; iLevel < finestLevel; iLevel++) {
        BoxArray baf = mf[iLevel + 1]->boxArray();
        baf.coarsen(amrData.RefRatio()[iLevel]);

        for (MFIter mfi(*mf[iLevel]); mfi.isValid(); ++mfi) {
          FArrayBox& myFab = (*mf[iLevel])[mfi];

          int idx = mfi.index();

          std::vector<std::pair<int, Box>> isects =
            baf.intersections(mf[iLevel]->boxArray()[idx]);

          for (int ii = 0; ii < isects.size(); ii++)
            myFab.setVal(0, isects[ii].second, ix_isVar, 1);
        }
      }

      // Evaluate the pdf for each pair
      if (verbose)
        std::cout << "Evaluating pdfs..." << std::endl;
      int iPair = 0;
      for (int var1 = 0; var1 < nVars_jpdf; var1++) {
        for (int var2 = var1 + 1; var2 < nVars_jpdf; var2++) {
          if (verbose)
            std::cout << "   + " << whichVar[var1] << "-" << whichVar[var2]
                      << std::endl;
          Real* bin = binArray[iPair].dataPtr();
          Real* binX1 = binX1Array[iPair].dataPtr();
          Real* binX2 = binX2Array[iPair].dataPtr();
          Vector<Real*> binCondMeans(nVars_condMean);
          for (int v = 0; v < nVars_condMean; v++) {
            binCondMeans[v] = binCondMeanArray[iPair][v].dataPtr();
          }
          Real *binAv, *binAvX1, *binAvX2;
          Vector<Real*> binAvCondMeans(nVars_condMean);
          if (do_average) {
            binAv = binAvArray[iPair].dataPtr();
            binAvX1 = binAvX1Array[iPair].dataPtr();
            binAvX2 = binAvX2Array[iPair].dataPtr();
            for (int v = 0; v < nVars_condMean; v++) {
              binAvCondMeans[v] = binAvCondMeanArray[iPair][v].dataPtr();
            }
          }
          for (int iLevel = 0; iLevel < nLevels; iLevel++) {
            if (verbose)
              std::cout << "      Level " << iLevel << std::endl;
            int v1l = 0;
            int v1g = 0;
            int v2l = 0;
            int v2g = 0;
            for (MFIter ntmfi(*mf[iLevel]); ntmfi.isValid(); ++ntmfi) {
              const FArrayBox& myFab = (*mf[iLevel])[ntmfi];
              const Real* dx = amrData.DxLevel()[iLevel].dataPtr();
              const Real* cPtr = myFab.dataPtr(cVar); // Conditioning
              const Real* v1Ptr = myFab.dataPtr(var1);
              const Real* v2Ptr = myFab.dataPtr(var2);
              Vector<const Real*> condVarPtrs(nVars_condMean);
              for (int v = 0; v < nVars_condMean; v++) {
                condVarPtrs[v] = myFab.dataPtr(ix_condMean[v]);
              }
              const Real* isPtr = myFab.dataPtr(ix_isVar); // Intersect

              const Box& bx = ntmfi.validbox();
              const int* lo = bx.loVect();
              const int* hi = bx.hiVect();

              const int ix = hi[0] - lo[0] + 1;
              const int jx = hi[1] - lo[1] + 1;
              Real Vol = dx[0] * dx[1];
#if (AMREX_SPACEDIM == 3)
              const int kx = hi[2] - lo[2] + 1;
              ;
              Vol *= dx[2];
              for (int k = 0; k < kx; k++) {
                Real z = probLo[2] + dx[2] * (0.5 + (Real)(k + lo[2]));
#endif
                for (int j = 0; j < jx; j++) {
                  Real y = probLo[1] + dx[1] * (0.5 + (Real)(j + lo[1]));
                  for (int i = 0; i < ix; i++) {
                    Real x = probLo[0] + dx[0] * (0.5 + (Real)(i + lo[0]));
#if (AMREX_SPACEDIM == 3)
                    int cell = (k * jx + j) * ix + i;
#else
                  int cell = j * ix + i;
#endif
                    if (isPtr[cell] > 0) { // Intersect

                      int contribute = 1;

                      if (do_conditioning > 0) {

                        double cVal = cPtr[cell];

                        if (norm_cVal == 1)
                          cVal = (cVal - cNormMin) / (cNormMax - cNormMin);

                        if (do_conditioning == 2)
                          cVal = cVal * (1. - cVal);

                        if (cVal < cMin || cVal > cMax)
                          contribute = 0;
                      }

                      if (contribute == 1) {
                        int v1i = (int)(nBins * (v1Ptr[cell] - vMin[var1]) /
                                        (vMax[var1] - vMin[var1]));
                        if (v1i < 0) {
                          v1l++;
                          v1i = 0;
                        }
                        if (v1i >= nBins) {
                          v1g++;
                          v1i = nBins - 1;
                        }
                        int v2i = (int)(nBins * (v2Ptr[cell] - vMin[var2]) /
                                        (vMax[var2] - vMin[var2]));
                        if (v2i < 0) {
                          v2l++;
                          v2i = 0;
                        }
                        if (v2i >= nBins) {
                          v2g++;
                          v2i = nBins - 1;
                        }
                        bin[v1i * nBins + v2i] += Vol;
                        binX1[v1i * nBins + v2i] += Vol * v1Ptr[cell];
                        binX2[v1i * nBins + v2i] += Vol * v2Ptr[cell];
                        for (int v = 0; v < nVars_condMean; v++) {
                          binCondMeans[v][v1i * nBins + v2i] +=
                            Vol * condVarPtrs[v][cell];
                        }
                        if (do_average) {
                          binAv[v1i * nBins + v2i] += Vol;
                          binAvX1[v1i * nBins + v2i] += Vol * v1Ptr[cell];
                          binAvX2[v1i * nBins + v2i] += Vol * v2Ptr[cell];
                          for (int v = 0; v < nVars_condMean; v++) {
                            binAvCondMeans[v][v1i * nBins + v2i] +=
                              Vol * condVarPtrs[v][cell];
                          }
                        }
                      }
                    }
                  } // i
                } // j
#if (AMREX_SPACEDIM == 3)
              } // k
#endif
            } // MFI
            ParallelDescriptor::ReduceIntSum(v1l);
            ParallelDescriptor::ReduceIntSum(v1g);
            ParallelDescriptor::ReduceIntSum(v2l);
            ParallelDescriptor::ReduceIntSum(v2g);
            if (verbose) {
              if (v1l)
                std::cout << "v1i<0:      " << v1l << std::endl;
              if (v1g)
                std::cout << "v1i>=nBins: " << v1g << std::endl;
              if (v2l)
                std::cout << "v2i<0:      " << v2l << std::endl;
              if (v2g)
                std::cout << "v2i>=nBins: " << v2g << std::endl;
            }
          } // Level
          iPair++;
        }
      }
      if (verbose)
        std::cout << "   ...done." << std::endl;

      // Reduce data
      for (int iPair = 0; iPair < nPairs; iPair++) {
        Real* bin = binArray[iPair].dataPtr();
        Real* binX1 = binX1Array[iPair].dataPtr();
        Real* binX2 = binX2Array[iPair].dataPtr();
        ParallelDescriptor::ReduceRealSum(
          bin, binArray[iPair].size(), ParallelDescriptor::IOProcessorNumber());
        ParallelDescriptor::ReduceRealSum(
          binX1, binX1Array[iPair].size(),
          ParallelDescriptor::IOProcessorNumber());
        ParallelDescriptor::ReduceRealSum(
          binX2, binX2Array[iPair].size(),
          ParallelDescriptor::IOProcessorNumber());

        for (int v = 0; v < nVars_condMean; v++) {
          Real* binCondMean = binCondMeanArray[iPair][v].dataPtr();
          ParallelDescriptor::ReduceRealSum(
            binCondMean, binCondMeanArray[iPair][v].size(),
            ParallelDescriptor::IOProcessorNumber());
        }
      }

      // Output the data to file
      if (ParallelDescriptor::IOProcessor()) {

        if (outSuffix != "") {
          std::string oFile = infile + outSuffix;
          if (!UtilCreateDirectory(oFile, 0755))
            CreateDirectoryFailed(oFile);
        }

        domainVol = 1;
        for (int dd = 0; dd < AMREX_SPACEDIM; ++dd) {
          domainVol *= probHi[dd] - probLo[dd];
        }

        Real small = 1.e-7;

        int iPair = 0;

        for (int var1 = 0; var1 < nVars_jpdf; var1++) {
          // Change in var1 between each bin
          Real dv1 = (vMax[var1] - vMin[var1]) / (Real)nBins;

          for (int var2 = var1 + 1; var2 < nVars_jpdf; var2++) {
            // Change in var2 between each bin
            Real dv2 = (vMax[var2] - vMin[var2]) / (Real)nBins;

            // Pull out the pointer for this pair
            Real* bin = binArray[iPair].dataPtr();
            Real* binX1 = binX1Array[iPair].dataPtr();
            Real* binX2 = binX2Array[iPair].dataPtr();

            Vector<Real*> binCondMeans(nVars_condMean);
            for (int v = 0; v < nVars_condMean; v++) {
              binCondMeans[v] = binCondMeanArray[iPair][v].dataPtr();
            }

            // Divide Xi by bin vol to average (has to come first)
            for (int v1i = 0, i = 0; v1i < nBins; v1i++) {
              Real v1 = vMin[var1] + dv1 * (0.5 + (Real)v1i);
              for (int v2i = 0; v2i < nBins; v2i++, i++) {
                Real v2 = vMin[var2] + dv2 * (0.5 + (Real)v2i);
                Real div = bin[i];
                if (div > 0) {
                  binX1[i] /= div;
                  binX2[i] /= div;
                  for (int v = 0; v < nVars_condMean; v++) {
                    binCondMeans[v][i] /= div;
                  }
                } else { // This is the bit that breaks the derivation of the 1D
                         // pdfs
                  binX1[i] = v1;
                  binX2[i] = v2;
                }
              }
            }

            // Normalise so that sum(bin) = 1. When conditioning is active the
            // total contributed volume is less than domainVol, so divide by the
            // actual accumulated volume rather than the full domain volume.
            Real normVol = domainVol;
            if (do_conditioning > 0) {
              Real totalVol = 0;
              for (int i = 0; i < nBins * nBins; i++)
                totalVol += bin[i];
              if (totalVol > 0)
                normVol = totalVol;
            }
            for (int i = 0; i < nBins * nBins; i++)
              bin[i] /= normVol;

            // Let's write some files...
            std::string filename;
            FILE* file;

            if (output_gnuplot) {
              // Output data for gnuplot
              // Format: x y pdf
              filename = infile + outSuffix + "/Pdf_" + whichVarOut[var1] +
                         "_" + whichVarOut[var2] + ".gpd";
              std::cout << "Opening file " << filename << std::endl;
              file = fopen(filename.c_str(), "w");
              for (int v1i = 0; v1i < nBins; v1i++) {
                Real v1 = vMin[var1] + dv1 * (0.5 + (Real)v1i);
                for (int v2i = 0; v2i < nBins; v2i++) {
                  Real v2 = vMin[var2] + dv2 * (0.5 + (Real)v2i);
                  fprintf(file, "%e %e %e\n", v1, v2, bin[v1i * nBins + v2i]);
                }
              }
              fclose(file);
            }

            if (output_matlab) {
              // Output data for matlab (5 files)
              // Format: pdf matrix
              filename = infile + outSuffix + "/Pdf_" + whichVarOut[var1] +
                         "_" + whichVarOut[var2] + ".dat";
              std::cout << "Opening file " << filename << std::endl;
              file = fopen(filename.c_str(), "w");
              for (int v1i = 0; v1i < nBins; v1i++) {
                for (int v2i = 0; v2i < nBins; v2i++)
                  fprintf(file, "%e ", bin[v1i * nBins + v2i]);
                fprintf(file, "\n");
              }
              fclose(file);
              // Format: var1 range
              filename =
                infile + outSuffix + "/Pdf_" + whichVarOut[var1] + "_x.dat";
              std::cout << "Opening file " << filename << std::endl;
              file = fopen(filename.c_str(), "w");
              for (int v1i = 0; v1i < nBins; v1i++)
                fprintf(file, "%e\n", vMin[var1] + dv1 * (0.5 + (Real)v1i));
              fclose(file);
              // Format: var2 range
              filename =
                infile + outSuffix + "/Pdf_" + whichVarOut[var2] + "_x.dat";
              std::cout << "Opening file " << filename << std::endl;
              file = fopen(filename.c_str(), "w");
              for (int v2i = 0; v2i < nBins; v2i++)
                fprintf(file, "%e\n", vMin[var2] + dv2 * (0.5 + (Real)v2i));
              fclose(file);
              // PdfX1
              filename = infile + outSuffix + "/PdfX1_" + whichVarOut[var1] +
                         "_" + whichVarOut[var2] + ".dat";
              std::cout << "Opening file " << filename << std::endl;
              file = fopen(filename.c_str(), "w");
              for (int v1i = 0; v1i < nBins; v1i++) {
                for (int v2i = 0; v2i < nBins; v2i++)
                  fprintf(file, "%e ", binX1[v1i * nBins + v2i]);
                fprintf(file, "\n");
              }
              fclose(file);
              // PdfX2
              filename = infile + outSuffix + "/PdfX2_" + whichVarOut[var1] +
                         "_" + whichVarOut[var2] + ".dat";
              std::cout << "Opening file " << filename << std::endl;
              file = fopen(filename.c_str(), "w");
              for (int v1i = 0; v1i < nBins; v1i++) {
                for (int v2i = 0; v2i < nBins; v2i++)
                  fprintf(file, "%e ", binX2[v1i * nBins + v2i]);
                fprintf(file, "\n");
              }
              fclose(file);
              for (int v = 0; v < nVars_condMean; v++) {
                // Format: condMean matrix
                filename = infile + outSuffix + "/condMean_" +
                           whichVarOut[ix_condMean[v]] + "_on_" +
                           whichVarOut[var1] + "_" + whichVarOut[var2] + ".dat";
                std::cout << "Opening file " << filename << std::endl;
                file = fopen(filename.c_str(), "w");
                for (int v1i = 0; v1i < nBins; v1i++) {
                  for (int v2i = 0; v2i < nBins; v2i++)
                    fprintf(file, "%e ", binCondMeans[v][v1i * nBins + v2i]);
                  fprintf(file, "\n");
                }
                fclose(file);
              }
            }

            if (output_tecplot) {
              // Output data for tecplot
              // Format: Header, each node, quadrilateral indices between nodes
              filename = infile + outSuffix + "/Pdf_" + whichVarOut[var1] +
                         "_" + whichVarOut[var2] + ".tpd";
              std::cout << "Opening file " << filename << std::endl;
              file = fopen(filename.c_str(), "w");
              fprintf(
                file, "VARIABLES = %s %s logpdf pdf\n", whichVar[var1].c_str(),
                whichVar[var2].c_str());
              fprintf(
                file, "ZONE N=%i E=%i F=FEPOINT ET=QUADRILATERAL\n",
                nBins * nBins, (nBins - 1) * (nBins - 1));
              for (int v1i = 0; v1i < nBins; v1i++) {
                Real v1 = vMin[var1] + dv1 * (0.5 + (Real)v1i);
                for (int v2i = 0; v2i < nBins; v2i++) {
                  Real v2 = vMin[var2] + dv2 * (0.5 + (Real)v2i);
                  Real p = bin[v1i * nBins + v2i];
                  fprintf(file, "%e %e %e %e\n", v1, v2, log(p + small), p);
                }
              }
              for (int v1i = 0; v1i < nBins - 1; v1i++) {
                for (int v2i = 0; v2i < nBins - 1; v2i++) {
                  int i1 = v1i * nBins + v2i + 1;
                  int i2 = (v1i + 1) * nBins + v2i + 1;
                  int i3 = (v1i + 1) * nBins + (v2i + 1) + 1;
                  int i4 = v1i * nBins + (v2i + 1) + 1;
                  fprintf(file, "%i %i %i %i\n", i1, i2, i3, i4);
                }
              }
              fclose(file);
            }

            if (output_fab) {
              // Output data in fab format
              filename = infile + outSuffix + "/Pdf_" + whichVarOut[var1] +
                         "_" + whichVarOut[var2] + ".fab";
              std::cout << "Opening file " << filename << std::endl;
              std::ofstream os(filename.c_str(), std::ios::binary);
              FArrayBox fab(
                Box(
                  IntVect::TheZeroVector(),
                  IntVect(AMREX_D_DECL(nBins - 1, nBins - 1, 0))),
                4);
              std::cout << "box: " << fab.box() << std::endl;
              for (int v1i = 0; v1i < nBins; v1i++) {
                Real v1 = vMin[var1] + dv1 * (0.5 + (Real)v1i);
                for (int v2i = 0; v2i < nBins; v2i++) {
                  Real v2 = vMin[var2] + dv2 * (0.5 + (Real)v2i);
                  Real p = bin[v1i * nBins + v2i];
                  IntVect iv(AMREX_D_DECL(v1i, v2i, 0));

                  fab(iv, 0) = v1;
                  fab(iv, 1) = v2;
                  fab(iv, 2) = log(p + small);
                  fab(iv, 3) = p;
                }
              }
              fab.writeOn(os);
              os.close();
            }

            if (output_scatter) {
              // Output non-zero bins as scatter plot
              // Format: x y
              filename = infile + outSuffix + "/Scatter_" + whichVarOut[var1] +
                         "_" + whichVarOut[var2] + ".dat";
              std::cout << "Opening file " << filename << std::endl;
              file = fopen(filename.c_str(), "w");
              for (int v1i = 0; v1i < nBins; v1i++) {
                Real v1 = vMin[var1] + dv1 * (0.5 + (Real)v1i);
                for (int v2i = 0; v2i < nBins; v2i++) {
                  Real v2 = vMin[var2] + dv2 * (0.5 + (Real)v2i);
                  if (bin[v1i * nBins + v2i] > 0)
                    fprintf(file, "%e %e\n", v1, v2);
                }
              }
              fclose(file);
            }

            iPair++;
          } // var2
        } // var1
      } // IOProcessor

      if (output_plotfile) {
        // Let's write the data as a "plotfile"
        //
        // Put the data into a multifab
        //
        Real small = 1.e-7;
        int ngrow(0);
        // Vector<int> pmap(1+1);
        // pmap[0] = ParallelDescriptor::IOProcessorNumber();
        // pmap[1] = ParallelDescriptor::MyProc();
        // DistributionMapping distMap(pmap);
        BoxArray ba(Box(
          IntVect::TheZeroVector(),
          IntVect(AMREX_D_DECL(nBins - 1, nBins - 1, 0))));
        DistributionMapping distMap(ba);
        MultiFab omf(ba, distMap, 2 * nPairs, ngrow);
        omf.setVal(0);
        for (MFIter ntmfi(omf); ntmfi.isValid(); ++ntmfi) {
          for (int iPair = 0; iPair < nPairs; iPair++) {
            Real* bin = binArray[iPair].dataPtr();
            Real* fab = omf[ntmfi].dataPtr(iPair);
            for (int v1i = 0; v1i < nBins; v1i++)
              for (int v2i = 0; v2i < nBins; v2i++)
                fab[v2i * nBins + v1i] = bin[v1i * nBins + v2i];
          }
          for (int iPair = 0; iPair < nPairs; iPair++) {
            Real* bin = binArray[iPair].dataPtr();
            Real* fab = omf[ntmfi].dataPtr(iPair + nPairs);
            for (int v1i = 0; v1i < nBins; v1i++)
              for (int v2i = 0; v2i < nBins; v2i++)
                fab[v2i * nBins + v1i] = log(small + bin[v1i * nBins + v2i]);
          }
        }
        //
        // Now do all the work for a plot file
        //
        std::string pltfile;
        if (outSuffix != "")
          pltfile = infile + outSuffix;
        else
          pltfile = infile + "jpdf";

        if (ParallelDescriptor::IOProcessor())
          if (!UtilCreateDirectory(pltfile, 0755))
            CreateDirectoryFailed(pltfile);
        ParallelDescriptor::Barrier();

        std::string HeaderFileName = pltfile + "/Header";

        static const std::string the_plot_file_type("NavierStokes-V1.1");

        if (ParallelDescriptor::IOProcessor()) {

          std::ofstream os;
          os.open(HeaderFileName.c_str());

          int old_prec = os.precision(15);

          // The plot file type
          os << the_plot_file_type << '\n';
          // The number of variables
          os << 2 * nPairs << '\n';
          // The variable names
          for (int var1 = 0; var1 < nVars_jpdf; var1++) {
            for (int var2 = var1 + 1; var2 < nVars_jpdf; var2++) {
              std::string variableName =
                "Pdf_" + whichVar[var1] + "_" + whichVar[var2];
              os << variableName << '\n';
            }
          }
          for (int var1 = 0; var1 < nVars_jpdf; var1++) {
            for (int var2 = var1 + 1; var2 < nVars_jpdf; var2++) {
              std::string variableName =
                "Pdf_" + whichVar[var1] + "_" + whichVar[var2] + " (log)";
              os << variableName << '\n';
            }
          }
          // The number of space dimensions
          os << "2" << '\n'; // Hardwire
          // Time
          os << amrData.Time() << '\n';
          // Finest level
          os << "0" << '\n'; // Hardwire
          // Domain
          os << "0 0\n";
          os << "1 1\n";
          // Refinement ratios
          os << '\n'; // Hardwire
          // Cell sizes
          os << "((0,0) (" << nBins - 1 << "," << nBins - 1 << ") (0,0))"
             << '\n';
          // Time steps
          // os << amrData.Step << '\n'; // FIXME
          os << "0" << '\n';
          // dx
          os << 1.0 / nBins << " " << 1.0 / nBins << '\n';
          // CoordSys & bndry
          os << "0\n0\n"; // Hardwire
          //
          // Now for the specific grids
          //
          // Level x- Number of grids - Time
          os << "0 1 " << amrData.Time() << '\n';
          // Time steps
          // os << amrData.Step << '\n'; // FIXME
          os << "0" << '\n';
          // Grid domain
          os << "0 1\n";
          os << "0 1\n";
          // Where to find the data
          os << "Level_0/Cell" << '\n';
          // Now let's add a bit so we know what the axes are
          for (int v = 0; v < nVars_jpdf; v++)
            os << vMin[v] << " " << vMax[v] << '\n';
          // That's all folks
          os.close();
        }
        // Build the directory to hold the MultiFab at this level.
        // The name is relative to the directory containing the Header file.
        //
        static const std::string BaseName = "/Cell";
        std::string Level = "Level_0";
        //
        // Now for the full pathname of that directory.
        //
        std::string FullPath = pltfile;
        if (!FullPath.empty() && FullPath[FullPath.length() - 1] != '/')
          FullPath += '/';
        FullPath += Level;
        //
        // Write the multifab
        //
        if (ParallelDescriptor::IOProcessor())
          if (!UtilCreateDirectory(FullPath, 0755))
            CreateDirectoryFailed(FullPath);
        ParallelDescriptor::Barrier();
        //
        // Use the Full pathname when naming the MultiFab.
        //
        std::string TheFullPath = FullPath;
        TheFullPath += BaseName;
        VisMF::Write(omf, TheFullPath);
      } // Output plotfile

    } // iPlot

    if (do_average) {
      // Reduce data
      for (int iPair = 0; iPair < nPairs; iPair++) {
        Real* binAv = binAvArray[iPair].dataPtr();
        Real* binAvX1 = binAvX1Array[iPair].dataPtr();
        Real* binAvX2 = binAvX2Array[iPair].dataPtr();
        ParallelDescriptor::ReduceRealSum(
          binAv, binAvArray[iPair].size(),
          ParallelDescriptor::IOProcessorNumber());
        ParallelDescriptor::ReduceRealSum(
          binAvX1, binAvX1Array[iPair].size(),
          ParallelDescriptor::IOProcessorNumber());
        ParallelDescriptor::ReduceRealSum(
          binAvX2, binAvX2Array[iPair].size(),
          ParallelDescriptor::IOProcessorNumber());

        for (int v = 0; v < nVars_condMean; v++) {
          Real* binAvCondMean = binAvCondMeanArray[iPair][v].dataPtr();
          ParallelDescriptor::ReduceRealSum(
            binAvCondMean, binAvCondMeanArray[iPair][v].size(),
            ParallelDescriptor::IOProcessorNumber());
        }
      }

      // Output the data to file
      if (ParallelDescriptor::IOProcessor()) {

        std::string infile = "JPDFAverage";
        std::string oFile = infile + outSuffix;
        if (!UtilCreateDirectory(oFile, 0755))
          CreateDirectoryFailed(oFile);

        Real small = 1.e-7;

        int iPair = 0;

        for (int var1 = 0; var1 < nVars_jpdf; var1++) {
          // Change in var1 between each bin
          Real dv1 = (vMax[var1] - vMin[var1]) / (Real)nBins;

          for (int var2 = var1 + 1; var2 < nVars_jpdf; var2++) {
            // Change in var2 between each bin
            Real dv2 = (vMax[var2] - vMin[var2]) / (Real)nBins;

            // Pull out the pointer for this pair
            Real* binAv = binAvArray[iPair].dataPtr();
            Real* binAvX1 = binAvX1Array[iPair].dataPtr();
            Real* binAvX2 = binAvX2Array[iPair].dataPtr();

            Vector<Real*> binAvCondMeans(nVars_condMean);
            for (int v = 0; v < nVars_condMean; v++) {
              binAvCondMeans[v] = binAvCondMeanArray[iPair][v].dataPtr();
            }

            // Divide Xi by bin vol to average (has to come first)
            for (int v1i = 0, i = 0; v1i < nBins; v1i++) {
              Real v1 = vMin[var1] + dv1 * (0.5 + (Real)v1i);
              for (int v2i = 0; v2i < nBins; v2i++, i++) {
                Real v2 = vMin[var2] + dv2 * (0.5 + (Real)v2i);
                Real div = binAv[i];
                if (div > 0) {
                  binAvX1[i] /= div;
                  binAvX2[i] /= div;
                  for (int v = 0; v < nVars_condMean; v++) {
                    binAvCondMeans[v][i] /= div;
                  }
                } else {
                  binAvX1[i] = v1;
                  binAvX2[i] = v2;
                }
              }
            }

            // Normalise so that sum(binAv) = 1. When conditioning is active,
            // divide by the actual accumulated volume across all files.
            Real normAvVol = domainVol * (Real)nPlotFiles;
            if (do_conditioning > 0) {
              Real totalAvVol = 0;
              for (int i = 0; i < nBins * nBins; i++)
                totalAvVol += binAv[i];
              if (totalAvVol > 0)
                normAvVol = totalAvVol;
            }
            for (int i = 0; i < nBins * nBins; i++)
              binAv[i] /= normAvVol;

            // Let's write some files...
            std::string filename;
            FILE* file;

            if (output_gnuplot) {
              // Output data for gnuplot
              // Format: x y pdf
              filename = infile + outSuffix + "/Pdf_" + whichVarOut[var1] +
                         "_" + whichVarOut[var2] + ".gpd";
              std::cout << "Opening file " << filename << std::endl;
              file = fopen(filename.c_str(), "w");
              for (int v1i = 0; v1i < nBins; v1i++) {
                Real v1 = vMin[var1] + dv1 * (0.5 + (Real)v1i);
                for (int v2i = 0; v2i < nBins; v2i++) {
                  Real v2 = vMin[var2] + dv2 * (0.5 + (Real)v2i);
                  fprintf(file, "%e %e %e\n", v1, v2, binAv[v1i * nBins + v2i]);
                }
              }
              fclose(file);
            }

            if (output_matlab) {
              // Output data for matlab (3 files)
              // Format: pdf matrix
              filename = infile + outSuffix + "/Pdf_" + whichVarOut[var1] +
                         "_" + whichVarOut[var2] + ".dat";
              std::cout << "Opening file " << filename << std::endl;
              file = fopen(filename.c_str(), "w");
              for (int v1i = 0; v1i < nBins; v1i++) {
                for (int v2i = 0; v2i < nBins; v2i++)
                  fprintf(file, "%e ", binAv[v1i * nBins + v2i]);
                fprintf(file, "\n");
              }
              fclose(file);
              // Format: var1 range
              filename =
                infile + outSuffix + "/Pdf_" + whichVarOut[var1] + "_x.dat";
              std::cout << "Opening file " << filename << std::endl;
              file = fopen(filename.c_str(), "w");
              for (int v1i = 0; v1i < nBins; v1i++)
                fprintf(file, "%e\n", vMin[var1] + dv1 * (0.5 + (Real)v1i));
              fclose(file);
              // Format: var2 range
              filename =
                infile + outSuffix + "/Pdf_" + whichVarOut[var2] + "_x.dat";
              std::cout << "Opening file " << filename << std::endl;
              file = fopen(filename.c_str(), "w");
              for (int v2i = 0; v2i < nBins; v2i++)
                fprintf(file, "%e\n", vMin[var2] + dv2 * (0.5 + (Real)v2i));
              fclose(file);
              // PdfX1
              filename = infile + outSuffix + "/PdfX1_" + whichVarOut[var1] +
                         "_" + whichVarOut[var2] + ".dat";
              std::cout << "Opening file " << filename << std::endl;
              file = fopen(filename.c_str(), "w");
              for (int v1i = 0; v1i < nBins; v1i++) {
                for (int v2i = 0; v2i < nBins; v2i++)
                  fprintf(file, "%e ", binAvX1[v1i * nBins + v2i]);
                fprintf(file, "\n");
              }
              fclose(file);
              // PdfX2
              filename = infile + outSuffix + "/PdfX2_" + whichVarOut[var1] +
                         "_" + whichVarOut[var2] + ".dat";
              std::cout << "Opening file " << filename << std::endl;
              file = fopen(filename.c_str(), "w");
              for (int v1i = 0; v1i < nBins; v1i++) {
                for (int v2i = 0; v2i < nBins; v2i++)
                  fprintf(file, "%e ", binAvX2[v1i * nBins + v2i]);
                fprintf(file, "\n");
              }
              fclose(file);
              for (int v = 0; v < nVars_condMean; v++) {
                // Format: condMean matrix
                filename = infile + outSuffix + "/condMean_" +
                           whichVarOut[ix_condMean[v]] + "_on_" +
                           whichVarOut[var1] + "_" + whichVarOut[var2] + ".dat";
                std::cout << "Opening file " << filename << std::endl;
                file = fopen(filename.c_str(), "w");
                for (int v1i = 0; v1i < nBins; v1i++) {
                  for (int v2i = 0; v2i < nBins; v2i++)
                    fprintf(file, "%e ", binAvCondMeans[v][v1i * nBins + v2i]);
                  fprintf(file, "\n");
                }
                fclose(file);
              }
            }

            if (output_tecplot) {
              // Output data for tecplot
              // Format: Header, each node, quadrilateral indices between nodes
              filename = infile + outSuffix + "/Pdf_" + whichVarOut[var1] +
                         "_" + whichVarOut[var2] + ".tpd";
              std::cout << "Opening file " << filename << std::endl;
              file = fopen(filename.c_str(), "w");
              fprintf(
                file, "VARIABLES = %s %s logpdf pdf\n", whichVar[var1].c_str(),
                whichVar[var2].c_str());
              fprintf(
                file, "ZONE N=%i E=%i F=FEPOINT ET=QUADRILATERAL\n",
                nBins * nBins, (nBins - 1) * (nBins - 1));
              for (int v1i = 0; v1i < nBins; v1i++) {
                Real v1 = vMin[var1] + dv1 * (0.5 + (Real)v1i);
                for (int v2i = 0; v2i < nBins; v2i++) {
                  Real v2 = vMin[var2] + dv2 * (0.5 + (Real)v2i);
                  Real p = binAv[v1i * nBins + v2i];
                  fprintf(file, "%e %e %e %e\n", v1, v2, log(p + small), p);
                }
              }
              for (int v1i = 0; v1i < nBins - 1; v1i++) {
                for (int v2i = 0; v2i < nBins - 1; v2i++) {
                  int i1 = v1i * nBins + v2i + 1;
                  int i2 = (v1i + 1) * nBins + v2i + 1;
                  int i3 = (v1i + 1) * nBins + (v2i + 1) + 1;
                  int i4 = v1i * nBins + (v2i + 1) + 1;
                  fprintf(file, "%i %i %i %i\n", i1, i2, i3, i4);
                }
              }
              fclose(file);
            }

            if (output_fab) {
              // Output data in fab format
              filename = infile + outSuffix + "/Pdf_" + whichVarOut[var1] +
                         "_" + whichVarOut[var2] + ".fab";
              std::cout << "Opening file " << filename << std::endl;
              std::ofstream os(filename.c_str());
              FArrayBox fab(
                Box(
                  IntVect::TheZeroVector(),
                  IntVect(AMREX_D_DECL(nBins - 1, nBins - 1, 0))),
                4);
              for (int v1i = 0; v1i < nBins; v1i++) {
                Real v1 = vMin[var1] + dv1 * (0.5 + (Real)v1i);
                for (int v2i = 0; v2i < nBins; v2i++) {
                  Real v2 = vMin[var2] + dv2 * (0.5 + (Real)v2i);
                  Real p = binAv[v1i * nBins + v2i];
                  IntVect iv(AMREX_D_DECL(v1i, v2i, 0));
                  fab(iv, 0) = v1;
                  fab(iv, 1) = v2;
                  fab(iv, 2) = log(p + small);
                  fab(iv, 3) = p;
                }
              }
              fab.writeOn(os);
              os.close();
            }

            if (output_scatter) {
              // Output non-zero bins as scatter plot
              // Format: x y
              filename = infile + outSuffix + "/Scatter_" + whichVarOut[var1] +
                         "_" + whichVarOut[var2] + ".dat";
              std::cout << "Opening file " << filename << std::endl;
              file = fopen(filename.c_str(), "w");
              for (int v1i = 0; v1i < nBins; v1i++) {
                Real v1 = vMin[var1] + dv1 * (0.5 + (Real)v1i);
                for (int v2i = 0; v2i < nBins; v2i++) {
                  Real v2 = vMin[var2] + dv2 * (0.5 + (Real)v2i);
                  if (binAv[v1i * nBins + v2i] > 0)
                    fprintf(file, "%e %e\n", v1, v2);
                }
              }
              fclose(file);
            }

            iPair++;
          } // var2
        } // var1
      } // IOProcessor
    }
  }
  Finalize();
  return 0;
}
