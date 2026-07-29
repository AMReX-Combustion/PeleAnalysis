#include <cstring>
#include <string>
#include <iostream>

#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include <AMReX_DataServices.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_PlotFileUtil.H>

using namespace amrex;

#include <optimalEstimatorANN.H>

static void
print_usage(int, char* argv[])
{
  std::cerr
    << "Evaluates a network trained by optimalEstimatorTraining at every "
       "cell of an\n"
       "AMReX plotfile and writes the target, the conditional estimate and "
       "the squared\n"
       "residual, whose volume average is the irreducible error.\n\n"

    << "Usage:\n"
    << "  " << argv[0]
    << " infile=FILE features=\"VAR1 ...\" targets=\"VAR1 ...\" "
       "neurons=\"N1 ...\" [OPTIONS]\n\n"

    << "Required arguments:\n"
    << "  infile=FILE             AMReX plotfile to evaluate the estimator "
       "on\n"
    << "  features=\"VAR1 ...\"     Conditioning variables; must match the "
       "training run\n"
    << "  targets=\"VAR1 ...\"      Estimated variables; must match the "
       "training run\n"
    << "  neurons=\"N1 ...\"        Hidden layers; must match the training "
       "run\n\n"

    << "Options:\n"
    << "  model_path=PATH         Network to load; \".pt\" is appended "
       "(DEF: optimal_estimator)\n"
    << "  minmax_path=PATH        Normalisation bounds; \".bin\" is appended "
       "(DEF: minmax)\n"
    << "  outfile=NAME            Output plotfile (DEF: <infile>_OE)\n"
    << "  finestLevel=N           Finest AMR level used (DEF: finest in "
       "file)\n"
    << "  is_per=I J [K]          Periodicity flags per direction (DEF: 1 1 "
       "1)\n"
    << "  num_threads=N           libtorch threads (DEF: 1 under MPI)\n"
    << "  -h, --help              Show this help message\n\n"

    << "The architecture is not stored in the checkpoint, so features, "
       "targets and\n"
    << "neurons must repeat the values used for training.\n"
    << "Visit PeleAnalysis/Src/InputSamples for examples or refer to "
       "the documentation.\n";

  std::exit(1);
}

std::string
getFileRoot(const std::string& infile)
{
  std::vector<std::string> tokens = Tokenize(infile, std::string("/"));
  while (!tokens.empty() && tokens.back().empty()) {
    tokens.pop_back();
  }
  return tokens.empty() ? infile : tokens.back();
}

int
main(int argc, char* argv[])
{
  Initialize(argc, argv);
  {
    if (argc < 2) {
      print_usage(argc, argv);
    } else if (
      (std::strcmp(argv[1], "-h") == 0) ||
      (std::strcmp(argv[1], "--help") == 0)) {
      print_usage(argc, argv);
    }

    ParmParse pp;

    if (pp.contains("help"))
      print_usage(argc, argv);

    // Open plotfile header and create an amrData object pointing into it
    std::string plotFileName;
    pp.get("infile", plotFileName);
    DataServices::SetBatchMode();
    Amrvis::FileType fileType(Amrvis::NEWPLT);
    DataServices dataServices(plotFileName, fileType);
    if (!dataServices.AmrDataOk()) {
      DataServices::Dispatch(DataServices::ExitRequest, NULL);
    }
    AmrData& amrData = dataServices.AmrDataRef();

    // Set up input field data names, and destination components to load data
    // upon read.
    int nFeatures = pp.countval("features");
    int nTargets = pp.countval("targets");
    int nCompIn = nFeatures + nTargets;
    Vector<std::string> inNames(nCompIn);
    Vector<int> destFillComps(nCompIn);
    Vector<std::string> features(nFeatures);
    pp.getarr("features", features);
    Vector<std::string> targets(nTargets);
    pp.getarr("targets", targets);
    for (int n = 0; n < nFeatures; n++) {
      inNames[n] = features[n];
      destFillComps[n] = n;
    }
    for (int n = 0; n < nTargets; n++) {
      inNames[n + nFeatures] = targets[n];
      destFillComps[n + nFeatures] = n + nFeatures;
    }

    // Loop over AMR levels in the plotfile, read the data and do work
    int finestLevel = amrData.FinestLevel();
    pp.query("finestLevel", finestLevel);
    finestLevel = std::min(finestLevel, amrData.FinestLevel());
    Vector<int> is_per(AMREX_SPACEDIM, 1);
    pp.queryarr("is_per", is_per, 0, AMREX_SPACEDIM);
    Print() << "Periodicity assumed for this case: ";
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
      Print() << is_per[idim] << " ";
    }
    Print() << "\n";

    std::string path = "optimal_estimator";
    pp.query("model_path", path);
    path += ".pt";
    int Nlev = finestLevel + 1;
    const int nGrow = 0;
    int nCompOut = 3 * nTargets;
    Vector<std::string> outNames(nCompOut);
    for (int n = 0; n < nTargets; n++) {
      outNames[n] = targets[n];
      outNames[n + nTargets] = targets[n] + "_cond_";
      for (int nf = 0; nf < nFeatures; nf++) {
        outNames[n + nTargets] += features[nf];
        if (nf != nFeatures - 1) {
          outNames[n + nTargets] += ",";
        }
      }
      outNames[n + 2 * nTargets] = "irr_" + outNames[n + nTargets];
    }

    // declare Vector of MultiFab, one MultiFab for each level
    Vector<MultiFab> outdata(Nlev);
    // also need to hold the geometry of each level
    Vector<Geometry> geoms(Nlev);

    RealBox rb(&(amrData.ProbLo()[0]), &(amrData.ProbHi()[0]));
    int n_layers = pp.countval("neurons");
    Vector<int> neurons(n_layers);
    pp.getarr("neurons", neurons);
    // set pytorch data type (default is float or torch::kFloat32)
    auto dtype0 = torch::kDouble;
    auto model = std::make_shared<Net>(nFeatures, neurons, nTargets);
    torch::load(model, path);

    auto tensoropt = torch::TensorOptions().dtype(dtype0);

    // Inference is a single forward pass per box; no autograd graph is needed.
    torch::NoGradGuard no_grad;
    model->eval();

    int num_threads = -1;
    pp.query("num_threads", num_threads);
    if (num_threads > 0) {
      torch::set_num_threads(num_threads);
    } else if (ParallelDescriptor::NProcs() > 1) {
      torch::set_num_threads(1);
    }
    Vector<Real> f_max(nFeatures, -1e100);
    Vector<Real> f_min(nFeatures, 1e100);
    Vector<Real> t_max(nTargets, -1e100);
    Vector<Real> t_min(nTargets, 1e100);

    std::string minmax_path = "minmax";
    pp.query("minmax_path", minmax_path);
    minmax_path += ".bin";

    std::ifstream file(minmax_path, std::ios::binary);
    if (!file.good()) {
      amrex::Abort("Could not open minmax file " + minmax_path);
    }
    file.read((char*)f_min.dataPtr(), sizeof(Real) * nFeatures);
    file.read((char*)f_max.dataPtr(), sizeof(Real) * nFeatures);
    file.read((char*)t_min.dataPtr(), sizeof(Real) * nTargets);
    file.read((char*)t_max.dataPtr(), sizeof(Real) * nTargets);
    if (!file) {
      amrex::Abort(
        "Short read from " + minmax_path +
        ": it does not match the requested number of features and targets");
    }
    file.close();

    for (int nv = 0; nv < nFeatures; nv++) {
      Print() << "f_max[" << nv << "] = " << f_max[nv] << std::endl;
      Print() << "f_min[" << nv << "] = " << f_min[nv] << std::endl;
    }
    for (int nv = 0; nv < nTargets; nv++) {
      Print() << "t_max[" << nv << "] = " << t_max[nv] << std::endl;
      Print() << "t_min[" << nv << "] = " << t_min[nv] << std::endl;
    }

    // now we have a trained model which should be the same on all processes

    Long nOutOfRange = 0;
    Long nSamples = 0;

    for (int lev = 0; lev < Nlev; ++lev) {
      // Get the array of boxes for this level
      const BoxArray ba = amrData.boxArray(lev);
      // Distribution mapping i.e. how are boxes distributed across processors
      const DistributionMapping dm(ba);
      // set up a multifab for the outdata with the same boxarray, distribution
      // mapping, components and ngrow
      outdata[lev] = MultiFab(ba, dm, nCompOut, nGrow);
      // we'll load each level as we go, so only need to fill one multifab at a
      // time
      MultiFab indata(ba, dm, nCompIn, nGrow);
      int coord = 0; // indicates cartesian
      geoms[lev] =
        Geometry(amrData.ProbDomain()[lev], &rb, coord, &(is_per[0]));

      Print() << "Reading data for level " << lev << std::endl;
      amrData.FillVar(indata, lev, inNames, destFillComps); // magic IO call
      Print() << "Data has been read for level " << lev << std::endl;

      MultiFab::Copy(
        outdata[lev], indata, nFeatures, 0, nTargets,
        nGrow); // copy the target data

      // Iterate over the multiFabs - distributes each box to a process and
      // iterates over it
      for (MFIter mfi(indata); mfi.isValid(); ++mfi) {
        // box for this iteration
        const Box& bx = mfi.validbox();
        const IntVect bx_lo = bx.smallEnd();

        const IntVect nbox = bx.size();
#if AMREX_SPACEDIM == 2
        int ncell = nbox[0] * nbox[1];
#else
        int ncell = nbox[0] * nbox[1] * nbox[2];
#endif
        // arrays for in and out data
        Array4<Real> const& inbox = indata.array(mfi);
        Array4<Real> const& outbox_t = outdata[lev].array(mfi, 0);
        Array4<Real> const& outbox_oe = outdata[lev].array(mfi, nTargets);
        Array4<Real> const& outbox_irr = outdata[lev].array(mfi, 2 * nTargets);
        // create array to feed to model
        // The checkpoint is always written in double precision, so the network
        // is evaluated in double regardless of how amrex::Real is configured.
        std::vector<double> trainingVec(ncell * nFeatures);
        double* trainingPtr = trainingVec.data();
        nSamples += (Long)ncell * nFeatures;
        amrex::LoopOnCpu(
          amrex::lbound(bx), amrex::ubound(bx), [&](int i, int j, int k) {
            int ii = i - bx_lo[0];
            int jj = j - bx_lo[1];
            int index = jj * nbox[0] + ii;
#if AMREX_SPACEDIM == 3
            int kk = k - bx_lo[2];
            index += kk * nbox[0] * nbox[1];
#endif
            for (int n = 0; n < nFeatures; n++) {
              const double xn = -1.0 + 2.0 * (inbox(i, j, k, n) - f_min[n]) /
                                         (f_max[n] - f_min[n]);
              // Features outside the training range put the network into
              // extrapolation, where the saturated tanh flattens the estimate.
              if (xn < -1.0 || xn > 1.0) {
                nOutOfRange++;
              }
              trainingPtr[index * nFeatures + n] = xn;
            }
          });
        torch::Tensor local_tensor =
          torch::from_blob(trainingPtr, {ncell, nFeatures}, tensoropt);
        torch::Tensor output = model->forward(local_tensor).contiguous();
        const double* outputPtr = output.data_ptr<double>();
        amrex::LoopOnCpu(
          amrex::lbound(bx), amrex::ubound(bx), [&](int i, int j, int k) {
            int ii = i - bx_lo[0];
            int jj = j - bx_lo[1];
            int index = jj * nbox[0] + ii;
#if AMREX_SPACEDIM == 3
            int kk = k - bx_lo[2];
            index += kk * nbox[0] * nbox[1];
#endif
            for (int n = 0; n < nTargets; n++) {
              outbox_oe(i, j, k, n) =
                t_min[n] + 0.5 * (t_max[n] - t_min[n]) *
                             (outputPtr[index * nTargets + n] + 1.0);
              outbox_irr(i, j, k, n) =
                (outbox_t(i, j, k, n) - outbox_oe(i, j, k, n)) *
                (outbox_t(i, j, k, n) - outbox_oe(i, j, k, n));
            }
          });
      }
      Print() << "Derive finished for level " << lev << std::endl;
    }

    ParallelDescriptor::ReduceLongSum(nOutOfRange);
    ParallelDescriptor::ReduceLongSum(nSamples);
    if (nOutOfRange > 0) {
      Print() << "\n*** WARNING: " << nOutOfRange << " of " << nSamples << " ("
              << 100.0 * (Real)nOutOfRange / (Real)std::max(Long(1), nSamples)
              << "%) feature values fall outside the range seen during "
                 "training.\n"
              << "    The network extrapolates there and the estimate is not "
                 "trustworthy.\n"
              << "    This is expected when inferring on a different snapshot "
                 "than the one\n"
              << "    used for training; retrain on a set that brackets the "
                 "inference data.\n\n";
    }

    std::string outfile = getFileRoot(plotFileName) + "_OE";
    pp.query("outfile", outfile);
    Print() << "Writing new data to " << outfile << std::endl;
    Vector<int> isteps(Nlev, 0);
    Vector<IntVect> refRatios(Nlev - 1);
    for (int lev = 0; lev < Nlev - 1; ++lev) {
      const int rr = amrData.RefRatio()[lev];
      refRatios[lev] = IntVect{AMREX_D_DECL(rr, rr, rr)};
    }
    amrex::WriteMultiLevelPlotfile(
      outfile, Nlev, GetVecOfConstPtrs(outdata), outNames, geoms, 0.0, isteps,
      refRatios);
  }
  Finalize();
  return 0;
}
