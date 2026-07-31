#include <optimalEstimatorANN.H> //contains and libraries for prescribed architecture

#include <cstring>
#include <algorithm>
#include <type_traits>
#include <vector>

using namespace amrex;

// amrex::Real as a torch dtype, so staged host data can be wrapped directly.
static constexpr torch::ScalarType amrexDtype =
  std::is_same_v<Real, double> ? torch::kFloat64 : torch::kFloat32;

static void
print_usage(int, char* argv[])
{
  std::cerr
    << "Trains a neural network to approximate the conditional mean\n"
       "<targets | features>, the optimal estimator, on an AMReX plotfile.\n"
       "Use optimalEstimatorInfer to evaluate the trained network.\n\n"

    << "Usage:\n"
    << "  " << argv[0]
    << " infile=FILE features=\"VAR1 ...\" targets=\"VAR1 ...\" "
       "neurons=\"N1 ...\" [OPTIONS]\n\n"

    << "Required arguments:\n"
    << "  infile=FILE             AMReX plotfile holding the training data\n"
    << "  features=\"VAR1 ...\"     Conditioning variables (network inputs)\n"
    << "  targets=\"VAR1 ...\"      Variables whose conditional mean is "
       "sought\n"
    << "  neurons=\"N1 ...\"        Neurons per hidden layer, one entry per "
       "layer\n\n"

    << "Options:\n"
    << "  model_path=PATH         Where to write the network; \".pt\" is "
       "appended (DEF: optimal_estimator)\n"
    << "  minmax_path=PATH        Where to write the normalisation bounds; "
       "\".bin\" is appended (DEF: minmax)\n"
    << "  minLevel=N              Coarsest AMR level used; cells covered by a "
       "finer level are skipped (DEF: 0)\n"
    << "  finestLevel=N           Finest AMR level used (DEF: finest in "
       "file)\n"
    << "  split=F                 Fraction of boxes used for training "
       "(DEF: 0.7)\n"
    << "  nEpochs=N               Upper bound on training epochs (DEF: "
       "1000)\n"
    << "  batch_size=N            Samples per mini-batch (DEF: 16384)\n"
    << "  learning_rate=F         Initial Adam step size (DEF: 1e-3)\n"
    << "  alpha=F                 L2 weight decay; biases the estimate "
       "(DEF: 0)\n"
    << "  volume_weight=0|1       Weight each sample by its cell volume "
       "(DEF: 1)\n"
    << "  use_double=0|1          Train in double precision (DEF: 0)\n"
    << "  num_threads=N           libtorch threads (DEF: 1 under MPI)\n\n"

    << "Convergence options:\n"
    << "  minEpochs=N             Earliest epoch that may early-stop (DEF: "
       "100)\n"
    << "  patience=N              Stop after N epochs without improvement; "
       "0 disables (DEF: 50)\n"
    << "  min_delta=F             Smallest improvement in R2 that counts "
       "(DEF: 1e-3)\n"
    << "  lr_patience=N           Reduce the rate after N flat epochs; 0 "
       "disables (DEF: 20)\n"
    << "  lr_factor=F             Learning-rate multiplier on a plateau "
       "(DEF: 0.5)\n"
    << "  min_lr=F                Learning-rate floor (DEF: 1e-6)\n"
    << "  print_every=N           Epoch summary interval (DEF: 1)\n"
    << "  -h, --help              Show this help message\n\n"

    << "The weights of the best-scoring epoch are written, always in double "
       "precision.\n"
    << "Visit PeleAnalysis/Src/InputSamples for examples or refer to "
       "the documentation.\n";

  std::exit(1);
}

// Sum the gradients of every parameter across all ranks and rescale to the
// mean, so that each rank applies the same update as a serial run over the
// union of the data. Must be called after backward() and before step().
static void
allReduceGradients(const std::shared_ptr<Net>& model)
{
  const int nProcs = ParallelDescriptor::NProcs();
  if (nProcs == 1) {
    return;
  }

  for (const auto& p : model->parameters()) {
    const auto& g = p.grad();
    if (!g.defined()) {
      continue;
    }
    // The reduction writes through the raw buffer, so it must be contiguous.
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
      g.is_contiguous(), "Non-contiguous gradient cannot be reduced in place");

    if (g.scalar_type() == torch::kFloat64) {
      ParallelDescriptor::ReduceRealSum(g.data_ptr<double>(), g.numel());
    } else {
      ParallelDescriptor::ReduceRealSum(g.data_ptr<float>(), g.numel());
    }
    g.div_(static_cast<double>(nProcs));
  }
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

    const int nProcs = ParallelDescriptor::NProcs();

    // libtorch otherwise sizes its thread pool from the visible core count,
    // which oversubscribes badly once several ranks share a node.
    int num_threads = -1;
    pp.query("num_threads", num_threads);
    if (num_threads > 0) {
      torch::set_num_threads(num_threads);
    } else if (nProcs > 1) {
      torch::set_num_threads(1);
    }

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
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
      nFeatures > 0 && nTargets > 0,
      "At least one feature and one target must be given");

    int nEpochs = 1000;
    pp.query("nEpochs", nEpochs);
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
    int minLevel = 0;
    pp.query("finestLevel", finestLevel);
    finestLevel = std::min(finestLevel, amrData.FinestLevel());
    pp.query("minLevel", minLevel);
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
      minLevel >= 0 && minLevel <= finestLevel,
      "minLevel must lie between 0 and finestLevel");

    std::string path = "optimal_estimator";
    pp.query("model_path", path);
    path += ".pt";
    int Nlev = finestLevel + 1;
    const int nGrow = 0;
    int n_layers = pp.countval("neurons");
    Vector<int> neurons(n_layers);
    pp.getarr("neurons", neurons);
    torch::manual_seed(42); // identical model init on every rank

    // Training is done in single precision by default: the optimal estimator
    // is a statistical quantity and does not need the extra digits, while
    // float32 roughly doubles the throughput of the matrix products. The
    // checkpoint is always written in double precision (see below), so
    // optimalEstimatorInfer is unaffected by this choice.
    int use_double = 0;
    pp.query("use_double", use_double);
    const torch::ScalarType trainDtype =
      use_double ? torch::kFloat64 : torch::kFloat32;

    Vector<Real> f_max(nFeatures, std::numeric_limits<Real>::lowest());
    Vector<Real> f_min(nFeatures, std::numeric_limits<Real>::max());
    Vector<Real> t_max(nTargets, std::numeric_limits<Real>::lowest());
    Vector<Real> t_min(nTargets, std::numeric_limits<Real>::max());

    Vector<MultiFab> indata(Nlev);

    // A coarse cell that lies under a finer level is not an independent
    // sample: it covers the same physical region as the fine cells that
    // replace it, and its value is whatever the plotfile happens to carry
    // there. Taking every cell of every level between minLevel and finestLevel
    // would therefore feed each refined region to the network once per level
    // covering it, weighting it by the depth of refinement rather than by its
    // volume. makeFineMask flags exactly those covered cells, so only the
    // uncovered ones - the same set AMReX itself uses for a volume average -
    // are kept below.
    Vector<iMultiFab> validMask(Nlev);

    for (int lev = minLevel; lev < Nlev; lev++) {
      // Get the array of boxes for this level
      const BoxArray ba = amrData.boxArray(lev);
      // Distribution mapping i.e. how are boxes distributed across processors
      const DistributionMapping dm(ba);

      indata[lev].define(ba, dm, nCompIn, nGrow);
      Print() << "Reading data for level " << lev << std::endl;
      amrData.FillVar(indata[lev], lev, inNames, destFillComps); // magic IO
                                                                 // call
      Print() << "Data has been read for level " << lev << std::endl;

      if (lev < finestLevel) {
        // 1 where the cell is uncovered and usable, 0 where a finer level
        // takes over. The finer level is the one requested here, not
        // necessarily the finest in the file: with finestLevel=N the level-N
        // data is used whole, exactly as if the file stopped there.
        const int rr = amrData.RefRatio()[lev];
        validMask[lev] = makeFineMask(
          ba, dm, amrData.boxArray(lev + 1), IntVect(AMREX_D_DECL(rr, rr, rr)),
          /*crse_value=*/1, /*fine_value=*/0);
      } else {
        validMask[lev].define(ba, dm, 1, nGrow);
        validMask[lev].setVal(1);
      }
    }

    // Normalisation bounds over the uncovered cells only, so that they
    // describe exactly the data the network is trained on. The number of
    // uncovered cells per local box is collected in the same sweep, in the
    // order the staging loop below walks them.
    Long nValidLocal = 0;
    Long nCoveredLocal = 0;
    Vector<Long> boxCells;
    for (int lev = minLevel; lev < Nlev; lev++) {
      for (MFIter mfi(indata[lev]); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        Array4<const Real> const& fab = indata[lev].const_array(mfi);
        Array4<const int> const& msk = validMask[lev].const_array(mfi);

        Long nValidBox = 0;
        amrex::LoopOnCpu(
          amrex::lbound(bx), amrex::ubound(bx), [&](int i, int j, int k) {
            if (msk(i, j, k) == 0) {
              nCoveredLocal++;
              return;
            }
            nValidBox++;
            for (int nv = 0; nv < nFeatures; nv++) {
              f_max[nv] = std::max(f_max[nv], fab(i, j, k, nv));
              f_min[nv] = std::min(f_min[nv], fab(i, j, k, nv));
            }
            for (int nv = 0; nv < nTargets; nv++) {
              const Real v = fab(i, j, k, nv + nFeatures);
              t_max[nv] = std::max(t_max[nv], v);
              t_min[nv] = std::min(t_min[nv], v);
            }
          });
        boxCells.push_back(nValidBox);
        nValidLocal += nValidBox;
      }
    }
    // The loop above is per rank; make the bounds global.
    ParallelDescriptor::ReduceRealMax(f_max.dataPtr(), nFeatures);
    ParallelDescriptor::ReduceRealMin(f_min.dataPtr(), nFeatures);
    ParallelDescriptor::ReduceRealMax(t_max.dataPtr(), nTargets);
    ParallelDescriptor::ReduceRealMin(t_min.dataPtr(), nTargets);

    {
      Long nValidGlobal = nValidLocal;
      Long nCoveredGlobal = nCoveredLocal;
      ParallelDescriptor::ReduceLongSum(nValidGlobal);
      ParallelDescriptor::ReduceLongSum(nCoveredGlobal);
      AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        nValidGlobal > 0,
        "Every cell between minLevel and finestLevel is covered by a finer "
        "level; there is nothing to train on.");
      Print() << "Cells on levels " << minLevel << "-" << finestLevel << ": "
              << nValidGlobal + nCoveredGlobal << ", of which " << nCoveredGlobal
              << " are covered by a finer level and are skipped, leaving "
              << nValidGlobal << " samples." << std::endl;
    }

    for (int nv = 0; nv < nFeatures; nv++) {
      Print() << "f_min[" << nv << "] = " << f_min[nv] << ", f_max[" << nv
              << "] = " << f_max[nv] << std::endl;
      AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        f_max[nv] > f_min[nv],
        "Feature is constant over the data set; normalisation would divide by "
        "zero. Remove it from the feature list.");
    }
    for (int nv = 0; nv < nTargets; nv++) {
      Print() << "t_min[" << nv << "] = " << t_min[nv] << ", t_max[" << nv
              << "] = " << t_max[nv] << std::endl;
      AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        t_max[nv] > t_min[nv],
        "Target is constant over the data set; normalisation would divide by "
        "zero.");
    }

    // Split the local boxes into a training and a validation set. Splitting on
    // whole boxes rather than individual cells keeps spatially adjacent - and
    // therefore strongly correlated - cells off both sides of the split.
    Real split = 0.7;
    pp.query("split", split);
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
      split > 0.0 && split < 1.0, "split must lie strictly between 0 and 1");

    const int nBoxesLocal = boxCells.size();

    // A box that a finer level covers completely carries no samples and must
    // stay out of the split, or it could be handed the entire validation set
    // and leave nothing to validate on.
    Vector<int> boxOrder;
    for (int i = 0; i < nBoxesLocal; i++) {
      if (boxCells[i] > 0) {
        boxOrder.push_back(i);
      }
    }
    const int nUsableBoxes = boxOrder.size();
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
      nUsableBoxes >= 2,
      "This rank owns fewer than two boxes holding uncovered cells and cannot "
      "be split into a training and a validation set. Run with fewer MPI "
      "ranks.");

    // Deterministic, but decorrelated between ranks.
    std::mt19937 g(12345 + ParallelDescriptor::MyProc());
    std::shuffle(boxOrder.begin(), boxOrder.end(), g);

    const int nTrainBoxes =
      std::max(1, std::min(nUsableBoxes - 1, (int)(split * nUsableBoxes)));
    std::vector<char> isTrainBox(nBoxesLocal, 0);
    for (int i = 0; i < nTrainBoxes; i++) {
      isTrainBox[boxOrder[i]] = 1;
    }

    // A conditional mean is an average over volume, but a cell is one sample
    // whatever its size, so on a multi-level set the fit would be pulled
    // towards the refined regions: they contribute r^DIM samples where the
    // unrefined ones contribute a single, physically much larger, cell. Giving
    // every sample the volume of its cell as a weight in the loss restores the
    // volume average, and the network then converges on the conditional mean
    // rather than on a cell-count-weighted approximation to it. The weights are
    // expressed relative to a finest-level cell; on a single level they are all
    // 1 and the weighted loss reduces exactly to the unweighted one.
    int volume_weight = 1;
    pp.query("volume_weight", volume_weight);
    Vector<Real> levWeight(Nlev, 1.0);
    if (volume_weight) {
      Real vFinest = 1.0;
      for (int idim = 0; idim < AMREX_SPACEDIM; idim++) {
        vFinest *= amrData.DxLevel()[finestLevel][idim];
      }
      for (int lev = minLevel; lev < Nlev; lev++) {
        Real v = 1.0;
        for (int idim = 0; idim < AMREX_SPACEDIM; idim++) {
          v *= amrData.DxLevel()[lev][idim];
        }
        levWeight[lev] = v / vFinest;
        if (lev < finestLevel) {
          Print() << "Level " << lev << " cells carry weight " << levWeight[lev]
                  << " in the loss (volume relative to a level-" << finestLevel
                  << " cell)" << std::endl;
        }
      }
    } else if (minLevel < finestLevel) {
      Print() << "\n*** WARNING: volume_weight=0 with minLevel < finestLevel.\n"
              << "    Every cell counts once regardless of its size, so the "
                 "fit is weighted by\n"
              << "    cell count instead of by volume and is not the "
                 "conditional mean.\n\n";
    }

    // Stage all cells into contiguous host buffers, normalised onto [-1,1].
    // Building one large block rather than one small tensor per box is what
    // makes the mini-batches big enough to be worth dispatching.
    std::vector<Real> trainFeat, trainTarg, trainW, valFeat, valTarg, valW;
    trainFeat.reserve(nValidLocal * nFeatures);
    trainTarg.reserve(nValidLocal * nTargets);
    trainW.reserve(nValidLocal);

    int boxCount = 0;
    for (int lev = minLevel; lev < Nlev; lev++) {
      const Real w = levWeight[lev];
      for (MFIter mfi(indata[lev]); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        Array4<Real> const& fab = indata[lev].array(mfi);
        Array4<const int> const& msk = validMask[lev].const_array(mfi);

        const bool train = isTrainBox[boxCount++] != 0;
        std::vector<Real>& fOut = train ? trainFeat : valFeat;
        std::vector<Real>& tOut = train ? trainTarg : valTarg;
        std::vector<Real>& wOut = train ? trainW : valW;

        amrex::LoopOnCpu(
          amrex::lbound(bx), amrex::ubound(bx), [&](int i, int j, int k) {
            if (msk(i, j, k) == 0) {
              return; // covered by a finer level, already sampled there
            }
            wOut.push_back(w);
            for (int n = 0; n < nFeatures; n++) {
              fOut.push_back(
                -1.0 +
                2.0 * (fab(i, j, k, n) - f_min[n]) / (f_max[n] - f_min[n]));
            }
            for (int n = 0; n < nTargets; n++) {
              tOut.push_back(
                -1.0 + 2.0 * (fab(i, j, k, n + nFeatures) - t_min[n]) /
                         (t_max[n] - t_min[n]));
            }
          });
      }
    }

    // The plotfile data is no longer needed; release it before allocating the
    // tensors so the two copies never coexist.
    indata.clear();
    validMask.clear();

    const int64_t nTrainLocal = (int64_t)trainTarg.size() / nTargets;
    const int64_t nValLocal = (int64_t)valTarg.size() / nTargets;

    auto toTensor = [&](std::vector<Real>& host, int64_t nrows, int64_t ncols) {
      auto t = torch::from_blob(
        host.data(), {nrows, ncols}, torch::TensorOptions().dtype(amrexDtype));
      // Force a copy so the tensor owns its storage even when trainDtype ==
      // amrexDtype, where .to() returns the from_blob tensor unchanged and
      // `out` would alias the buffer freed on the next line.
      auto out =
        t.to(trainDtype, /*non_blocking=*/false, /*copy=*/true).contiguous();
      std::vector<Real>().swap(host);
      return out;
    };

    torch::Tensor train_x = toTensor(trainFeat, nTrainLocal, nFeatures);
    torch::Tensor train_y = toTensor(trainTarg, nTrainLocal, nTargets);
    // One column, so that it broadcasts over the target columns below.
    torch::Tensor train_w = toTensor(trainW, nTrainLocal, 1);
    torch::Tensor val_x = toTensor(valFeat, nValLocal, nFeatures);
    torch::Tensor val_y = toTensor(valTarg, nValLocal, nTargets);
    torch::Tensor val_w = toTensor(valW, nValLocal, 1);

    // batch_size now counts SAMPLES, not the box edge length it used to mean.
    int64_t batch_size = 16384;
    pp.query("batch_size", batch_size);
    if (batch_size <= 0) {
      batch_size = nTrainLocal;
    }
    if (batch_size < 256) {
      Print() << "\n*** WARNING: batch_size = " << batch_size
              << " is very small.\n"
              << "    batch_size now counts SAMPLES (cells). It used to be a "
                 "box edge\n"
              << "    length, so an old input with batch_size=32 meant "
              << (AMREX_SPACEDIM == 2 ? 1024 : 32768) << " samples in "
              << AMREX_SPACEDIM << "D.\n"
              << "    Small batches are dominated by framework overhead; "
                 "prefer >= 4096.\n\n";
    }

    // Every rank must take the same number of optimiser steps or the gradient
    // reductions will not line up. The per-epoch reshuffle means the tail
    // dropped on data-rich ranks is a different one each epoch.
    int nSteps = (int)((nTrainLocal + batch_size - 1) / batch_size);
    ParallelDescriptor::ReduceIntMin(nSteps);
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
      nSteps > 0, "Some rank has no training data. Run with fewer MPI ranks.");

    {
      Long nTrainGlobal = nTrainLocal;
      Long nValGlobal = nValLocal;
      ParallelDescriptor::ReduceLongSum(nTrainGlobal);
      ParallelDescriptor::ReduceLongSum(nValGlobal);
      Print() << "Training samples: " << nTrainGlobal
              << ", validation samples: " << nValGlobal << "\n"
              << "Batch size: " << batch_size << ", steps per epoch: " << nSteps
              << ", ranks: " << nProcs
              << ", torch threads: " << torch::get_num_threads() << std::endl;
    }

    // Variance of the targets in normalised space, used to turn the raw MSE
    // into a scale-free coefficient of determination. A model that only ever
    // predicts the unconditional mean scores R2 = 0. Weighted like the loss it
    // normalises, so that R2 keeps that meaning on a multi-level set.
    Real targVar = 1.0;
    {
      Real s = (val_w * val_y).sum().item<double>() +
               (train_w * train_y).sum().item<double>();
      Real s2 = (val_w * val_y.pow(2)).sum().item<double>() +
                (train_w * train_y.pow(2)).sum().item<double>();
      Real w = (Real)nTargets * (val_w.sum().item<double>() +
                                 train_w.sum().item<double>());
      ParallelDescriptor::ReduceRealSum(s);
      ParallelDescriptor::ReduceRealSum(s2);
      ParallelDescriptor::ReduceRealSum(w);
      const Real mean = s / w;
      targVar = std::max(s2 / w - mean * mean, Real(1e-30));
    }

    Real learning_rate = 1e-3;
    pp.query("learning_rate", learning_rate);

    // L2 penalty. Defaults to zero: regularisation shrinks the estimate and so
    // inflates the irreducible error, which is the quantity being measured.
    Real alpha = 0.0;
    pp.query("alpha", alpha);

    auto model = std::make_shared<Net>(nFeatures, neurons, nTargets);
    model->to(trainDtype);

    // Adam's built-in weight decay adds the penalty gradient directly, instead
    // of building an autograd graph over every parameter on every step.
    auto optimizer = torch::optim::Adam(
      model->parameters(),
      torch::optim::AdamOptions(learning_rate).weight_decay(alpha));

    // Convergence control
    int minEpochs = 100;
    int patience = 50;
    Real min_delta = 1e-3; // in R2 units
    int lr_patience = 20;
    Real lr_factor = 0.5;
    Real min_lr = 1e-6;
    int print_every = 1;
    pp.query("minEpochs", minEpochs);
    pp.query("patience", patience);
    pp.query("min_delta", min_delta);
    pp.query("lr_patience", lr_patience);
    pp.query("lr_factor", lr_factor);
    pp.query("min_lr", min_lr);
    pp.query("print_every", print_every);

    const Real min_delta_abs = min_delta * targVar;

    std::vector<torch::Tensor> best_state;
    auto snapshot = [&]() {
      torch::NoGradGuard ng;
      best_state.clear();
      for (const auto& p : model->parameters()) {
        best_state.push_back(p.detach().clone());
      }
    };
    auto restore = [&]() {
      torch::NoGradGuard ng;
      auto ps = model->parameters();
      for (size_t i = 0; i < ps.size(); i++) {
        ps[i].copy_(best_state[i]);
      }
    };

    auto setLearningRate = [&](Real lr) {
      for (auto& group : optimizer.param_groups()) {
        static_cast<torch::optim::AdamOptions&>(group.options()).lr(lr);
      }
    };

    // Volume-weighted mean squared error over the whole (distributed) set.
    // Accumulating the weighted sum of squares and dividing by the global sum
    // of weights keeps the result independent of how the data happens to be
    // spread over the ranks.
    auto evaluate = [&](const torch::Tensor& x, const torch::Tensor& y,
                        const torch::Tensor& w) {
      torch::NoGradGuard ng;
      model->eval();
      const int64_t n = x.size(0);
      Real sse = 0.0;
      Real wsum = 0.0;
      for (int64_t s = 0; s < n; s += batch_size) {
        const int64_t m = std::min(batch_size, n - s);
        auto wb = w.narrow(0, s, m);
        auto out = model->forward(x.narrow(0, s, m));
        sse += (wb * (out - y.narrow(0, s, m)).pow(2)).sum().item<double>();
        wsum += (Real)nTargets * wb.sum().item<double>();
      }
      ParallelDescriptor::ReduceRealSum(sse);
      ParallelDescriptor::ReduceRealSum(wsum);
      return wsum > 0.0 ? sse / wsum : Real(0);
    };

    Real best_val = std::numeric_limits<Real>::max();
    int best_epoch = -1;
    int bad_epochs = 0;
    int lr_bad_epochs = 0;
    Real lr = learning_rate;
    int epoch = 0;

    for (epoch = 0; epoch < nEpochs; ++epoch) {
      model->train();

      auto perm = torch::randperm(
        nTrainLocal, torch::TensorOptions().dtype(torch::kLong));

      Real epoch_training_loss = 0.0;
      for (int s = 0; s < nSteps; ++s) {
        const int64_t off = (int64_t)s * batch_size;
        const int64_t m = std::min(batch_size, nTrainLocal - off);
        auto idx = perm.narrow(0, off, m);

        auto xb = train_x.index_select(0, idx);
        auto yb = train_y.index_select(0, idx);
        auto wb = train_w.index_select(0, idx);

        optimizer.zero_grad();
        auto out = model->forward(xb);
        // Weighted MSE: sum(w * e^2) / (sum(w) * nTargets). With uniform
        // weights this is identical to torch::mse_loss(out, yb).
        auto loss = (wb * (out - yb).pow(2)).sum() / (wb.sum() * nTargets);
        loss.backward();
        allReduceGradients(model);
        optimizer.step();

        epoch_training_loss += loss.item<double>();
      }
      epoch_training_loss /= std::max(1, nSteps);
      ParallelDescriptor::ReduceRealSum(epoch_training_loss);
      epoch_training_loss /= (Real)nProcs;

      const Real val_loss = evaluate(val_x, val_y, val_w);
      const Real r2 = 1.0 - val_loss / targVar;

      // Every rank sees the same reduced val_loss, so they all take the same
      // branch here and stay in lockstep.
      if (val_loss < best_val - min_delta_abs) {
        best_val = val_loss;
        best_epoch = epoch;
        bad_epochs = 0;
        lr_bad_epochs = 0;
        snapshot();
      } else {
        bad_epochs++;
        lr_bad_epochs++;
        if (best_epoch < 0) {
          // Keep something to fall back on even if we never improve.
          best_val = val_loss;
          best_epoch = epoch;
          snapshot();
        }
      }

      if (
        print_every > 0 && (epoch % print_every == 0 || epoch == nEpochs - 1)) {
        Print() << "Epoch [" << epoch + 1 << "/" << nEpochs
                << "], Training Loss: " << epoch_training_loss
                << ", Validation Loss: " << val_loss << ", R2: " << r2
                << ", lr: " << lr << std::endl;
      }

      if (lr_patience > 0 && lr_bad_epochs >= lr_patience && lr > min_lr) {
        lr = std::max(lr * lr_factor, min_lr);
        setLearningRate(lr);
        lr_bad_epochs = 0;
        Print() << "  validation loss plateaued, reducing learning rate to "
                << lr << std::endl;
      }

      if (epoch + 1 >= minEpochs && patience > 0 && bad_epochs >= patience) {
        Print() << "Early stopping at epoch " << epoch + 1 << ": no "
                << "improvement of more than " << min_delta << " in R2 for "
                << patience << " epochs." << std::endl;
        break;
      }
    }

    if (best_epoch >= 0) {
      restore();
      Print() << "Restored weights from epoch " << best_epoch + 1
              << " (validation loss " << best_val << ", R2 "
              << 1.0 - best_val / targVar << ")." << std::endl;
    }

    if (ParallelDescriptor::IOProcessor()) {
      // Always checkpoint in double precision so that optimalEstimatorInfer,
      // which runs in double, can load the file regardless of trainDtype.
      model->to(torch::kFloat64);
      torch::save(model, path);

      std::string minmax_path = "minmax";
      pp.query("minmax_path", minmax_path);
      minmax_path += ".bin";
      writeMinMax(
        minmax_path, features, targets, neurons, f_min, f_max, t_min, t_max);
      Print() << "Wrote " << path << " and " << minmax_path << std::endl;
    }
  }
  Finalize();
  return 0;
}
