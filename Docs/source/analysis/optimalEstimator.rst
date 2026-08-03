optimalEstimator
================

Description
-----------

The *optimal estimator* quantifies how much of a quantity :math:`\phi` can, in
principle, be predicted from a chosen set of variables :math:`\mathbf{q}`. It is
defined as the conditional expectation

.. math::

   \langle \phi \mid \mathbf{q} \rangle = \mathrm{E}\!\left[\phi \mid \mathbf{q}\right],

which is the estimator that minimises the mean square error among *all*
functions of :math:`\mathbf{q}`. The residual that no such function can remove,

.. math::

   \varepsilon_{\mathrm{irr}}^2
     = \left\langle \left( \phi - \langle \phi \mid \mathbf{q} \rangle \right)^2 \right\rangle,

is the **irreducible error**. It is a property of the data set and the chosen
feature set alone — not of any particular model — and therefore provides a lower
bound on the error of any closure written in terms of :math:`\mathbf{q}`. A large
irreducible error means the feature set is missing information; a small one means
a well-chosen model form could in principle succeed.

Estimating :math:`\mathrm{E}[\phi \mid \mathbf{q}]` by binning becomes infeasible
beyond two or three features. Following Berger et al. (2018), these tools instead
approximate the conditional mean with a small feed-forward neural network trained
under a mean-square-error loss, whose minimiser is exactly the conditional mean.

The average in :math:`\langle \cdot \rangle` is over volume, so on an AMR data
set each sample enters the loss weighted by the volume of its cell:

.. math::

   \mathcal{L} = \frac{\sum_c V_c \left( \phi_c - f(\mathbf{q}_c) \right)^2}
                      {\sum_c V_c} ,

with :math:`V_c` the cell volume and the sum running over the cells that are not
covered by a finer level. Minimising this converges on the conditional mean of
the *field*; an unweighted sum would converge on a cell-count average instead,
in which a refined region counts :math:`r^{\mathrm{DIM}}` times as heavily as
the same volume left unrefined. On a single level the two coincide.

The functionality is split across two executables:

``optimalEstimatorTraining``
   Reads a plotfile, extracts the feature and target fields, trains the network
   and writes the trained weights plus the normalisation constants to disk.

``optimalEstimatorInfer``
   Reads a plotfile and a previously trained network, evaluates the estimator at
   every cell and writes a new plotfile containing the target, the conditional
   estimate and the pointwise squared residual.

Training and inference may be run on the same plotfile (in-sample estimate of the
irreducible error) or on different snapshots of the same configuration
(out-of-sample estimate). Both executables must be given the *same* ``features``,
``targets`` and ``neurons`` arguments — the network architecture is not stored
inside the checkpoint and is rebuilt from these values before the weights are
loaded. All three are recorded in ``<minmax_path>.bin`` and checked when it is
read back (see *Normalisation*), so a mismatch aborts rather than producing a
wrong answer.


Building
--------

Both tools link against the PyTorch C++ API (libtorch) and are therefore not
built by the default CI matrix. In ``Src/GNUmakefile`` they are listed behind a
double ``##`` comment; select one explicitly on the command line:

.. code-block:: bash

   make -j EBASE=optimalEstimatorTraining DIM=2
   make -j EBASE=optimalEstimatorInfer    DIM=2

The build first checks whether ``python3 -c "import torch"`` succeeds and, if so,
reuses that installation's headers and libraries. On the JUPITER-HPC system (status 07/2026) use the following module:

.. code-block:: bash

   module load GCC/14.3.0 OpenMPI/5.0.8 PyTorch/2.9.1

If no importable ``torch`` is found, the build falls back to a standalone
libtorch under ``Tools/libtorch``, which can be fetched once with

.. code-block:: bash

   make pytorch EBASE=optimalEstimatorTraining

.. note::

   libtorch is distributed in both the pre-C++11 and the C++11 ABI. The build
   sets ``-D_GLIBCXX_USE_CXX11_ABI=1``, so the C++11-ABI package must be used.
   Linker errors mentioning ``std::__cxx11::basic_string`` almost always mean the
   wrong package was picked up.


Usage
-----

.. code-block:: bash

   optimalEstimatorTraining infile=<s> features=<s s ...> targets=<s s ...> neurons=<i i ...> [options]
   optimalEstimatorInfer    infile=<s> features=<s s ...> targets=<s s ...> neurons=<i i ...> [options]


Input File
----------

Training:

.. code-block:: none

   infile        = plt00000
   features      = temp density
   targets       = I_R(prog)
   neurons       = 6 12 8
   nEpochs       = 1000
   minLevel      = 0
   batch_size    = 32
   learning_rate = 1e-3
   alpha         = 0.01
   model_path    = optimal_estimator
   minmax_path   = minmax

Inference:

.. code-block:: none

   infile      = plt00000
   features    = temp density
   targets     = I_R(prog)
   neurons     = 6 12 8
   is_per      = 1 0
   model_path  = optimal_estimator
   minmax_path = minmax
   outfile     = plt00000_OE


Shared Parameters
~~~~~~~~~~~~~~~~~

``infile``
   Path to the input AMReX plotfile. It must contain every field named in
   ``features`` and ``targets``.

``features``
   Names of the conditioning variables :math:`\mathbf{q}`, e.g.
   ``progVar Z``. The network input layer has one node per feature.

``targets``
   Names of the fields :math:`\phi` whose conditional mean is sought, e.g.
   ``I_R(progVar)``. Several targets may be estimated at once; they share one
   network with a multi-output layer.

``neurons``
   Number of neurons in each hidden layer, given as a space-separated list.
   ``neurons = 32 64 16`` builds a network with three hidden layers. Hidden
   layers use a ``tanh`` activation and the output layer is linear, following
   Berger et al. (2018). Weights are initialised with Xavier uniform sampling
   and a fixed seed, so a given configuration trains reproducibly.

``model_path``
   Path where the trained network is written (training) or read from
   (inference). The extension ``.pt`` is appended automatically.
   Default: ``optimal_estimator``.

``minmax_path``
   Path of the binary file holding the feature and target normalisation bounds.
   The extension ``.bin`` is appended automatically. Default: ``minmax``.

``finestLevel``
   Finest AMR level to include. Default: plotfile finest level.


Training Parameters
~~~~~~~~~~~~~~~~~~~

``nEpochs``
   Number of training epochs. Default: ``1000``.

``minLevel``
   Coarsest AMR level to include in the training set. Default: ``0``.

   Cells of a level that are covered by the next finer level are excluded
   automatically, so no physical region enters the sample more than once
   regardless of how ``minLevel`` and ``finestLevel`` are chosen. The mask comes
   from ``amrex::makeFineMask``, so the retained cells are exactly the ones
   AMReX itself uses for a volume average. The number of skipped cells is
   reported at start-up::

      Cells on levels 0-2: 1310720, of which 786432 are covered by a finer
      level and are skipped, leaving 524288 samples.

   What remains is weighted by cell volume (see ``volume_weight``), so the two
   grid-dependent biases of a multi-level training set — a region sampled once
   per level covering it, and a refined region outweighing an unrefined one of
   the same size — are both removed.

``batch_size``
   Number of **samples** (cells) per mini-batch. Default: ``16384``.

   .. warning::

      This parameter changed meaning. It used to be a box edge length passed to
      ``BoxArray::maxSize``, so an old input file with ``batch_size = 32``
      requested batches of :math:`32^{\mathrm{DIM}}` cells — 1024 in 2D, 32768
      in 3D. Reusing such an input verbatim now asks for 32-sample batches,
      which is several hundred times slower. The tool prints a warning when
      ``batch_size < 256``.

   Very small batches spend nearly all of their time in framework overhead
   rather than arithmetic. Values in the range :math:`10^4`–:math:`10^5` give
   the best throughput. Set ``-1`` to use the whole local training set as a
   single batch (full-batch gradient descent).

``split``
   Fraction of boxes used for training; the remainder is held out for
   validation and early stopping. Default: ``0.7``.

   The partition is drawn over the whole data set with a fixed seed, keyed on
   each box's level and its index into that level's ``BoxArray`` — quantities
   that come from the plotfile header and so do not depend on the distribution
   mapping. The same boxes therefore end up in the validation set whatever the
   rank count, and the fraction is exact globally rather than applied to each
   rank's own box count. The counts are reported at start-up::

      Split 76 boxes into 53 for training and 23 for validation.

   Boxes wholly covered by a finer level hold no samples and take no part in the
   split.

``volume_weight``
   Set to ``0`` to give every cell the same weight in the loss instead of the
   volume of its cell. Default: ``1``.

   The weights are volumes relative to a cell of ``finestLevel``, so a level
   whose cells are :math:`r` times as wide carries weight
   :math:`r^{\mathrm{DIM}}` and ``finestLevel`` itself carries 1. They are
   reported at start-up::

      Level 0 cells carry weight 16 in the loss (volume relative to a level-2 cell)
      Level 1 cells carry weight 4 in the loss (volume relative to a level-2 cell)

   The same weights are applied to the validation loss and to the target
   variance behind :math:`R^2`, so both keep their meaning. On a single level
   (``minLevel = finestLevel``) every weight is 1 and the setting has no effect
   whatsoever. Turning it off on a multi-level set makes the fit a cell-count
   average rather than the conditional mean and is only useful for reproducing
   an older run; the tool prints a warning when you do.

``learning_rate``
   Initial step size of the Adam optimiser. Default: ``1e-3``.

``alpha``
   Strength of the :math:`L_2` weight decay applied by Adam. Default: ``0``.

   .. warning::

      Regularisation shrinks the network output towards zero and therefore
      *biases the estimated conditional mean*, which inflates the reported
      irreducible error — the very quantity being measured. Since the optimal
      estimator is a property of the data rather than a predictive model,
      overfitting is far less of a concern here than bias, which is why the
      default is zero. Early stopping on the validation loss already guards
      against fitting noise. If you do enable it, verify that the result is
      insensitive to the value.

``use_double``
   Set to ``1`` to train in double precision. Default: ``0`` (single
   precision), which roughly doubles throughput and is entirely adequate for a
   statistical quantity. This affects training only — the checkpoint is always
   written in double precision, so ``optimalEstimatorInfer`` is unaffected.

``num_threads``
   Number of threads libtorch may use for the matrix products. Default: the
   libtorch default when running on one rank, and ``1`` under MPI, where
   letting every rank size its own pool from the visible core count would
   oversubscribe the node badly.


Convergence Control
~~~~~~~~~~~~~~~~~~~

Training stops when the validation loss stops improving, and the weights from
the best epoch — not the last — are the ones written to disk.

Progress is measured with the coefficient of determination

.. math::

   R^2 = 1 - \frac{\mathrm{MSE}_\mathrm{val}}{\mathrm{Var}(\phi)},

which is dimensionless and bounded above by 1, so a single threshold works
across every target and every case. :math:`R^2 = 0` corresponds to an estimator
no better than the unconditional mean. Both the validation MSE and the variance
are volume-weighted, so the ratio is unaffected by the weighting.

``nEpochs``
   Hard upper bound on the number of epochs. Default: ``1000``.

``minEpochs``
   Minimum number of epochs before early stopping may trigger, so that Adam's
   initial transient is not mistaken for a plateau. Default: ``100``.

``patience``
   Stop after this many consecutive epochs without an improvement larger than
   ``min_delta``. Set to ``0`` to disable early stopping and always run
   ``nEpochs``. Default: ``50``.

``min_delta``
   Smallest improvement in :math:`R^2` that counts as progress.
   Default: ``1e-3``.

``lr_patience``, ``lr_factor``, ``min_lr``
   After ``lr_patience`` epochs without improvement the learning rate is
   multiplied by ``lr_factor``, down to a floor of ``min_lr``. Reducing the
   step size on a plateau usually reaches a lower final loss in fewer total
   epochs. Set ``lr_patience = 0`` to keep the rate fixed.
   Defaults: ``20``, ``0.5``, ``1e-6``.

``print_every``
   Print the epoch summary every N epochs. Default: ``1``.


Inference Parameters
~~~~~~~~~~~~~~~~~~~~

``is_per``
   Periodicity flags in each spatial direction (0: non-periodic, 1: periodic).
   Used only to build the output plotfile geometry. Default: ``1 1 1``.

``outfile``
   Name of the output plotfile. Default: the input file's basename with ``_OE``
   appended.


Normalisation
-------------

Both features and targets are linearly mapped onto :math:`[-1, 1]` using the
minimum and maximum found over the cells actually used for training — every
included AMR level of the *training* file, minus the cells covered by a finer
level:

.. math::

   \tilde{x} = -1 + 2\,\frac{x - x_{\min}}{x_{\max} - x_{\min}} .

These bounds are written to ``<minmax_path>.bin`` and read back during inference
so that the estimator is de-normalised consistently. The file is
self-describing: a header records the format version, the feature and target
names and the hidden-layer widths, followed by the bounds themselves:

.. code-block:: none

   char[8]  magic "PeleOEMM"
   int32    format version
   int32    number of features, int32 number of targets
   per name, features then targets: int32 length, then that many chars
   int32    number of hidden layers, then that many int32 widths
   double   f_min[nF], f_max[nF], t_min[nT], t_max[nT]

``optimalEstimatorInfer`` checks all of it against its own ``features``,
``targets`` and ``neurons`` arguments — names and layer widths, in order — and
aborts on any disagreement::

   amrex::Abort::0::minmax.bin was trained with features 'progVar Z', but
   'Z progVar' were requested. Training and inference must be given the same
   list, in the same order.

This is the only place the architecture is checked at all: ``torch::load`` does
*not* object to a checkpoint whose layers differ from the network it is loaded
into, so before this check a mistyped ``neurons`` ran to completion and wrote a
plotfile full of nonsense with a zero exit code.

The payload is always double, like the network checkpoint, so the two
executables may be built with different ``PRECISION`` without difficulty.

Files written before the header was introduced are still read — they are a bare
dump of the four arrays in the writing build's own precision — but nothing in
them can be verified, so a warning is printed recommending a retrain. Both
formats assume the reader and writer share an endianness, as AMReX plotfiles
themselves largely do.

.. note::

   When inference is run on a different snapshot than training, features may fall
   outside the training range and are then mapped outside :math:`[-1, 1]`. The
   network extrapolates there, and because the hidden activations saturate the
   estimate flattens rather than diverging. ``optimalEstimatorInfer`` counts
   these samples and prints a warning with the affected fraction; a percentage
   large enough to matter means the training set does not bracket the inference
   data and the estimator should be retrained on a set that does.


Training Procedure
------------------

The uncovered cells of every included level are gathered box by box, normalised,
and the resulting mini-batches are shuffled once and split 70 % / 30 % into a
training and a validation set. The split is performed on whole boxes rather than
individual cells, which keeps spatially adjacent — and therefore strongly
correlated — cells from appearing on both sides of the split. Boxes that a finer
level covers completely hold no samples and take no part in the split.

Each epoch reshuffles the training samples, loops over the mini-batches
performing an Adam update per batch, then evaluates the validation loss and
prints a summary:

.. code-block:: none

   Epoch [1/1000], Training Loss: 0.0421, Validation Loss: 0.0438, R2: 0.79, lr: 0.001

Both losses are volume-weighted mean-square errors in the *normalised* target
space, so they are dimensionless and comparable across targets. Both are
computed as a weighted sum of squared errors divided by the global sum of
weights, which makes them independent of how the data happens to be distributed
over the MPI ranks, and makes the two numbers on the line above the same
statistic.

Training ends either at ``nEpochs`` or when the validation loss plateaus (see
*Convergence Control* above). The weights are then rolled back to the best
epoch, and the network is written to ``<model_path>.pt`` with the normalisation
bounds in ``<minmax_path>.bin``.

Under MPI each rank normalises its loss by the weight of the mini-batch summed
over *all* ranks, so its gradient is that rank's share of the global gradient
and the shares are simply summed after every backward pass. A run on N ranks
therefore takes the same optimiser steps as a serial run over the union of the
data, whatever the spread in how much data — or how much *weight* — each rank
holds.

.. note::

   Rescaling the summed gradients by :math:`1/N` instead, which is what
   averaging per-rank mean gradients amounts to, is only equivalent when every
   rank's mini-batch carries the same weight. On a multi-level set it does not:
   the distribution mapping is built independently per level, so ranks hold
   uncorrelated shares of coarse and fine boxes. The difference is not a
   rescaling that the optimiser absorbs but a different objective — it
   up-weights the ranks holding the least total weight, and in the limit of one
   level per rank it cancels the volume weighting outright, returning the fit to
   the cell-count average that ``volume_weight`` exists to avoid.

All ranks take the same number of steps per epoch, which is what keeps them
issuing the same sequence of collectives; because the training set is reshuffled
every epoch, the tail that this drops on data-rich ranks is a different one each
time.


Output
------

``optimalEstimatorInfer`` writes a plotfile with ``3 * nTargets`` components. For
each target :math:`\phi` with feature list :math:`q_1,\dots,q_n`:

- ``<target>`` — the target field, copied unchanged from the input
- ``<target>_cond_<q1>-...-<qn>`` — the optimal estimator
  :math:`\langle \phi \mid \mathbf{q}\rangle`
- ``irr_<target>_cond_<q1>-...-<qn>`` — the pointwise squared residual
  :math:`(\phi - \langle \phi \mid \mathbf{q}\rangle)^2`

The irreducible error itself is the *volume average* of the third field and is
obtained by post-processing the output, for example with ``integral``:

.. code-block:: bash

   integral infile=plt00000_OE vars='irr_I_R(progVar)_cond_progVar-Z' \
            integralDimension=2 avg=1

Normalising this by the variance of the target gives a dimensionless measure
between 0 (the features fully determine the target) and 1 (they carry no
information at all).

The output inherits the domain geometry and box structure of the input plotfile.


Typical Workflow
----------------

A complete study compares the irreducible error of competing feature sets for the
same target, so that the feature set retaining the most information can be
identified. Assuming ``plt_train`` and ``plt_test`` are two snapshots that already
contain the progress variables and mixture fraction (see ``progVar`` and
``computeMixtureFraction``):

.. code-block:: bash

   # 1. Train the estimator on the training snapshot, finest level only
   optimalEstimatorTraining infile=plt_train \
       features='progVar' 'Z' targets='I_R(progVar)' \
       neurons=32 64 16 nEpochs=1000 minLevel=2 batch_size=16384 \
       learning_rate=1e-3 \
       model_path=progVar-Z/optimal_estimator \
       minmax_path=progVar-Z/minmax 2>&1 | tee progVar-Z/training.log

   # 2. Evaluate it on the test snapshot
   optimalEstimatorInfer infile=plt_test \
       features='progVar' 'Z' targets='I_R(progVar)' \
       neurons=32 64 16 is_per=1 0 \
       model_path=progVar-Z/optimal_estimator \
       minmax_path=progVar-Z/minmax

   # 3. Volume-average the squared residual to get the irreducible error
   integral infile=plt_test_OE \
       vars='irr_I_R(progVar)_cond_progVar-Z' integralDimension=2 avg=1

   # 4. Optionally, look at the residual conditioned on the progress variable
   combinePlts infiles=plt_test plt_test_OE outfile=plt_test_results \
       vars='progVar' 'irr_I_R(progVar)_cond_progVar-Z' is_per=1 0
   conditionalMean infile=plt_test_results doBin=true binComp=0 avgComps=1 \
       binMin=0 binMax=1 nBins=256 outSuffix='_conditionalMean'

Repeating steps 1–3 for each candidate feature set and comparing the resulting
averages is the core of the analysis. Because the network is only an approximation
to the conditional mean, the absolute value of a single irreducible error should
be treated with care; *differences* between feature sets computed with an
identical architecture and training budget are far more robust.


Notes and Current Limitations
-----------------------------

- **CPU only.** Both tools run on the host. ``USE_CUDA`` is not supported.

- **Every rank needs at least one training box.** The split is global, so a rank
  may legitimately draw no *validation* boxes — harmless, since it then simply
  contributes nothing to either global sum — but a rank left with no *training*
  data has no batch to take its optimiser step from and the tool aborts. With a
  small number of large boxes this limits how many ranks are useful; splitting
  the plotfile into more boxes (or using fewer ranks) resolves it.

- **Throughput.** Training cost is ``nEpochs`` times the number of mini-batches.
  Make the job request more than one core — libtorch parallelises the matrix
  products over threads (``--cpus-per-task`` on Slurm) — and keep ``batch_size``
  large enough that the batches are worth dispatching. When sweeping many
  feature combinations, running the independent trainings concurrently (a Slurm
  job array, say) scales far better than adding ranks to a single one.

- **Mixed resolutions.** Covered cells are masked out and the rest are weighted
  by volume, so a multi-level training set is no longer biased by the grid in
  either of those two ways. What remains is that a coarse cell carries a
  *filtered* value while a fine cell carries a resolved one, and no weighting
  makes the conditional mean of the one equal to the conditional mean of the
  other. Training on a single level (``minLevel = finestLevel``) sidesteps it.

- **Degenerate fields.** A feature or target that is constant over the whole data
  set cannot be normalised; the tool aborts with a message naming the problem.


References
----------

A. Moreau, O. Teytaud, J.-P. Bertoglio, "Optimal estimation for large-eddy
simulation of turbulence and application to the analysis of subgrid models",
Physics of Fluids **18**, 105101 (2006).

B. Berger, Lukas, Konstantin Kleinheinz, Antonio Attili, Fabrizio Bisetti, 
Heinz Pitsch, and Michael E. Mueller. "Numerically accurate computational 
techniques for optimal estimator analyses of multi-parameter models." 
Combustion Theory and Modelling 22, no. 3 (2018): 480-504.
