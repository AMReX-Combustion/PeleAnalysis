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
loaded.


Building
--------

Both tools link against the PyTorch C++ API (libtorch) and are therefore not
built by the default CI matrix. In ``Src/GNUmakefile`` they are listed behind a
double ``##`` comment; select one explicitly on the command line:

.. code-block:: bash

   make -j EBASE=optimalEstimatorTraining DIM=2
   make -j EBASE=optimalEstimatorInfer    DIM=2

The build first checks whether ``python3 -c "import torch"`` succeeds and, if so,
reuses that installation's headers and libraries. On an HPC system this is
usually the right path:

.. code-block:: bash

   module load GCC OpenMPI PyTorch

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

   .. warning::

      With the default ``minLevel = 0`` the training set contains coarse cells
      that lie underneath finer grids, so the same physical region enters the
      sample more than once and the coarse (filtered) values bias the
      conditional mean. On refined data sets, set ``minLevel = finestLevel`` so
      that a single, uniformly resolved level is used.

``batch_size``
   Maximum box edge length, in cells, used to split the plotfile boxes
   (``BoxArray::maxSize``). Each resulting box becomes one mini-batch, so the
   actual number of samples per batch is ``batch_size**AMREX_SPACEDIM``, not
   ``batch_size``. A value of ``32`` in 2D therefore gives batches of 1024 cells.
   Default: ``-1`` (keep the plotfile's own box decomposition).

``learning_rate``
   Step size of the Adam optimiser. Default: ``1e-3``.

``alpha``
   Weight of the :math:`L_2` weight-decay term in the loss,

   .. math::

      \mathcal{L} = (1-\alpha)\,\mathrm{MSE} + \alpha \sum_i w_i^2 .

   Default: ``0.01``.

   .. warning::

      Regularisation shrinks the network output towards zero and therefore
      *biases the estimated conditional mean*, which inflates the reported
      irreducible error — the very quantity being measured. Since the optimal
      estimator is a property of the data rather than a predictive model,
      overfitting is far less of a concern here than bias. Use the smallest
      value that still trains stably, and verify that the result is insensitive
      to it. Note also that the penalty is an unnormalised sum over all weights,
      so its magnitude grows with network size and a value tuned for one
      architecture does not transfer to another.


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
minimum and maximum found over all included AMR levels of the *training* file:

.. math::

   \tilde{x} = -1 + 2\,\frac{x - x_{\min}}{x_{\max} - x_{\min}} .

These bounds are written to ``<minmax_path>.bin`` as a flat binary record of
``Real`` values in the order ``f_min, f_max, t_min, t_max``, and are read back
during inference so that the estimator is de-normalised consistently.

Because the file carries no header, the number of features and targets and the
floating-point precision must match between training and inference. Building one
executable with ``PRECISION = DOUBLE`` and the other with ``PRECISION = FLOAT``
will silently produce nonsense.

.. note::

   When inference is run on a different snapshot than training, features may fall
   outside the training range and are then mapped outside :math:`[-1, 1]`. The
   network extrapolates there, and because the hidden activations saturate the
   estimate flattens rather than diverging. Check the printed bounds against the
   inference data set before interpreting the result.


Training Procedure
------------------

The cells of every included level are gathered box by box, normalised, and the
resulting mini-batches are shuffled once and split 70 % / 30 % into a training
and a validation set. The split is performed on whole boxes rather than
individual cells, which keeps spatially adjacent — and therefore strongly
correlated — cells from appearing on both sides of the split.

Each epoch loops over the training batches performing an Adam update per batch,
then evaluates the validation loss, and prints both:

.. code-block:: none

   Epoch [1/1000], Training Loss: 0.0421, Validation Loss: 0.0438

Both losses are mean-square errors in the *normalised* target space, so they are
dimensionless and comparable across targets. A useful reference point is that a
model which predicts the unconditional mean of a target scores an MSE equal to
that target's normalised variance; a converged optimal estimator must do
substantially better than this, otherwise the chosen features carry essentially
no information about the target.

After the final epoch the network is written to ``<model_path>.pt`` and the
normalisation bounds to ``<minmax_path>.bin``.


Output
------

``optimalEstimatorInfer`` writes a plotfile with ``3 * nTargets`` components. For
each target :math:`\phi` with feature list :math:`q_1,\dots,q_n`:

- ``<target>`` — the target field, copied unchanged from the input
- ``<target>_cond_<q1>,...,<qn>`` — the optimal estimator
  :math:`\langle \phi \mid \mathbf{q}\rangle`
- ``irr_<target>_cond_<q1>,...,<qn>`` — the pointwise squared residual
  :math:`(\phi - \langle \phi \mid \mathbf{q}\rangle)^2`

The irreducible error itself is the *volume average* of the third field and is
obtained by post-processing the output, for example with ``integral``:

.. code-block:: bash

   integral infile=plt00000_OE vars='irr_I_R(progVar)_cond_progVar,Z' \
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
       neurons=32 64 16 nEpochs=1000 minLevel=2 batch_size=32 \
       learning_rate=1e-3 alpha=1e-6 \
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
       vars='irr_I_R(progVar)_cond_progVar,Z' integralDimension=2 avg=1

   # 4. Optionally, look at the residual conditioned on the progress variable
   combinePlts infiles=plt_test plt_test_OE outfile=plt_test_results \
       vars='progVar' 'irr_I_R(progVar)_cond_progVar,Z' is_per=1 0
   conditionalMean infile=plt_test_results doBin=true binComp=0 avgComps=1 \
       binMin=0 binMax=1 nBins=256 outSuffix='_conditionalMean'

Repeating steps 1–3 for each candidate feature set and comparing the resulting
averages is the core of the analysis. Because the network is only an approximation
to the conditional mean, the absolute value of a single irreducible error should
be treated with care; *differences* between feature sets computed with an
identical architecture and training budget are far more robust.


Notes and Current Limitations
-----------------------------

- **Run in serial.** Although the tools build with ``USE_MPI=TRUE``, the
  gradients computed on the individual ranks are not synchronised before the
  optimiser update, so every rank trains its own model on its own subset of the
  data and only the I/O rank's model is written out. Until this is addressed,
  ``optimalEstimatorTraining`` must be run with a single rank. Inference is
  embarrassingly parallel and unaffected.

- **CPU only.** The CUDA code paths are present but untested and do not currently
  build; use ``USE_CUDA=FALSE``.

- **Performance.** Training cost is dominated by the number of optimiser steps,
  which is ``nEpochs`` times the number of boxes. A small ``batch_size`` produces
  very many tiny batches whose cost is dominated by framework overhead rather
  than arithmetic. Prefer larger ``batch_size`` values (giving batches of
  :math:`10^4`–:math:`10^5` cells) and make sure the job requests more than one
  core, since libtorch parallelises the matrix products over threads
  (``--cpus-per-task`` on Slurm, or ``OMP_NUM_THREADS``).

- **Fixed epoch count.** Training always runs the full ``nEpochs``; there is no
  early stopping, and the best-scoring weights are not retained. Inspect the
  printed validation loss to confirm that it has plateaued and is not rising
  again before trusting the result.

- **Degenerate fields.** A feature or target that is constant over the whole data
  set makes ``x_max - x_min`` vanish and the normalisation produces ``NaN``.


References
----------

A. Moreau, O. Teytaud, J.-P. Bertoglio, "Optimal estimation for large-eddy
simulation of turbulence and application to the analysis of subgrid models",
Physics of Fluids **18**, 105101 (2006).

The network architecture (``tanh`` hidden layers with a linear output layer)
follows "Berger et al. (2018)" as cited in the source; the full reference still
needs to be filled in here.
