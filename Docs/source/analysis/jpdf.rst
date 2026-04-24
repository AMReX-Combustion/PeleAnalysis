.. highlight:: bash


jpdf
****

Given a plotfile and a list of variable names, compute joint probability density functions
(JPDFs) for every unique pair of variables, yielding :math:`n(n-1)/2` distributions.
Results are written to a directory named ``<infile><outSuffix>`` and can be output in
multiple formats (AMReX plotfile, gnuplot, MATLAB, Tecplot, FAB, or scatter).

The tool also supports computing **2D conditional means** of additional variables on the
JPDF axes, temporal averaging over a sequence of plotfiles, and optional conditioning to
restrict the data entering the PDF to a sub-range of a progress variable.

Usage: ::

   ./jpdf2d.gnu.MPI.ex infile=<plt> vars=<v1> <v2> [options]

Example: ::

   ./jpdf2d.gnu.MPI.ex ./InputSamples/jpdf.inp


Tool Options
############

::

   #------------------- IO CONTROL -----------------------------------------------------------
   infile = plt00000 plt00001                 # Plot file(s); multiple files can be listed for averaging
   outSuffix = _jpdf                          # DEF: "jpdf"; Suffix appended to infile path for output directory
   finestLevel = 2                            # DEF: finest level of plot file; Sets the finest level to read.

   #------------------- OUTPUT FORMATS -------------------------------------------------------
   output_plotfile = 1                        # [0, 1], DEF: 1; Write JPDFs as AMReX plotfile
   output_matlab   = 0                        # [0, 1], DEF: 0; Write JPDFs in MATLAB format (.dat)
   output_gnuplot  = 0                        # [0, 1], DEF: 0; Write JPDFs in gnuplot format (.gpd)
   output_tecplot  = 0                        # [0, 1], DEF: 0; Write JPDFs in Tecplot format (.tpd)
   output_fab      = 0                        # [0, 1], DEF: 0; Write JPDFs as FAB binary (.fab)
   output_scatter  = 0                        # [0, 1], DEF: 0; Write non-zero bins as scatter plot (.dat)

   #------------------- VARIABLES ------------------------------------------------------------
   vars = Y(H2O) temp                         # List of variables; JPDF computed for each unique pair (min. 2)
   nBins = 64                                 # DEF: 64; Number of bins per direction in JPDF

   #useminmax1 = 0.0 1.0                      # Override auto min/max for variable 1 (index starts at 1)
   #useminmax2 = 0.0 3000.0                   # Override auto min/max for variable 2

   #------------------- 2D CONDITIONAL MEAN -------------------------------------------------
   # condMean_vars computes the mean of each listed variable conditioned simultaneously on
   # both JPDF axes (v1 and v2). For each bin (i,j) of the v1-v2 plane, the tool computes
   #   <condMean_var | v1_i, v2_j>  =  sum( condMean_var * dV ) / sum( dV )
   # where the sums run over all cells whose v1 and v2 values fall in that bin.
   # This is a 2D conditional mean; use the conditionalMean tool for 1D conditioning.
   # Variables listed here may overlap with vars; duplicates are not loaded twice.
   condMean_vars = Y(OH)                      # Variable(s) to average conditioned on the 2D JPDF grid

   #------------------- AVERAGING ------------------------------------------------------------
   do_average = 1                             # [0, 1], DEF: 0; Accumulate JPDFs across all infiles and write
                                              # a single time-averaged result to JPDFAverage<outSuffix>

   #------------------- CONDITIONING --------------------------------------------------------
   # Restricts which cells contribute to the PDF. This is independent of condMean_vars.
   do_conditioning = 0                        # [0, 1, 2], DEF: 0; 0: no conditioning
                                              #   1: only cells where cMin <= cVar <= cMax contribute
                                              #   2: effective value is c(1-c); normalisation is forced
   #cVar = 0                                  # DEF: 0; Index in 'vars' of the conditioning variable
   #cMin = 0.0                                # DEF: 0.0; Min value of conditioning variable
   #cMax = 1.0                                # DEF: 1.0; Max value of conditioning variable
   #norm_cVal = 0                             # [0, 1], DEF: 0; Normalise cVar to [0,1] before cMin/cMax check
   #cNormMin = 0.0                            # DEF: 0.0; Min value for normalisation
   #cNormMax = 1.0                            # DEF: 1.0; Max value for normalisation


Details
#######

For each pair of variables ``(v1, v2)`` the data is binned onto a uniform ``nBins x nBins``
grid spanning ``[vMin, vMax]`` in each direction. The min and max are determined
automatically from the data across all AMR levels, but can be overridden with
``useminmax<i>``. Values that fall outside the bin range are clamped to the edge bin and
reported as out-of-range counts. The PDF is normalised so that its integral equals unity:

.. math::

   \int \int P(v_1, v_2) \, dv_1 \, dv_2 = 1

2D conditional mean
~~~~~~~~~~~~~~~~~~~

When ``condMean_vars`` is set, the tool computes, for each bin :math:`(i,j)` of the
:math:`v_1`–:math:`v_2` plane, the volume-weighted mean of each listed variable over all
cells that fall in that bin:

.. math::

   \langle \phi \mid v_{1,i},\, v_{2,j} \rangle
   = \frac{\displaystyle\sum_{\text{cells} \in (i,j)} \phi \, \Delta V}
          {\displaystyle\sum_{\text{cells} \in (i,j)} \Delta V}

This is a **2D conditional mean**: the conditioning is on two variables simultaneously,
which makes it a natural companion to the JPDF. It differs from the
``conditionalMean`` tool, which conditions on a single variable (1D binning).

Variables listed in ``condMean_vars`` may overlap with ``vars``; any duplicates are
detected automatically and the data is not loaded twice. Output files are named
``condMean_<condVar>_on_<v1>_<v2>.<ext>``.

Averaging over multiple plotfiles
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When ``do_average = 1`` and more than one plotfile is supplied in ``infile``, the
per-file PDFs and conditional means are accumulated. A single time-averaged result is
written to the directory ``JPDFAverage<outSuffix>`` once all files have been processed.
The averaged PDF is normalised by the total volume summed across all files.

Conditioning
~~~~~~~~~~~~

``do_conditioning`` restricts the cells that contribute to the PDF and to any conditional
means. It is independent of ``condMean_vars``.

* **1** – only cells where ``cMin <= cVar <= cMax`` contribute. If ``norm_cVal = 1``, the
  conditioning variable is first mapped to :math:`[0,1]` using ``cNormMin`` and
  ``cNormMax``.
* **2** – the effective conditioning value is :math:`c(1-c)`; normalisation is forced.
  ``cVar`` selects which entry in ``vars`` provides :math:`c`.

Output formats
~~~~~~~~~~~~~~

All formats write one file per variable pair per plotfile. In MATLAB mode an additional
file per conditional-mean variable is written. The AMReX plotfile format stores both the
PDF and its natural logarithm as separate components and encodes the bin axis ranges in
the plotfile header.

Testing
#######

A functional test suite for this tool lives in ``Tests/jpdf/``. It covers
all output formats, 2D conditional means, conditioning modes, temporal
averaging, and MPI correctness in both 2D and 3D. Run with::

   cd Tests/jpdf
   ./run_tests.sh

See :doc:`/testing/jpdf` for the full test matrix and expected outputs.
