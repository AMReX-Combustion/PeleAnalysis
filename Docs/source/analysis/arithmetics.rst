arithmetics
===========

Description
-----------
``arithmetics`` applies a simple arithmetic operation (``+``, ``-``, ``*``, or ``/``)
to two field variables from an AMReX plotfile. The operator is selected at runtime
via the ``operator`` parameter. The tool reads a plotfile, computes the result of
``A OP B``, and writes a new plotfile containing all original variables plus the
derived field.

Usage
-----
.. code-block:: bash

   arithmetics3d.gnu.MPI.OMP.ex infile=FILE inVarAName=NAME inVarBName=NAME outVarName=NAME operator=add [OPTIONS]

Parameters
~~~~~~~~~~

``infile``
   Path to the AMReX plotfile.

``inVarAName``
   Name of operand A (must exist in the plotfile).

``inVarBName``
   Name of operand B (must exist in the plotfile).

``outVarName``
   Name of the resulting output variable.

``operator``
   Arithmetic operator to apply: ``add``, ``subtract``, ``multiply``, or ``divide``.

``outfile``
   Output plotfile name.
   Default: ``<infile>_<operator>``.

``finestLevel``
   Finest AMR level to process.
   Default: plotfile finest level.

``is_per``
   Periodicity flags in each spatial direction (0: non-periodic, 1: periodic).
   Default: ``1 1 1``.

``verbose``
   Enable verbose AMReX output if present.

``checkDivByZero``
   When ``operator=divide``, check whether operand B contains any zero values
   before performing the division.  Set to ``0`` to skip the check.
   Default: ``1`` (check enabled).

``coord``
   Coordinate system identifier passed to the output ``Geometry``.
   ``0`` = Cartesian (default).  Other values follow the AMReX convention.

``n_files``
   Maximum number of binary files used to write the output plotfile data
   (AMReX ``VisMF::SetNOutFiles``). Lower this to reduce the number of files
   created for large parallel post-processing runs. AMReX clamps the value to
   the number of MPI ranks, so a serial run always writes a single data file.
   Default: the AMReX default.

Output
------
A new AMReX plotfile containing:

- all original input variables (copied unchanged from ``infile``)
- the derived field ``outVarName = A OP B`` appended as the last component

The refinement ratios are read from the input plotfile (``amrData.RefRatio()``)
so that multi-level AMR data is written correctly even when the refinement
ratio is not 2.

Typical Applications
--------------------
- Computing dimensionless quantities (e.g. mixture fraction from species mass
  fractions)
- Weighting a scalar by a density or volume fraction
- Differencing two plotfiles of the same field at different times

Notes
-----
Both ``inVarAName`` and ``inVarBName`` must exist in the plotfile — the tool
aborts with an informative message if either is missing.

For ``operator=divide``, operand B is checked for zeros before the division
using a GPU-compatible ``ReduceOps<ReduceOpLogicalOr>`` reduction across all
AMR levels and MPI ranks.  Setting ``checkDivByZero=0`` skips this check; the
division then follows IEEE 754 semantics (producing ``Inf`` or ``NaN`` for
zero denominator cells without aborting).

The tool is compatible with MPI and GPU (CUDA/HIP/SYCL) builds.

Testing
-------
A self-contained test suite is provided in ``Tests/arithmetics/``.  See
:doc:`/testing/arithmetics` for the full test matrix.
