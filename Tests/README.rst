.. highlight:: bash

Tests
*****

Functional test suites for PeleAnalysis tools. Each subdirectory contains
synthetic plotfile generators, tool input files, an automated test runner,
and detailed documentation.

Running
#######

Each suite provides a ``run_tests.sh`` script that compiles the required
executables, generates synthetic plotfiles, and validates all test cases::

   cd Tests/jpdf
   ./run_tests.sh               # compile + run (requires AMReX build env)
   ./run_tests.sh --no-compile  # skip build, use existing executables in Src/

All test artefacts are written to a ``testrun/`` subdirectory inside the suite
and can be deleted freely.

Available suites
################

.. list-table::
   :header-rows: 1
   :widths: 20 35 45

   * - Directory
     - Tool
     - Coverage
   * - ``jpdf/``
     - Joint PDFs and 2D conditional means
     - 8 scenarios: all output formats, condMean, conditioning modes 0/1/2,
       norm_cVal, temporal averaging; 2D + 3D serial + 3D MPI (2 and 4 ranks)
   * - ``arithmetics/``
     - Binary field arithmetic (add, subtract, multiply, divide)
     - 15 assertions: all four operators, divide-by-zero detection/bypass,
       output naming, round-trip value verification, error handling, 3D MPI
       (2 and 4 ranks)

Adding a new suite
##################

1. Create ``Tests/<toolname>/`` with input files and a ``run_tests.sh``
   modelled on ``Tests/jpdf/run_tests.sh``.
2. Set ``SRC_DIR="$SCRIPT_DIR/../../Src"`` to point at the source tree.
3. Add a ``TESTING.md`` documenting expected outputs.
4. Register the suite in ``Docs/source/testing.rst``.
