.. highlight:: bash

Testing
*******

Functional test suites live in the ``Tests/`` directory at the repository
root. Each subdirectory covers one tool and is self-contained: it provides
synthetic input generators, tool input files, an automated runner, and a
``TESTING.md`` with per-test expected outputs.

Running a suite
###############

::

   cd Tests/<toolname>
   ./run_tests.sh               # compile executables then run all tests
   ./run_tests.sh --no-compile  # skip build, use existing executables in Src/

The runner prints a colour-coded PASS/FAIL line for each assertion and exits
with a non-zero status if any test fails. All artefacts land in
``Tests/<toolname>/testrun/`` and are excluded from version control.

Available suites
################

.. list-table::
   :header-rows: 1
   :widths: 10 35 55

   * - Suite
     - Tool
     - Coverage summary
   * - :doc:`jpdf <testing/jpdf>`
     - Joint PDFs and 2D conditional means
     - 8 scenarios; 2D + 3D serial + 3D MPI (2 and 4 ranks)
   * - :doc:`arithmetics <testing/arithmetics>`
     - Binary field arithmetic (add, subtract, multiply, divide)
     - 15 assertions; round-trip value verification; 3D serial + MPI

.. toctree::
   :hidden:

   testing/jpdf
   testing/arithmetics

Adding a new suite
##################

1. Create ``Tests/<toolname>/`` with input files and a ``run_tests.sh``
   modelled on ``Tests/jpdf/run_tests.sh``.
2. Set ``SRC_DIR="$SCRIPT_DIR/../../Src"`` so the runner finds the built
   executables.
3. Write a ``TESTING.md`` documenting expected outputs for each test.
4. Add a row to the table above and a subpage under ``Docs/source/testing/``.
