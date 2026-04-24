.. highlight:: bash

arithmetics
***********

Covers the :doc:`arithmetics </analysis/arithmetics>` tool (binary field
arithmetic on AMReX plotfiles).

**Location:** ``Tests/arithmetics/``

.. list-table::
   :header-rows: 1
   :widths: 5 35 60

   * - #
     - Feature
     - Pass criterion
   * - 1
     - ``operator=add``
     - Exit 0; output plotfile exists; Header nComp = 3; "sum" in Header
   * - 2
     - ``operator=subtract``
     - Exit 0; output plotfile exists; Header nComp = 3; "diff" in Header
   * - 3
     - ``operator=multiply``
     - Exit 0; output plotfile exists; Header nComp = 3; "prod" in Header
   * - 4
     - ``operator=divide`` (non-zero denominator)
     - Exit 0; output plotfile exists; Header nComp = 3; "quot" in Header
   * - 5
     - Divide-by-zero detected (``checkDivByZero=1``)
     - Exits non-zero; no output plotfile written
   * - 6
     - Divide-by-zero bypassed (``checkDivByZero=0``)
     - Exit 0; output plotfile written despite zero denominator cells
   * - 7
     - Default output naming (``infile_operator``)
     - ``plt_arith_add/`` created without specifying ``outfile``
   * - 7b
     - ``outfile`` override
     - User-specified directory created instead of default name
   * - 8
     - Add/subtract round-trip: (A + B) − B = A
     - Recovered component equals original varA = 4.0 in every cell
   * - 9
     - Multiply/divide round-trip: (A × B) ÷ B = A
     - Recovered component equals original varA = 4.0 in every cell
   * - 10
     - ``inVarAName`` not found in plotfile
     - Exits non-zero with informative message
   * - 11
     - ``inVarBName`` not found in plotfile
     - Exits non-zero with informative message
   * - 12
     - Invalid ``operator`` value
     - Exits non-zero with informative message
   * - MPI-2
     - Add with 2 MPI ranks
     - ``sum`` = 6.0 and ``varA`` = 4.0 in all cells
   * - MPI-4
     - Add with 4 MPI ranks
     - ``sum`` = 6.0 and ``varA`` = 4.0 in all cells
   * - MPI-div
     - Divide-by-zero with 4 MPI ranks
     - All ranks abort (``mpirun`` exits non-zero); no MPI deadlock

Value correctness (tests 8–9 and MPI tests) is verified by
``check_plt_value.py``, a self-contained Python 3 script that reads AMReX
binary FAB files directly and checks every cell against an expected value.
See ``Tests/arithmetics/TESTING.md`` for the full test matrix, expected
outputs, and known pitfalls.
