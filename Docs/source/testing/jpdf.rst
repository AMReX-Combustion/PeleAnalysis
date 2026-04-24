.. highlight:: bash

jpdf
****

Covers the :doc:`jpdf </analysis/jpdf>` tool (Joint PDFs and 2D conditional means).

**Location:** ``Tests/jpdf/``

.. list-table::
   :header-rows: 1
   :widths: 5 30 65

   * - #
     - Feature
     - Pass criterion
   * - 1
     - All six output formats (plotfile, gnuplot, MATLAB, Tecplot, FAB, scatter)
     - All file types present; PDF sum = 1.0
   * - 2
     - 2D conditional mean (independent variable)
     - ``condMean_var_cond_on_*`` file created; all non-zero entries ≈ 0.5
   * - 3
     - Duplicate variable in ``condMean_vars``
     - Tool runs without error; deduplication confirmed in verbose output
   * - 4
     - ``useminmax`` range override + bin clamping
     - Header encodes overridden axis range; verbose reports ``v1g > 0``
   * - 5
     - ``do_conditioning=1`` (range filter)
     - PDF restricted to conditioned cells; sum = 1.0
   * - 6
     - ``do_conditioning=2`` (c(1−c) filter)
     - Tail bins of conditioning variable excluded; sum = 1.0
   * - 7
     - ``norm_cVal=1`` (normalised conditioning)
     - Only normalised-range cells contribute; sum = 1.0
   * - 8
     - Temporal averaging (``do_average=1``)
     - ``JPDFAverage*/`` created; averaged PDF matches per-file PDF exactly

Each test runs in both **2D** (tests 1, 2, 5, 8) and **3D** serial, and
tests 1, 5, 8 are additionally validated with **MPI** at 2 and 4 ranks.
MPI results are compared element-wise to serial (max diff < 10⁻¹⁰).
