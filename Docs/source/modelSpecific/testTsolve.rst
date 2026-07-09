.. highlight:: bash

testTsolve
******************************************
Verify the consistency of the equation of state (EOS) temperature solver across an AMReX plot file produced by a reacting flow solver such as PeleLMeX. For each cell the tool computes the mixture enthalpy from the stored temperature and species mass fractions, then solves back for temperature from that enthalpy, and writes both the recovered temperature and the residual difference to a new plot file.

A non-zero residual field ``dtemp`` indicates a discrepancy between the temperature stored in the plot file and the value recovered by the EOS solver, which can arise from solver convergence issues, inconsistent thermodynamic data, or numerical drift in the simulation.

Usage: ::

   ./testTsolve.gnu.MPI.ex infile=FILE [OPTIONS]

Example: ::

   ./testTsolve.gnu.MPI.ex ./InputSamples/testTsolve.inp

.. note::

   The input plot file must contain mass fractions ``Y(<species>)`` for all species in the mechanism and the temperature field ``temp``. The output file is always written to ``<infile>_T`` and cannot be overridden.

.. note::

   The EOS temperature solve (``HY2T``) uses a fixed initial guess of 300 K. In regions where the true temperature is far from this initial guess, convergence behaviour may affect the recovered value and the residual.

Tool Options
#############
::

   #------------------- IO CONTROL -----------------------------------------------------------
   infile = plt00500                          # Input AMReX plot file (must contain Y(...) and temp)

`infile` specifies the AMReX plot file to read. The file must have been produced by a PeleLMeX-compatible solver and must contain all species mass fractions ``Y(<species>)`` and the temperature field ``temp``. The output plot file is automatically named ``<infile>_T``.
::

   #------------------- AMR Control ----------------------------------------------------------
   finestLevel = 2                            # DEF: finest level in file; Finest AMR level to process

`finestLevel` sets the finest AMR level to include in the computation. Defaults to the finest level present in the input file. All levels from 0 up to and including `finestLevel` are processed and written to the output.
::

   #------------------- Additional Flags -----------------------------------------------------
   verbose                                    # Enable verbose output during data loading

`verbose` enables additional console output from the AMReX data services layer during file reading. It takes no value and is activated simply by its presence in the input file or on the command line.
::

   #------------------- Auxiliary variables --------------------------------------------------
   Aux_Variables = density                    # DEF: none; variables copied unchanged to output

`Aux_Variables` lists variables that are copied unchanged from the input to the output plot file, appended after ``temp`` and ``dtemp``. This is useful for carrying through fields (such as ``density``) that the tool does not otherwise write. The variables must exist in the input plot file. Default: none.

Output
######
The output plot file ``<infile>_T`` contains two components at every grid point across all processed AMR levels:

.. list-table::
   :widths: 20 80
   :header-rows: 1

   * - Variable
     - Description
   * - ``temp``
     - Temperature recovered by solving ``HY2T(TY2H(T, Y), Y)``
   * - ``dtemp``
     - Residual: recovered temperature minus input temperature

A spatially uniform ``dtemp`` of zero confirms full EOS self-consistency. Large or spatially structured residuals indicate regions where the stored temperature is not consistent with the stored species composition under the compiled EOS.

.. note::

   ``testTsolve`` writes its output through the legacy ``WritePlotFile`` writer, which always writes one data file per MPI rank. The ``n_files`` option available in most other tools therefore does not apply here.

