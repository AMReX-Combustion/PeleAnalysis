.. highlight:: bash

plotYTtoLe
**********
Compute per-species Lewis numbers from an AMReX plot file produced by a reacting flow solver such as PeleLMeX. The tool reads species mass fractions, temperature, and density from each AMR level, calls the PelePhysics transport library to evaluate mixture transport coefficients, and derives the Lewis number for every species at every grid point. Results are written as a new AMReX plot file.

The Lewis number for species *n* is defined as:

.. math::

   \mathrm{Le}(n) = \frac{D_n}{\lambda / C_{p,\mathrm{mix}}}

where :math:`D_n` is the species mass diffusivity, :math:`\lambda` is the mixture thermal conductivity, and :math:`C_{p,\mathrm{mix}}` is the mixture heat capacity at constant pressure. A Lewis number of unity indicates that the species diffuses at the same rate as heat. Transport coefficients are evaluated using the PelePhysics transport library and the computation is GPU-accelerated when built with GPU support.

Usage: ::

   ./plotYTtoLe.gnu.MPI.ex infile=<s> [options]

Example: ::

   ./plotYTtoLe.gnu.MPI.ex ./InputSamples/plotYTtoLe.inp

.. note::

   The input plot file must contain the variables ``Y(<species>)`` for all species in the mechanism, ``temp``, and ``density``. If any of these are missing the tool will print a warning and produce incorrect results. The output file is always written to ``<infile>_Le`` and the name cannot be overridden.

Tool Options
#############
::

   #------------------- IO CONTROL -----------------------------------------------------------
   infile = plt00500                          # Input AMReX plot file

`infile` specifies the AMReX plot file to read. The file must have been produced by a PeleLMeX-compatible solver and must contain all species mass fractions ``Y(<species>)``, the temperature field ``temp``, and the density field ``density``. The output plot file is automatically named ``<infile>_Le`` and will contain one component ``Le(<species>)`` for each species in the chemical mechanism, in the same order as the mechanism species list.
::

   #------------------- AMR Control ----------------------------------------------------------
   finestLevel = 2                            # DEF: finest level in file; Finest AMR level to process

`finestLevel` sets the finest AMR level to include in the computation. Defaults to the finest level present in the input file. All levels from 0 up to and including `finestLevel` are processed and written to the output. Note that the refinement ratio between levels is assumed to be 2 in all spatial directions.
::

   #------------------- Additional Flags -----------------------------------------------------
   verbose                                    # Enable verbose output during data loading

The flag `verbose` enables additional console output from the AMReX data services layer during file reading, reporting grid and variable information as each level is loaded. This flag takes no value and is activated simply by its presence in the input file or on the command line.
::

   #------------------- Auxiliary variables --------------------------------------------------
   Aux_Variables = density                    # DEF: none; variables copied unchanged to output

`Aux_Variables` lists variables that are copied unchanged from the input to the output plot file, in addition to the computed ``Le(<species>)`` fields. This is useful for carrying through fields (such as ``density`` or ``temp``) that the tool does not otherwise write. The variables must exist in the input plot file. Default: none.
::

   #------------------- Output control -------------------------------------------------------
   n_files = 64                               # DEF: AMReX default; cap on the number of plotfile data files

`n_files` caps the number of binary files used to write the output plot file data (AMReX ``VisMF::SetNOutFiles``). Lower it to reduce the number of files created for large parallel post-processing runs. AMReX clamps the value to the number of MPI ranks, so a serial run always writes a single data file. Default: the AMReX default.

Output
######
The output plot file ``<infile>_Le`` is a standard multi-level AMReX plot file and can be visualised with any AMReX-compatible tool such as VisIt or ParaView (via the Boxlib/AMReX reader). It contains one variable per chemical species:

.. code-block:: none

   Le(H2)
   Le(H)
   Le(O)
   ...

The variables are ordered according to the species ordering defined in the compiled chemical mechanism (``mechanism.H``), followed by any ``Aux_Variables`` requested. The output inherits the domain geometry, coordinate system, and box structure from the input plot file.

Dependencies
############
This tool requires a PelePhysics-enabled build. The chemical mechanism and transport model are compiled in at build time and determine the species list, the number of output components, and the transport coefficient evaluation. The tool cannot be used with a generic AMReX build.
