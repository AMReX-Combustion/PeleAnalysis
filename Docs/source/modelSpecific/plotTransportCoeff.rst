.. highlight:: bash

plotTransportCoeff
******************
.. warning::
   The compiling process requires absolute paths to PelePhysics in the GNUmakefile.

Evaluate and output all mixture transport coefficients from an AMReX plot file produced by a reacting flow solver such as PeleLMeX. The tool reads species mass fractions, temperature, and density at each AMR level, calls the PelePhysics transport library to compute the full set of transport coefficients, and writes the results as a new AMReX plot file. This tool is closely related to :doc:`plotTYtoLe`, which derives Lewis numbers from the same transport coefficients; here the raw coefficients are written directly without further post-processing.

The following transport coefficients are computed and written for each grid point:

.. list-table::
   :widths: 25 75
   :header-rows: 1

   * - Variable
     - Description
   * - ``rhoD(<species>)``
     - Density-weighted mass diffusivity :math:`\rho D_n` for each species
   * - ``chi(<species>)``
     - Soret (thermal diffusion) coefficient :math:`\chi_n` for each species
   * - ``mu``
     - Dynamic (shear) viscosity of the mixture
   * - ``xi``
     - Bulk viscosity of the mixture
   * - ``lambda``
     - Thermal conductivity of the mixture

Transport coefficients are evaluated using the PelePhysics transport library and the computation is GPU-accelerated when built with GPU support.

Usage: ::

   ./plotTransportCoeff.gnu.MPI.ex infile=<s> [options]

Example: ::

   ./plotTransportCoeff.gnu.MPI.ex ./InputSamples/plotTransportCoeff.inp

.. note::

   The input plot file must contain the variables ``Y(<species>)`` for all species in the mechanism, ``temp``, and ``density``. If any of these are missing the tool will print a warning and produce incorrect results. The output file is always written to ``<infile>_D`` and the name cannot be overridden.

Tool Options
#############
::

   #------------------- IO CONTROL -----------------------------------------------------------
   infile = plt00500                          # Input AMReX plot file

`infile` specifies the AMReX plot file to read. The file must have been produced by a PeleLMeX-compatible solver and must contain all species mass fractions ``Y(<species>)``, the temperature field ``temp``, and the density field ``density``. The output plot file is automatically named ``<infile>_D`` and contains ``2 * NUM_SPECIES + 3`` components: density-weighted diffusivities and Soret coefficients for each species, followed by the three scalar mixture properties ``mu``, ``xi``, and ``lambda``.
::

   #------------------- AMR Control ----------------------------------------------------------
   finestLevel = 2                            # DEF: finest level in file; Finest AMR level to process

`finestLevel` sets the finest AMR level to include in the computation. Defaults to the finest level present in the input file. All levels from 0 up to and including `finestLevel` are processed and written to the output. The refinement ratio between levels is assumed to be 2 in all spatial directions.
::

   #------------------- Additional Flags -----------------------------------------------------
   verbose                                    # Enable verbose output during data loading

The flag `verbose` enables additional console output from the AMReX data services layer during file reading. It takes no value and is activated simply by its presence in the input file or on the command line.

Output
######
The output plot file ``<infile>_D`` is a standard multi-level AMReX plot file and can be visualised with any AMReX-compatible tool such as VisIt or ParaView (via the Boxlib/AMReX reader). The variables are written in the following order:

.. code-block:: none

   rhoD(H2)
   rhoD(H)
   rhoD(O)
   ...                    [one rhoD per species, in mechanism order]
   chi(H2)
   chi(H)
   chi(O)
   ...                    [one chi per species, in mechanism order]
   mu
   xi
   lambda

The output inherits the domain geometry, coordinate system, and box structure from the input plot file.

Relation to Other Tools
########################
- :doc:`plotTYtoLe` uses the same inputs and calls the same transport library, but derives Lewis numbers :math:`\mathrm{Le}(n) = D_n / (\lambda / C_{p,\mathrm{mix}})` rather than writing the raw coefficients.
- The ``rhoD`` values written by this tool are the density-weighted diffusivities :math:`\rho D_n`. To recover the mass diffusivity :math:`D_n` alone, divide by the ``density`` field from the original plot file.

Dependencies
############
This tool requires a PelePhysics-enabled build. The chemical mechanism and transport model are compiled in at build time and determine the species list, the number of output components, and the transport coefficient evaluation. The tool cannot be used with a generic AMReX build.
