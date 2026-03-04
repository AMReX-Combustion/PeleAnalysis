.. highlight:: bash

plotQPD
*******
Compute a Quasi-Planar Diagram (QPD) of atomic mass fluxes through a chemical mechanism from an AMReX plot file produced by a reacting flow solver such as PeleLMeX. The tool integrates forward and reverse reaction rates across the entire domain, maps them onto species-to-species edges for a chosen atom, normalises the fluxes, and writes the result as a plain-text edge-list file suitable for QPD visualisation tools.

QPDs are a standard diagnostic in combustion analysis, used to identify dominant reaction pathways, quantify how an atom (typically carbon) flows between species, and determine which intermediates play the most significant role in fuel breakdown. Each edge in the diagram connects two species and carries the net flux of the tracked atom between them, summed over all reactions that link that pair.

Normalisation is performed relative to the CH4 → CH3 edge (the primary fuel destruction step in methane combustion) by default, so all fluxes are expressed as fractions of that reference flux. An optional additional scale factor can be applied on top of this.

Usage: ::

   ./plotQPD.gnu.MPI.ex infile=<s> [options]

Example: ::

   ./plotQPD.gnu.MPI.ex ./InputSamples/plotQPD.inp

.. note::

   The input plot file must contain mole fractions ``X(<species>)`` for all species in the mechanism, the temperature field ``temp``, and the density field ``density``. Mass fractions ``Y(...)`` are **not** sufficient — this tool requires mole fractions. If any required field is missing the tool will print a warning and produce incorrect results.

.. warning::

   OpenMP threading is currently disabled in this tool and will cause an abort if enabled at build time. MPI parallelism is fully supported.

Tool Options
#############
::

   #------------------- IO CONTROL -----------------------------------------------------------
   infile      = plt00500                     # Input AMReX plot file
   QPDfileName = plt00500_QPD.dat             # DEF: <infile>_QPD.dat; Output QPD data file
   QPDlabel    = plt00500                     # DEF: <infile>; Label written to header of output file

`infile` specifies the AMReX plot file to read. The file must contain mole fractions ``X(<species>)`` for all mechanism species, ``temp``, and ``density``. `QPDfileName` sets the name of the output edge-list file; it defaults to ``<infile>_QPD.dat``. `QPDlabel` sets a descriptive label written to the first line of the output file header, which can be used by downstream visualisation tools to identify the data source; it defaults to the input file name.
::

   #------------------- AMR Control ----------------------------------------------------------
   finestLevel = 2                            # DEF: finest level in file; Finest AMR level to process

`finestLevel` sets the finest AMR level to include in the integration. Defaults to the finest level present in the input file. Data covered by finer levels is automatically zeroed out at coarser levels to prevent double-counting in the volume-weighted summation.
::

   #------------------- QPD Options ----------------------------------------------------------
   QPDatom   = C                              # DEF: C; Atom to track through the reaction network
   scaleNorm = 1.0                            # Optional additional scale factor for normalisation

`QPDatom` specifies which atom to track through the mechanism. The edge list is constructed by finding all pairs of species that exchange this atom through one or more reactions. Defaults to ``C`` (carbon), which is the standard choice for hydrocarbon flame analysis. Other common choices include ``H``, ``O``, or ``N``.

All fluxes are normalised relative to the net flux on the CH4 → CH3 edge (the dominant fuel consumption step). If `scaleNorm` is provided, the normalisation value is additionally multiplied by this factor, allowing fluxes to be expressed in physical units or rescaled to a different reference reaction.
::

   #------------------- Fuel Species Diagnostics ---------------------------------------------
   fuelSpec = CH4                             # If set, prints reaction partner breakdown to screen

If `fuelSpec` is set, the tool prints a detailed breakdown to the screen for every edge involving that species. For each such edge, the contribution of each reaction partner is listed along with the sum of positive and negative contributions. This is useful for identifying which co-reactants most strongly drive consumption or production of the fuel species. This output goes to stdout only and is not written to `QPDfileName`.
::

   #------------------- Additional Flags -----------------------------------------------------
   dump_edges                                 # Print all edges to screen before computation
   verbose                                    # Enable verbose AMReX data loading output

`dump_edges` causes the full list of atom-transfer edges to be printed to stdout before the flux computation begins. This is useful for verifying that the mechanism edge graph has been constructed correctly for the chosen atom. `verbose` enables additional AMReX data services output during file reading.

Output
######
The output file ``<infile>_QPD.dat`` (or the name set via `QPDfileName`) is a plain-text file with the following structure:

.. code-block:: none

   <QPDlabel>
   H2 H O O2 OH H2O ...     [space-separated species list]
   CH4 CH3 <Qf> <-Qr>
   CH4 CO  <Qf> <-Qr>
   CH3 CH2 <Qf> <-Qr>
   ...

Each data line contains the left species, the right species, the normalised forward flux, and the negated normalised reverse flux. The sign convention is such that a positive forward flux represents net transfer from left to right. All values are normalised to the CH4 → CH3 reference edge unless `scaleNorm` modifies the normalisation.

This format is compatible with standard QPD visualisation tools. The species list in the header allows downstream tools to map species names to indices.

Relation to Other Tools
########################
- This tool requires **mole fractions** ``X(...)`` as input, whereas :doc:`plotTYtoLe` and :doc:`plotTransportCoeff` require **mass fractions** ``Y(...)``. Ensure the correct variable names are present in the plot file before running.
- The reaction rate evaluation uses the Chemkin-style routines ``CKPX`` (pressure from mole fractions) and ``CKKFKR`` (forward and reverse rate constants), which are part of the compiled chemical mechanism.

Dependencies
############
This tool requires a PelePhysics-enabled build with a compiled chemical mechanism. The species list, reaction count, reaction map, and edge graph are all determined at compile time from ``mechanism.H``. The tool cannot be used with a generic AMReX build.
