.. highlight:: bash

subPlt
******

Extract a spatial subregion and/or variable subset from an AMReX plot file and write the result as a new, self-consistent plot file. The tool correctly handles multi-level AMR data by coarsening and refining the requested box across all levels to respect the refinement hierarchy, ensuring the output is a valid AMReX plot file that can be used directly in downstream tools.

Usage: ::

   ./subPlt.gnu.MPI.ex infile=<s> [outfile=<s>] [options]

Example: ::

   ./subPlt.gnu.MPI.ex ./InputSamples/subPlt.inp

Tool Options
#############

::

   #------------------- IO CONTROL -----------------------------------------------------------
   infile      = plt00500                     # Input AMReX plot file
   outfile     = plt00500_section             # DEF: <infile>_section; Output plot file name

`infile` specifies the AMReX plot file to read. `outfile` specifies the name of the extracted output plot file. If `outfile` is not provided, the output is named ``<infile>_section`` by default.

::

   #------------------- AMR Control ----------------------------------------------------------
   finestLevel = 2                            # DEF: finest level in file; Finest AMR level to include

`finestLevel` sets the finest AMR level to include in the extraction. Defaults to the finest level present in the input file. All levels from 0 up to and including `finestLevel` are written to the output. Reducing `finestLevel` can significantly reduce output file size for highly refined datasets.

::

   #------------------- Spatial Subsetting ---------------------------------------------------
   box = 32 0 0 64 128 128                    # DEF: full domain; lo and hi index bounds (ix,iy,iz iX,iY,iZ)

`box` defines the spatial subregion to extract as index-space bounds at the finest requested level. The format is ``ix iy iz iX iY iZ``, where ``(ix, iy, iz)`` is the low corner and ``(iX, iY, iZ)`` is the high corner. For ``AMREX_SPACEDIM=2``, only four values are required: ``ix iy iX iY``. The box is automatically coarsened and refined to produce consistent bounds at every AMR level. If `box` is not specified, the full problem domain is used.

::

   #------------------- Variable Selection ---------------------------------------------------
   comps = 0 1 4 7                            # Specific component indices to extract (overrides sComp/nComp)
   sComp = 0                                  # DEF: 0; Start component index (used if comps not set)
   nComp = 4                                  # DEF: all; Number of components to extract (used if comps not set)

There are two ways to select which variables to extract. The first is to provide an explicit list of zero-based component indices via `comps`. This allows selecting an arbitrary, non-contiguous set of variables (e.g. ``comps = 0 1 4 7``). If `comps` is not set, the contiguous range starting at `sComp` with length `nComp` is used instead. `sComp` defaults to 0 and `nComp` defaults to the total number of components in the input file, meaning all variables are extracted by default. `comps` and the `sComp`/`nComp` pair are mutually exclusive; if `comps` is provided it takes precedence.

::

   #------------------- Additional Flags -----------------------------------------------------
   verbose                                    # Enable verbose output during data loading

The flag `verbose` enables additional output during execution, reporting which variables are being filled on each AMR level. This is useful for monitoring progress on large plot files or for debugging unexpected output.
