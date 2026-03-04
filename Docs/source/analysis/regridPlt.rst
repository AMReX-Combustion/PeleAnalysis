.. highlight:: bash

regridPlt
******************************************
Regrid an AMReX plot file to a uniform box decomposition. The tool reads a plot file, rebuilds the ``BoxArray`` at each AMR level by splitting boxes to a specified maximum size, and writes the result as a new plot file. This is useful for improving load balance for downstream analysis tools, reducing per-box memory pressure, or standardising the grid layout across a set of plot files.

A subset of variables can be selected for output, and periodicity flags can be set independently of the original plot file.

Usage: ::

   ./regridPlt.gnu.MPI.ex infile=FILE outfile=FILE [OPTIONS]

Example: ::

   ./regridPlt.gnu.MPI.ex ./InputSamples/regridPlt.inp

Tool Options
#############
::

   #------------------- IO CONTROL -----------------------------------------------------------
   infile  = plt00500                         # Input AMReX plot file
   outfile = plt00500_rg                      # Output regridded plot file

`infile` and `outfile` are both required. `outfile` receives no default and must always be specified explicitly. The output is a standard multi-level AMReX plot file preserving the geometry, refinement ratios, and simulation time of the input.
::

   #------------------- AMR Control ----------------------------------------------------------
   finestLevel = 2                            # DEF: all levels in file; Finest AMR level to regrid

`finestLevel` sets the finest AMR level to include in the output. Defaults to the finest level present in the input file. All levels from 0 to `finestLevel` are regridded and written.
::

   #------------------- Grid Options ---------------------------------------------------------
   max_grid_size = 128                        # DEF: 128; Maximum box size in each dimension

`max_grid_size` controls the maximum extent of any single box in the output ``BoxArray`` in each spatial dimension. Any box in the original layout that exceeds this size is split into uniform tiles of at most ``max_grid_size`` cells. Smaller values produce more, smaller boxes and can improve parallelism and memory locality at the cost of increased overhead. The default of 128 is a reasonable starting point for most use cases.
::

   #------------------- Variable Selection ---------------------------------------------------
   comps = 0 1 4 7                            # Specific component indices to write (overrides sComp/nComp)
   sComp = 0                                  # DEF: 0; Start component index (used if comps not set)
   nComp = 4                                  # DEF: all; Number of components (used if comps not set)

Variable selection follows the same convention as other tools in this suite. An explicit list of zero-based component indices can be provided via `comps`. If `comps` is not set, the contiguous range starting at `sComp` with length `nComp` is used instead, defaulting to all components in the input file. `comps` takes precedence over `sComp`/`nComp` if both are provided.
::

   #------------------- Domain Options -------------------------------------------------------
   is_per = 1 1 0                             # DEF: 1 1 1; Periodicity flags per dimension

`is_per` sets the periodicity of the domain in each spatial dimension (1 = periodic, 0 = non-periodic). Provide one value per spatial dimension. Defaults to fully periodic (all 1s). Note that this flag is set independently of whatever periodicity was used in the original simulation — it affects only the geometry metadata written to the output file.
::

   #------------------- Additional Flags -----------------------------------------------------
   verbose                                    # Enable verbose output during data loading

`verbose` enables per-variable progress output during data loading, reporting each variable name as it is read from the input file and again after it is flushed from memory. Useful for monitoring progress on large plot files.
