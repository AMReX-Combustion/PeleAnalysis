.. highlight:: bash

avgToPlane
******************************************
Project an AMReX plot file along a chosen spatial direction and render the result as a PPM colour image. The tool sums all cell values along the selected axis to produce a 2D projection, maps the result through a user-supplied 256-colour binary palette, and writes a PPM image file. This provides a fast visual summary of 3D simulation data without requiring a full visualisation tool.

Usage: ::

   ./avgToPlane.gnu.MPI.ex infile=FILE palette=FILE [OPTIONS]

Example: ::

   ./avgToPlane.gnu.MPI.ex ./InputSamples/avgToPlane.inp

.. note::

   Only a single level (`finestLevel`) is used for the projection. Data from coarser levels is not included. The output file is always named ``<infile>.ppm`` and cannot be overridden.

.. note::

   `palette` is a required argument. The palette file must be a raw binary file containing exactly 256 red values, followed by 256 green values, followed by 256 blue values (768 bytes minimum). An optional fourth block of 256 alpha values may follow. Palette files compatible with VisIt or Amrvis can be used directly.

Tool Options
#############
::

   #------------------- IO CONTROL -----------------------------------------------------------
   infile  = plt00500                         # Input AMReX plot file
   palette = palettes/rainbow.pal             # Binary palette file (256 RGB triplets, required)

`infile` specifies the AMReX plot file to read. `palette` specifies the binary colour palette file used to map projected data values to RGB colours. Both are required. The output PPM image is automatically named ``<infile>.ppm``.
::

   #------------------- AMR Control ----------------------------------------------------------
   finestLevel = 2                            # DEF: finest level in file; Level used for projection

`finestLevel` sets the AMR level from which data is read for projection. Only this single level is used; data from other levels is ignored. Defaults to the finest level present in the input file.
::

   #------------------- Projection Options ---------------------------------------------------
   dir = 2                                    # DEF: 0; Direction to project along (0=x, 1=y, 2=z)

`dir` specifies the spatial axis along which the data is summed. All cells along this direction are accumulated into a single 2D image plane. Must be in the range ``[0, AMREX_SPACEDIM-1]``. Defaults to 0 (x-axis). The two remaining dimensions form the axes of the output image.
::

   #------------------- Spatial Subsetting ---------------------------------------------------
   box = 0 0 0 127 127 127                    # DEF: full domain; Index bounds at finestLevel (lo_x lo_y lo_z hi_x hi_y hi_z)

`box` restricts the projection to a spatial sub-region specified in index space at `finestLevel`. The format is ``lo_x lo_y lo_z hi_x hi_y hi_z`` for 3D, or ``lo_x lo_y hi_x hi_y`` for 2D. Defaults to the full problem domain.
::

   #------------------- Variable Selection ---------------------------------------------------
   comps = 3                                  # Specific component indices to render (overrides sComp/nComp)
   sComp = 0                                  # DEF: 0; Start component index (used if comps not set)
   nComp = 1                                  # DEF: all; Number of components (used if comps not set)

Variable selection follows the same convention as other tools in this suite. An explicit list of zero-based component indices can be provided via `comps`. If `comps` is not set, the contiguous range starting at `sComp` with length `nComp` is used. `comps` takes precedence if both are provided. Note that all selected components are projected and accumulated into a single result image; to render individual components separately, run the tool once per component.
::

   #------------------- Colour Scale ---------------------------------------------------------
   min = 0.0                                  # DEF: projected data minimum; Lower bound of colour scale
   max = 1.0                                  # DEF: projected data maximum; Upper bound of colour scale

`min` and `max` set the data range mapped to the full 256-level colour scale. Values below `min` are clamped to the first palette entry and values above `max` are clamped to the last. If not specified, the actual minimum and maximum of the projected data are used automatically. Setting these explicitly is recommended when comparing images across multiple plot files to ensure a consistent colour scale.
::

   #------------------- Grid Options ---------------------------------------------------------
   max_grid_size = 32                         # DEF: 32; Internal box size for parallel decomposition

`max_grid_size` controls the maximum box size used for the internal parallel decomposition of the data before projection. It does not affect the output image. Smaller values increase parallelism at the cost of overhead; the default of 32 is appropriate for most use cases.

Output
######
The output file ``<infile>.ppm`` is a standard 24-bit binary PPM (Portable Pixmap) image. It can be opened directly in most image viewers, or converted to PNG or other formats with standard tools:

.. code-block:: bash

   convert plt00500.ppm plt00500.png     # ImageMagick
   ffmpeg -i plt00500.ppm plt00500.png   # FFmpeg

The image dimensions are determined by the extent of the projected domain (or `box`) in the two non-projected directions.
