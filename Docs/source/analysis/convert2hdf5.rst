.. highlight:: bash

convert2hdf5
************

Convert a standard AMReX plot file to HDF5 format. All variables and AMR levels are written into a single HDF5 dataset using AMReX's parallel HDF5 I/O layer. Optional compression is supported via ZLIB (lossless) or ZFP (lossy), controlled through the H5Z-ZFP filter. The output file is compatible with AMReX's HDF5 plotfile reader and can be used with VisIt or ParaView via the Boxlib/AMReX reader.

Usage: ::

   ./convert2hdf5.gnu.MPI.ex infile=FILE [OPTIONS]

Example: ::

   ./convert2hdf5.gnu.MPI.ex ./InputSamples/convert2hdf5.inp

.. note::

   The output file is always written to ``<infile>_hdf5`` and the name cannot be overridden. All variables present in the input plot file are converted; there is no option to select a subset of components.

Build Configuration
####################
This tool requires HDF5 support to be enabled at compile time. The relevant settings are in ``GNUmakefile``:

.. code-block:: makefile

   ifeq ($(EBASE),convert2hdf5)
      USE_HDF5     = TRUE
      USE_HDF5_ZFP = TRUE          # Set to FALSE to disable ZFP lossy compression support
      HDF5_HOME    = /path/to/hdf5 # Replace with your HDF5 installation path
      ZFP_HOME     = /path/to/zfp  # Replace with your ZFP installation path (only needed if USE_HDF5_ZFP=TRUE)
      H5Z_HOME     = /path/to/h5z-zfp # Replace with your H5Z-ZFP installation path (only needed if USE_HDF5_ZFP=TRUE)
   endif

The three path variables must be set to the actual installation locations on your system. ``USE_HDF5_ZFP`` can be set to ``FALSE`` if ZFP lossy compression is not needed, in which case ``ZFP_HOME`` and ``H5Z_HOME`` do not need to be set. ``USE_HDF5 = TRUE`` is always required for this tool.

Tool Options
#############
::

   #------------------- IO CONTROL -----------------------------------------------------------
   infile = plt00500                          # Input AMReX plot file

`infile` specifies the AMReX plot file to read. All variables present in the file are read and converted. The output is written to ``<infile>_hdf5`` automatically.
::

   #------------------- AMR Control ----------------------------------------------------------
   finestLevel = 2                            # DEF: all levels in file; Finest AMR level to convert

`finestLevel` sets the finest AMR level to include in the conversion. Defaults to the finest level present in the input file. Coarsening the output by reducing `finestLevel` can significantly reduce file size for highly refined datasets.
::

   #------------------- HDF5 Options ---------------------------------------------------------
   hdf5_compression = ZLIB@9                  # DEF: ZLIB@9; Compression filter and level

`hdf5_compression` controls the compression filter applied when writing the HDF5 dataset. The default is ``ZLIB@9``, which applies lossless deflate compression at the maximum compression level. The available options are:

.. list-table::
   :widths: 30 70
   :header-rows: 1

   * - Value
     - Description
   * - ``ZLIB@N``
     - Lossless deflate (zlib) compression. ``N`` sets the compression level in the range ``[1, 9]``, where 1 is fastest and 9 gives the smallest file. Always available when ``USE_HDF5 = TRUE``.
   * - ``ZFP@rate``
     - Lossy ZFP compression at the specified bit rate. Smaller rates give higher compression with greater data loss. Requires ``USE_HDF5_ZFP = TRUE`` in ``GNUmakefile`` and the ZFP and H5Z-ZFP libraries.
   * - ``NONE``
     - No compression. Produces the largest output file but requires no additional libraries beyond HDF5.

.. warning::

   ZFP compression (``ZFP@rate``) is lossy. Data values in the output file will differ from the input by an amount determined by the specified bit rate. Do not use ZFP compression if exact data reproduction is required.

Output
######
The output file ``<infile>_hdf5`` is a multi-level AMReX HDF5 plotfile stored as a single dataset. It contains all variables from the input file across all converted AMR levels, and can be opened with any AMReX-compatible HDF5 reader. Run time is printed to stdout on completion.

Dependencies
############
This tool requires an HDF5-enabled AMReX build. ``USE_HDF5 = TRUE`` must be set in ``GNUmakefile`` and a compatible HDF5 installation must be available at ``HDF5_HOME``. ZFP compression additionally requires the ZFP and H5Z-ZFP libraries and ``USE_HDF5_ZFP = TRUE``.
