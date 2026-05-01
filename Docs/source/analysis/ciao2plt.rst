.. highlight:: bash

ciao2plt
********

Convert a CIAO (NGA2) HDF5 file to a standard AMReX single-level plotfile.
The output can be used directly with all other PeleAnalysis tools.
The tool discovers the mesh group name at runtime, reads time and periodicity
from the file when available, and can auto-detect all available field names.
Both 2D and 3D files are supported via the compile-time ``DIM`` flag.

Usage: ::

   ./ciao2plt3d.gnu.ex infile=FILE [OPTIONS]

Example: ::

   ./ciao2plt3d.gnu.ex ./InputSamples/ciao2plt.inp

Build Configuration
####################

This tool links against the HDF5 C++ library directly. Two variables control
the HDF5 installation used:

``HDF5_DIR``
   Root of the **serial** HDF5 installation. Used for ``USE_MPI=FALSE`` builds
   (the default).

``HDF5_PAR_DIR``
   Root of a **parallel** HDF5 installation (compiled with MPI support). Used
   for ``USE_MPI=TRUE`` builds. Falls back to ``HDF5_DIR`` if not set.

On most HPC systems, serial and parallel HDF5 use the same library names
(``libhdf5``, ``libhdf5_cpp``) but live under different module paths. On
systems that give the parallel build distinct names (e.g., ``libhdf5_par``),
override ``LIBRARIES`` manually in ``GNUmakefile``.

.. code-block:: makefile

   # Serial
   make EBASE=ciao2plt HDF5_DIR=/path/to/hdf5-serial

   # MPI
   make EBASE=ciao2plt USE_MPI=TRUE HDF5_PAR_DIR=/path/to/hdf5-parallel

Tool Options
#############
::

   #------------------- IO CONTROL -----------------------------------------------------------
   infile  = SJ_H2_000000.h5                  # Input CIAO HDF5 file (required)
   outfile = plt_SJ_H2_000000                 # DEF: plt_<stem>; output plotfile directory

``infile`` is the only required argument. ``outfile`` defaults to ``plt_`` followed
by the input filename stem (extension stripped, directory stripped), e.g.
``SJ_H2_000000.h5`` → ``plt_SJ_H2_000000``.
::

   #------------------- Field selection -------------------------------------------------------
   vars = T RHO OH H2O                        # DEF: all fields

``vars`` is a space-separated list of field names to read and include in the
output plotfile. When omitted the tool enumerates all datasets in the
``cv_data_real`` group and all scalar components in ``scalars/SC`` and converts
everything. If any requested name is not found the tool aborts with a list of
available fields.
::

   #------------------- Grid control ----------------------------------------------------------
   max_grid_size = 32                         # DEF: 32

``max_grid_size`` controls how the single-level domain is decomposed into
AMReX boxes for parallel I/O. Larger values produce fewer, larger boxes.
::

   #------------------- Periodicity fallback --------------------------------------------------
   per = 0 0 0                                # DEF: 0 0 0

``per`` sets the periodicity flags used to construct the AMReX ``Geometry``
object. This parameter is only read when the HDF5 file does not contain
``geometry/sd_info/xper|yper|zper``; when those datasets are present they
take precedence.
::

   #------------------- Coordinate system fallback --------------------------------------------
   coord_sys = 0                              # DEF: 0 (Cartesian); 1 = cylindrical/RZ

``coord_sys`` sets the AMReX coordinate system type. ``0`` is Cartesian
(default); ``1`` is cylindrical/RZ, which requires a ``DIM=2`` build.
This parameter is only read when the HDF5 file does not contain
``geometry/sd_info/icyl``; when that dataset is present it takes precedence.

HDF5 File Structure
####################

The tool expects the standard NGA2/CIAO HDF5 layout:

.. code-block:: text

   /
   ├── IO-information/           (skipped)
   └── <mesh_group>/             (name discovered at runtime)
       ├── data/
       │   ├── cv_data_real/     individual 3-D datasets, one per field
       │   ├── scalars/SC        4-D array (nscalars, nz, ny, nx) with
       │   │                     field names stored as attributes "Index 1" …
       │   └── globals_r0/time   simulation time (required, F64 scalar)
       └── geometry/
           ├── grid/x, y, z      node coordinates (N+1 values for N cells)
           └── sd_info/          optional; xper, yper, zper periodicity flags;
                              icyl coordinate type (0=Cartesian, 1=cylindrical)

Grid coordinates are node-based: a coordinate array of length *N*+1 defines
*N* cells. Field data in ``cv_data_real`` is stored in C order ``(nz, ny, nx)``
and may be single (``F32``) or double (``F64``) precision; the tool always
writes double-precision AMReX output.

.. note::

   The mesh group name (e.g. ``sd_mesh_SJ_H2``, ``sd_box``) is discovered
   automatically. Exactly one mesh group must be present; the tool aborts if
   zero or more than one are found.

2D Files
########

Compile with ``DIM=2`` to process files where the ``z`` grid has exactly 2
nodes (1 cell). The tool reads the ``x`` and ``y`` coordinates only and
creates a 2-D AMReX plotfile. It aborts if the file has more than 1 z-cell
when built in 2D mode.

.. code-block:: bash

   make EBASE=ciao2plt DIM=2
   ./ciao2plt2d.gnu.ex infile=data2d.h5

Output
######

A single-level AMReX plotfile directory (``plt_<stem>/`` by default)
containing a ``Header`` file and cell-centred data. The plotfile can be
opened with VisIt, ParaView (AMReX/Boxlib reader), or any PeleAnalysis
post-processing tool. Time is read from ``globals_r0/time`` and stored in
the plotfile header.
