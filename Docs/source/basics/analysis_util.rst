.. highlight:: c++

analysis_util
***************************

``analysis_util`` is a C++ utility library that wraps low-level AMReX
plotfile and MEF I/O routines behind a concise, MPI-safe interface.
All symbols live in the ``analysis_util`` namespace.
Include the header with:

::

   #include <analysis_util.H>

For new tools, make sure that the library is flagged for compilation in your ``GNUmakefile`` by adding the tool to the 
space separated list:

::

  ifeq ($(EBASE),$(filter $(EBASE), template <your_tool_name>)) 
    USE_UTILS = TRUE
  endif


Data Structures
---------------

The library defines two main data structures for holding plotfile and MEF data in memory.

PlotfileData
~~~~~~~~~~~~

Returned by ``read_plotfile``.  Holds everything needed to work with
a multi-level dataset and to call ``write_plotfile``.

.. list-table::
   :header-rows: 1
   :widths: 25 20 55

   * - Field
     - Type
     - Description
   * - ``mf``
     - ``Vector<MultiFab>``
     - One ``MultiFab`` per AMR level; components correspond to the requested variables.
   * - ``geoms``
     - ``Vector<Geometry>``
     - Domain geometry for each level.
   * - ``ref_ratios``
     - ``Vector<int>``
     - Refinement ratio between consecutive levels (length ``n_lev - 1``).
   * - ``time``
     - ``Real``
     - Simulation time stored in the plotfile.
   * - ``n_lev``
     - ``int``
     - Number of levels loaded.
   * - ``var_names``
     - ``Vector<string>``
     - Names of the loaded variables; mirrors the ``var_names`` argument passed to ``read_plotfile``.
       Must match ``mf[0].nComp()`` for ``write_plotfile(outfile, data)`` to work correctly.

.. hint::
    There is an AMReX struct called ``PlotFileData``, but it only holds metadata and does not include the actual data arrays.  The
    ``PlotfileData`` struct defined here is a separate, higher-level construct that includes the data arrays and is specific to this utility library.

MEFData
~~~~~~~

Returned by ``read_mef`` / passed to ``write_mef``.  Represents an
unstructured surface or polyline in Marc's Element Format.

.. list-table::
   :header-rows: 1
   :widths: 25 20 55

   * - Field
     - Type
     - Description
   * - ``title``
     - ``string``
     - Label stored on the first line of the file.
   * - ``var_names``
     - ``vector<string>``
     - Per-component variable names.
   * - ``nodes``
     - ``FArrayBox``
     - Node data: ``Box`` is ``(0:nNodes-1, 0:0, 0:0)`` with ``nComp`` components.
   * - ``connectivity``
     - ``Vector<int>``
     - Flat array of length ``n_elts * nodes_per_elt``; 1-based node indices.
   * - ``n_elts``
     - ``int``
     - Number of elements.
   * - ``nodes_per_elt``
     - ``int``
     - Nodes per element (e.g. 3 for triangles, 2 for line segments).


Functions
---------

read_plotfile
~~~~~~~~~~~~~

**What it does**

Opens an AMReX plotfile, reads the requested variables into one
``MultiFab`` per AMR level, and returns a ``PlotfileData`` struct.
MPI-safe: redistributes boxes when the domain has fewer boxes than
MPI ranks to avoid ``FillVar`` deadlocks.

**Inputs**

.. list-table::
   :header-rows: 1
   :widths: 30 20 50

   * - Parameter
     - Default
     - Description
   * - ``infile``
     - required
     - Path to the plotfile directory.
   * - ``var_names``
     - required
     - | Variables to load; each name must exist in the file (aborts otherwise).
       | This list should only contain unique values; if a name appears multiple 
       | times, it is loaded multiple times.
       | Currently, there is no option to simply load aall variables; you must 
       | explicitly list the ones you want. This is to avoid accidentally loading 
       | large arrays that you don't need.
   * - ``finest_level``
     - ``1000``
     - | Cap on the highest AMR level to read; clamped to the file's actual finest 
       | level.
   * - ``n_grow``
     - ``0``
     - Number of ghost cells to allocate in the returned ``MultiFab``.
   * - ``is_per``
     - ``{}`` (all 0)
     - Periodicity flags, one per spatial dimension.

**Output**

``PlotfileData`` — see structure description above.

**Example**

::

   Vector<std::string> vars = {"density", "velocity_x"};
   Vector<int> is_per(AMREX_SPACEDIM, 1);   // fully periodic domain
   auto data = analysis_util::read_plotfile("plt00100", vars, /*finestLevel=*/2,
                                            /*n_grow=*/1, is_per);
   Print() << "Loaded " << data.n_lev << " level(s), t = " << data.time << "\n";


----

write_plotfile
~~~~~~~~~~~~~~

**What it does**

Writes a multi-level set of ``MultiFab`` objects to an AMReX plotfile.
MPI-safe: adds ``Barrier`` calls before and after the collective write,
and redistributes boxes when the domain has fewer boxes than MPI ranks.

**Inputs**

.. list-table::
   :header-rows: 1
   :widths: 30 20 50

   * - Parameter
     - Default
     - Description
   * - ``outfile``
     - required
     - Output directory name (created by AMReX).
   * - ``mf``
     - required
     - ``Vector<MultiFab>`` to write, one per level.
   * - ``var_names``
     - required
     - Variable name for each component; length must equal ``mf[0].nComp()``.
   * - ``geoms``
     - required
     - Domain geometry for each level.
   * - ``time``
     - ``0.0``
     - Simulation time to embed in the plotfile header.
   * - ``ref_ratios``
     - ``{}`` (all 2)
     - Refinement ratio between consecutive levels (length ``n_lev - 1``).

.. hint::
    A convenience overload is available that takes a single ``PlotfileData`` struct as input: ``write_plotfile(outfile, data)`` expands to
    ``write_plotfile(outfile, data.mf, data.var_names, data.geoms, data.time, data.ref_ratios)``.

**Output**

None (writes to disk).

**Example**

::

   const std::string outfile = "plt_result";
   analysis_util::write_plotfile(outfile, data.mf, {"density", "velocity_x"},
                                 data.geoms, data.time, data.ref_ratios);


----

read_mef
~~~~~~~~

**What it does**

Reads a binary MEF file from disk and returns a ``MEFData`` struct
containing the title, variable names, node data, and element connectivity.

**Inputs**

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Parameter
     - Description
   * - ``infile``
     - Path to the MEF file.

**Output**

``MEFData`` — see structure description above.

**Example**

::

   auto mef = analysis_util::read_mef("surface.mef");
   Print() << "Title: " << mef.title << "\n";
   Print() << mef.n_elts << " elements, " << mef.nodes_per_elt << " nodes/elt\n";


----

write_mef
~~~~~~~~~

**What it does**

Writes a ``MEFData`` struct to disk in binary MEF format.
Only the IO rank writes; all ranks synchronize via a ``Barrier``
before returning so that a subsequent ``read_mef`` on any rank sees
the complete file.

**Inputs**

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Parameter
     - Description
   * - ``outfile``
     - Output file path.
   * - ``data``
     - ``MEFData`` to write.

**Output**

None (writes to disk).

**Example**

::

   analysis_util::MEFData mef;
   mef.title        = "my surface";
   mef.var_names    = {"x", "y", "z", "temperature"};
   mef.n_elts       = nTri;
   mef.nodes_per_elt = 3;
   // ... fill mef.nodes and mef.connectivity ...
   analysis_util::write_mef("surface_out.mef", mef);


----

get_covered_mf
~~~~~~~~~~~~~~

**What it does**

Builds an integer mask ``MultiFab`` at each AMR level.
A cell is marked ``1`` (uncovered) if it is not covered by a finer
level, and ``0`` (covered) otherwise.
Used internally by ``integrate`` to avoid double-counting.

**Inputs**

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Parameter
     - Description
   * - ``mf``
     - ``Vector<MultiFab>`` — one MultiFab per AMR level; grid structure is taken from its ``boxArray`` and ``DistributionMap``.
   * - ``ref_ratios``
     - ``Vector<int>`` — refinement ratio between consecutive levels (length ``n_lev - 1``); pass empty for a single-level domain.

.. hint::
   A ``PlotfileData`` convenience overload is also available: ``get_covered_mf(data)``
   expands to ``get_covered_mf(data.mf, data.ref_ratios)``.

**Output**

``Vector<iMultiFab>`` — one mask MultiFab per level.

**Example**

::

   Vector<int> no_ratios;   // single-level domain
   auto mask = analysis_util::get_covered_mf(mf, no_ratios);
   // or using PlotfileData directly:
   auto mask = analysis_util::get_covered_mf(data);
   // mask[lev] == 1 where lev is the finest data at that location


----

integrate
~~~~~~~~~

**What it does**

Computes the volume integral of one or more components across all AMR
levels, correctly excluding fine-covered cells.
Partial integrals along selected axes are supported: e.g. integrating
over x and z leaves a 1-D array indexed by y.
Result is reduced across all MPI ranks before returning.

**Inputs**

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Parameter
     - Description
   * - ``a_mf``
     - ``Vector<MultiFab>`` — data to integrate (one per level).
   * - ``geoms``
     - ``Vector<Geometry>`` — domain geometry per level (used to compute cell volumes).
   * - ``ref_ratios``
     - ``Vector<int>`` — refinement ratio between consecutive levels; pass empty for single-level.
   * - ``scomp``
     - First component index to integrate.
   * - ``ncomp``
     - Number of components to integrate, starting at ``scomp``.
   * - ``axes_to_integrate``
     - ``Vector<int>`` of axis indices to sum over (0=x, 1=y, 2=z). Pass all three for a scalar result.

.. hint::
   A ``PlotfileData`` convenience overload is also available: ``integrate(data, scomp, ncomp, axes)``
   expands to ``integrate(data.mf, data.geoms, data.ref_ratios, scomp, ncomp, axes)``.

**Output**

``unique_ptr<Gpu::ManagedVector<Real>>`` of length
``ncomp * (Nx if x not integrated, else 1) * (Ny if y not integrated, else 1) * (Nz if z not integrated, else 1)``,
where Nx/Ny/Nz are the cell counts of the finest level.

**Example — scalar volume integral**

::

   // Using PlotfileData overload (most common):
   auto result = analysis_util::integrate(data, /*scomp=*/0, /*ncomp=*/1,
                                          {0, 1, 2});
   Print() << "Integral = " << (*result)[0] << "\n";

   // Or with explicit parameters:
   auto result = analysis_util::integrate(data.mf, data.geoms, data.ref_ratios,
                                          0, 1, {0, 1, 2});

**Example — line-average along y**

::

   // Integrate over x and z; result is an array indexed by y
   auto profile = analysis_util::integrate(data, 0, 1, {0, 2});
   // (*profile)[j] is the x-z integral at y-index j


----

gradient
~~~~~~~~

**What it does**

Computes the cell-centered gradient of ``ncomp`` scalar fields using
AMReX's ``MLMG`` / ``MLPoisson`` solver with 4th-order accuracy.
Boundary conditions are set automatically from the periodicity and
an optional ``sym_dir`` ParmParse flag for symmetric (reflect-odd)
directions.

**Inputs**

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Parameter
     - Description
   * - ``a_mf``
     - ``Vector<MultiFab>`` — source data; **must have at least 1 ghost cell**.
   * - ``geoms``
     - ``Vector<Geometry>`` — domain geometry per level.
   * - ``scomp``
     - First component to differentiate.
   * - ``ncomp``
     - Number of scalar fields to differentiate.
   * - ``grad_mf``
     - ``Vector<MultiFab>`` — output; must have ``AMREX_SPACEDIM * ncomp`` components.

.. hint::
   A ``PlotfileData`` convenience overload is also available: ``gradient(data, scomp, ncomp, grad_mf)``
   expands to ``gradient(data.mf, data.geoms, scomp, ncomp, grad_mf)``.

**Output**

None (fills ``grad_mf`` in-place).
Components are ordered as ``(df0/dx, df0/dy, df0/dz, df1/dx, ...)``.

**Example**

::

   // Read with n_grow=1 so the solver has ghost cells
   auto data = analysis_util::read_plotfile(infile, {"temperature"},
                                            finestLevel, /*n_grow=*/1, is_per);

   Vector<MultiFab> grad_mf(data.n_lev);
   for (int lev = 0; lev < data.n_lev; ++lev)
     grad_mf[lev].define(data.mf[lev].boxArray(),
                         data.mf[lev].DistributionMap(),
                         AMREX_SPACEDIM, 0);
   // PlotfileData overload:
   analysis_util::gradient(data, /*scomp=*/0, /*ncomp=*/1, grad_mf);
   // grad_mf[lev] now holds (dT/dx, dT/dy [, dT/dz])


----

get_file_root
~~~~~~~~~~~~~

**What it does**

Strips the directory prefix from a file path, returning only the
final component.

**Input**

Path string (e.g. ``"runs/case1/plt00100"``).

**Output**

Root filename string (e.g. ``"plt00100"``).

**Example**

::

   std::string root = analysis_util::get_file_root("runs/case1/plt00100");
   // root == "plt00100"
   const std::string outfile = root + "_processed";


----

find_var_index
~~~~~~~~~~~~~~

**What it does**

Searches a list of variable names for a given name and returns its
index.  Aborts with a helpful message if the name is not found
(unless ``abort_if_not_found = false``, in which case it returns ``-1``).

**Inputs**

.. list-table::
   :header-rows: 1
   :widths: 35 15 50

   * - Parameter
     - Default
     - Description
   * - ``var_names``
     - required
     - List of names to search (e.g. ``data.var_names``).
   * - ``name``
     - required
     - Variable name to find.
   * - ``abort_if_not_found``
     - ``true``
     - If ``true``, aborts when the name is absent; if ``false``, returns ``-1``.

**Output**

``int`` — 0-based component index, or ``-1`` if not found and
``abort_if_not_found = false``.

**Example**

::

   int iTemp = analysis_util::find_var_index(data.var_names, "temperature");
   // Use iTemp as scomp in integrate() or as an array index into data.mf


----

parse_title / parse_var_names
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**What they do**

Low-level stream parsers for reading the header section of MEF-style
binary files.

- ``parse_title(is)`` — reads and returns one line from ``is`` as the title string.
- ``parse_var_names(is)`` — reads the next line from ``is`` and splits it on
  spaces and commas, returning the tokens as a ``vector<string>``.

These are called internally by ``read_mef``; you only need them when
writing a custom file reader that follows the same header convention.

**Example**

::

   std::ifstream ifs("custom.mef");
   std::string title        = analysis_util::parse_title(ifs);
   auto        var_names    = analysis_util::parse_var_names(ifs);


