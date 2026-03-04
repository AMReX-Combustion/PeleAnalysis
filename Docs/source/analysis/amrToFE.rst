.. highlight:: bash


amrToFE
*******

.. warning::

   This documentation is deprecated or the tool erroneous. 

``amrToFE`` reads an AMReX plotfile and produces an unstructured finite-element style
representation consisting of:

* **Nodes**: either cell-corner nodes (default) or cell-center nodes (``connect_cc=1``)
* **Elements**: quadrilaterals in 2D (4 nodes) or bricks/hexahedra in 3D (8 nodes)
* **State variables**: a selectable subset of plotfile components written per node

The primary use case is exporting AMR data to tools that expect a FE connectivity list
(e.g. Tecplot FE zones).

Output formats
==============

``outType=tec`` (default)
  Writes a Tecplot FE zone.

  * **ASCII** (``.dat``): always available.
  * **Binary** (``.plt``): only if compiled with ``USE_TEC_BIN_IO`` and linked against TecIO.
    Enable via ``doBin=1``.

``outType=flt``
  Writes a compact binary file (``.flt``) containing header lines (text) followed by
  binary data blocks:

  1. line 1: ``<infile> time = <time>``
  2. line 2: variable list (same string that would be written to Tecplot)
  3. line 3: ``<nElts> <nodesPerElt>`` (``nodesPerElt`` is 4 in 2D, 8 in 3D)
  4. node data block: AMReX ``FArrayBox::writeOn`` for the node/state array
  5. connectivity block: raw ``int`` array of length ``nElts*nodesPerElt``

Usage
=====

::

  ./amrToFE*.ex infile=plt00010 [options]

The executable also accepts a ParmParse input file as first argument:

::

  ./amrToFE*.ex inputs_amrToFE.inp

Options
=======

Input / output
--------------

``infile`` (required)
  Path to the AMReX plotfile.

``outfile`` (optional)
  Output file name. Default depends on ``outType``:

  * ``outType=tec``  -> ``<infile>.dat`` (or ``<infile>.plt`` if ``doBin=1`` and supported)
  * ``outType=flt``  -> ``<infile>.flt``

``outType`` (optional)
  ``tec`` (default) or ``flt``.

``doBin`` (optional, Tecplot only)
  If ``USE_TEC_BIN_IO`` is available at compile time, set ``doBin=1`` to write Tecplot
  binary ``.plt`` instead of ASCII ``.dat``.

``verbose`` (optional)
  If ``1``, enables verbose AMReX ``AmrData`` output.

Variable selection
------------------

By default, all plotfile components are exported. You can restrict the variables in
two ways:

``comps`` (optional)
  Explicit integer list of component indices to export, e.g. ``comps="0 3 7"``.
  If provided, this overrides ``sComp`` and ``nComp``.

``sComp`` (optional)
  Start component index (default ``0``).

``nComp`` (optional)
  Number of components starting at ``sComp`` (default: all available components).

Spatial subsetting / AMR control
--------------------------------

``box`` (optional)
  Sub-box selection in *index space* at level 0.
  Provide ``2*AMREX_SPACEDIM`` integers: lower-left and upper-right corners.

  *2D example:* ``box="iLo jLo  iHi jHi"``

  The selected region is intersected with the level-0 problem domain.

``finestLevel`` (optional)
  Finest AMR level to include (default: the plotfile's finest level).
  Levels above this are ignored.

``nGrowPer`` (optional)
  If ``>0``, extends the level-0 selection by ``nGrowPer`` cells *only at periodic
  boundaries* (useful to include periodic wrap-around data near domain edges).
  Default: ``0``.

  When using ``nGrowPer>0`` you must provide the standard AMReX geometry parameters,
  for example:

  * ``geometry.coord_sys`` (0 Cartesian, 1 RZ)
  * ``geometry.is_periodic`` (per-axis)
  * ``geometry.prob_lo``
  * ``geometry.prob_hi``

Connectivity mode
-----------------

``connect_cc`` (optional, default ``1``)
  Controls how nodes and values are constructed:

  * ``connect_cc=1``: **Flattened** structure by connecting **cell centers**.
  * ``connect_cc=0``: Create nodes at **all cell corners** and copy the cell-centered
    values out to the nodes.

Help
----

``help``
  Print the built-in usage text and exit.

Examples
========

Tecplot ASCII (default)
-----------------------

::

  ./amrToFE2d.gnu.MPI.ex infile=plt00050 outType=tec outfile=plt00050.dat

Tecplot binary (requires USE_TEC_BIN_IO)
----------------------------------------

::

  ./amrToFE3d.gnu.MPI.ex infile=plt00050 outType=tec doBin=1  # -> plt00050.plt

Select a variable subset and a sub-box
--------------------------------------

::

  ./amrToFE2d.gnu.MPI.ex infile=plt00050 comps="0 2 5" box="0 0  255 255"

Write FLT
---------

::

  ./amrToFE2d.gnu.MPI.ex infile=plt00050 outType=flt outfile=plt00050.flt
