.. highlight:: bash

******************************************
generateTestPlt
******************************************
Generate a synthetic multilevel AMReX plot file populated with one or more analytically defined fields. Fields can be simple constants, planar step or smooth transitions, circular/spherical regions, annular rings or spherical shells, or sinusoidal patterns. Refinement regions can be specified geometrically (bounding box) or by field-value criteria. This tool is primarily used to create test and validation data for other tools in the PeleAnalysis suite without requiring a full simulation.

Usage: ::

   ./generateTestPlt.gnu.ex [OPTIONS]

Example: ::

   ./generateTestPlt.gnu.ex ./InputSamples/generateTestPlt.inp


Tool Options
#############
::

   #------------------- Domain -----------------------------------------------------------
   geometry.prob_lo = 0.0 0.0 0.0        # Physical lower bound of the domain (required)
   geometry.prob_hi = 1.0 1.0 1.0        # Physical upper bound of the domain (required)
   geometry.is_periodic = 1 1 1          # Periodicity flags per dimension (required)
   geometry.coord_sys = 0                # DEF: 0; Coordinate system (0=Cartesian, 1=cylindrical/RZ — 2D builds only)

These parameters follow the standard AMReX geometry convention and are all required except ``coord_sys``. ``coord_sys = 1`` (cylindrical/RZ) requires a 2D build (``DIM=2``) and produces a Header with ``spacedim=2`` matching PeleLMeX 2D RZ output.
::

   #------------------- Grid ------------------------------------------------------------
   amr.n_cell = 64 64 64                 # Number of cells per dimension on level 0 (required)
   amr.max_grid_size = 32                # DEF: 32; Maximum box size for domain decomposition
   amr.ngrow = 0                         # DEF: 0; Number of ghost cells

``amr.n_cell`` sets the level-0 resolution and is required. ``max_grid_size`` controls the parallel decomposition and does not affect the output.
::

   #------------------- AMR / Refinement ------------------------------------------------
   amr.max_level = 2                     # DEF: 0; Maximum refinement level (0 = single-level)
   amr.ref_ratio = 2 2                   # DEF: 2; Refinement ratio, one value per level transition
   amr.grid_eff = 0.7                    # DEF: 0.7; Clustering efficiency threshold
   amr.blocking_factor = 8               # DEF: 8; Minimum box size during grid generation
   amr.n_error_buf = 2 2                 # DEF: 1; Buffer cells per level around tagged regions

``amr.max_level`` enables multilevel output. When greater than 0, ``amr.refinement_indicators`` must list one or more criteria that determine where finer grids are placed. ``amr.ref_ratio`` accepts one integer per level transition; if fewer values are given the last value is repeated. ``amr.n_error_buf`` expands each tagged region by that many cells (in coarse index space) before clustering into boxes.

Refinement Indicators
######################
::

   amr.refinement_indicators = NAME [NAME ...]   # Space-separated list of indicator names

Each name in ``amr.refinement_indicators`` introduces a refinement criterion whose parameters are read under the ``amr.<NAME>.*`` namespace. A single plotfile run may combine multiple indicators of any type.

Common option for all indicator types::

   amr.<name>.max_level = N              # DEF: amr.max_level; Max refinement level for this criterion

Available criterion types
--------------------------

.. list-table::
   :widths: 30 70
   :header-rows: 1

   * - Criterion
     - Description
   * - ``in_box_lo`` / ``in_box_hi``
     - Refine all cells whose centres fall within the specified physical bounding box
   * - ``value_greater`` + ``field_name``
     - Refine cells where the named field evaluates above the threshold
   * - ``value_less`` + ``field_name``
     - Refine cells where the named field evaluates below the threshold
   * - ``adjacent_difference_greater`` + ``field_name``
     - Refine cells where the maximum absolute difference to any immediate neighbour exceeds the threshold

**Box-based criterion** ::

   amr.<name>.in_box_lo = 0.2 0.2 0.2   # Physical lower corner of the refinement region
   amr.<name>.in_box_hi = 0.8 0.8 0.8   # Physical upper corner of the refinement region

**Value-greater criterion** ::

   amr.<name>.value_greater = 5.0        # Refine where field > 5.0
   amr.<name>.field_name    = circle_step

**Value-less criterion** ::

   amr.<name>.value_less = 0.1           # Refine where field < 0.1
   amr.<name>.field_name = plane_step

**Adjacent-difference criterion** ::

   amr.<name>.adjacent_difference_greater = 0.4   # Refine at sharp interfaces
   amr.<name>.field_name                  = plane_step

.. note::

   Value-based criteria evaluate the analytic field expression at cell centres. They do not depend on data stored in the plotfile and can therefore be used at any refinement level without first filling a coarser level.

::

   #------------------- Output ----------------------------------------------------------
   plotfile_name = pltTestFile            # DEF: pltTestFile; Name of the output plot file
   n_files = 64                           # DEF: AMReX default; cap on the number of plotfile data files

The output is a standard AMReX multilevel plot file (one level when ``amr.max_level = 0``) readable by any tool in the PeleAnalysis suite or visualised with VisIt or ParaView. ``n_files`` caps the number of binary files used for the plotfile data (AMReX ``VisMF::SetNOutFiles``); AMReX clamps it to the number of MPI ranks, so a serial run always writes a single data file.
::

   #------------------- Fields ----------------------------------------------------------
   field.names = myField1 myField2 ...   # Space-separated list of field names (required)

``field.names`` defines the list of fields to generate. At least one field must be specified or the tool will abort. Each name becomes both the variable name in the plot file and the ParmParse prefix used to configure that field. For a field named ``myField``, all its parameters are read under the ``myField.*`` namespace.

Field Configuration
####################
Each field listed in ``field.names`` must have a ``type`` parameter and a set of type-specific parameters. All parameters within a field use the field name as a ParmParse prefix (e.g. ``myField.type``, ``myField.radius``, etc.).

Common parameters
-----------------
::

   myField.type = <type_string>           # Field type (required, see table below)
   myField.output_name = Y(OH)            # DEF: <field_name>; Variable name written to the plotfile header
   myField.value_inside  = 1.0            # DEF: 1.0; Value inside the region
   myField.value_outside = 0.0            # DEF: 0.0; Value outside the region

.. note::

   ``output_name`` lets you embed characters that ParmParse treats as key separators (e.g. ``/``) in the
   plotfile variable name without affecting the ParmParse prefix used to configure the field.  For example,
   setting ``myField.output_name = Y(OH)`` writes the variable as ``Y(OH)`` in the plotfile header while
   ``myField`` remains the valid ParmParse prefix.  This is the recommended way to create test fixtures for
   tools that handle slash-containing species names such as ``jpdf``.

Available field types
---------------------

.. list-table::
   :widths: 30 70
   :header-rows: 1

   * - Type string
     - Description
   * - ``constant``
     - Uniform value everywhere
   * - ``plane_step``
     - Step function at a single plane
   * - ``plane_smooth``
     - Smoothly blended transition at a single plane
   * - ``double_plane_step``
     - Step function between two parallel planes (slab)
   * - ``double_plane_smooth``
     - Smoothly blended slab between two parallel planes
   * - ``circle_step`` / ``sphere_step``
     - Step function inside a circle or sphere
   * - ``circle_smooth`` / ``sphere_smooth``
     - Smoothly blended circle or sphere
   * - ``ring_step`` / ``spherical_shell_step``
     - Step function in an annular ring or spherical shell
   * - ``ring_smooth`` / ``spherical_shell_smooth``
     - Smoothly blended annular ring or spherical shell
   * - ``cylinder_step``
     - Step function based on radial distance from an infinite axis line (3D Cartesian)
   * - ``cylinder_smooth``
     - Smoothly blended transition based on radial distance from an infinite axis line (3D Cartesian)
   * - ``sine``
     - Sinusoidal field: offset + amplitude * sin(2π f_x x) * cos(2π f_y y) * sin(2π f_z z)

.. note::

   ``circle_*`` and ``sphere_*`` are identical in behaviour — both compute a radial distance from a centre point using all spatial dimensions. The same applies to ``ring_*`` and ``spherical_shell_*``.

.. note::

   ``cylinder_*`` measures the radial distance from an *axis line* (not a point), making it suitable for jet or pipe geometries on 3D Cartesian grids. The axis direction is set with ``axis``; the ``center`` parameter gives any point on that axis line.

Smooth transitions use a *smootherstep* function (5th-order Hermite) centred on the specified position, over a width of ``smooth_width``.

Type-specific parameters
------------------------

**constant** ::

   myField.value = 1.0                    # DEF: 1.0; Uniform field value

**plane_step** / **plane_smooth** ::

   myField.axis         = 0               # DEF: 0; Normal axis (0=x, 1=y, 2=z)
   myField.position     = 0.5             # DEF: 0.5; Plane position in physical coordinates
   myField.value_left   = 1.0             # DEF: 1.0; Value for coord < position
   myField.value_right  = 0.0             # DEF: 0.0; Value for coord > position
   myField.smooth_width = 0.1             # DEF: 0.1; Transition width (plane_smooth only)

**double_plane_step** / **double_plane_smooth** ::

   myField.axis              = 0          # DEF: 0; Normal axis (0=x, 1=y, 2=z)
   myField.position          = 0.3        # DEF: 0.5; Position of first plane
   myField.position_second   = 0.7        # DEF: 0.6; Position of second plane
   myField.value_inside      = 1.0        # DEF: 1.0; Value between the two planes
   myField.value_outside     = 0.0        # DEF: 0.0; Value outside the two planes
   myField.smooth_width      = 0.05       # DEF: 0.1; Transition width at each plane (smooth only)

**circle_step** / **sphere_step** / **circle_smooth** / **sphere_smooth** ::

   myField.center       = 0.5 0.5 0.5    # DEF: 0.5 0.5 0.5; Centre of the circle/sphere
   myField.radius       = 0.25           # DEF: 0.25; Radius
   myField.value_inside = 1.0            # DEF: 1.0; Value inside the radius
   myField.value_outside = 0.0           # DEF: 0.0; Value outside the radius
   myField.smooth_width = 0.05           # DEF: 0.1; Transition width (smooth only)

**ring_step** / **spherical_shell_step** / **ring_smooth** / **spherical_shell_smooth** ::

   myField.center        = 0.5 0.5 0.5   # DEF: 0.5 0.5 0.5; Centre point
   myField.radius_inner  = 0.2           # DEF: 0.25; Inner radius
   myField.radius_outer  = 0.4           # DEF: 0.25; Outer radius
   myField.value_inside  = 1.0           # DEF: 1.0; Value between inner and outer radii
   myField.value_outside = 0.0           # DEF: 0.0; Value outside the shell
   myField.smooth_width  = 0.05          # DEF: 0.1; Transition width at each radius (smooth only)

**cylinder_step** / **cylinder_smooth** ::

   myField.axis          = 2             # DEF: 2; Cylinder axis direction (0=x, 1=y, 2=z)
   myField.center        = 0.5 0.5 0.5  # DEF: 0.5 0.5 0.5; Any point on the axis line
   myField.radius        = 0.25         # DEF: 0.25; Cylinder radius
   myField.value_inside  = 1.0          # DEF: 1.0; Value inside the cylinder
   myField.value_outside = 0.0          # DEF: 0.0; Value outside the cylinder
   myField.smooth_width  = 0.1          # DEF: 0.1; Transition width (cylinder_smooth only)

**sine** ::

   myField.frequency = 1.0 1.0 1.0       # DEF: 1.0 1.0 1.0; Spatial frequency per dimension
   myField.phase     = 0.0 0.0 0.0       # DEF: 0.0 0.0 0.0; Phase offset per dimension (radians)
   myField.amplitude = 1.0               # DEF: 1.0; Amplitude
   myField.offset    = 0.0               # DEF: 0.0; Additive offset

   # Result: offset + amplitude * sin(2π f_x x + φ_x) * cos(2π f_y y + φ_y) * sin(2π f_z z + φ_z)
