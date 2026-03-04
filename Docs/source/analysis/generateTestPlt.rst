.. highlight:: bash
******************************************
generateTestPlt
******************************************
Generate a synthetic single-level AMReX plot file populated with one or more analytically defined fields. Fields can be simple constants, planar step or smooth transitions, circular/spherical regions, annular rings or spherical shells, or sinusoidal patterns. This tool is primarily used to create test and validation data for other tools in the PeleAnalysis suite without requiring a full simulation.

Usage: ::

   ./generateTestPlt.gnu.ex [OPTIONS]

Example: ::

   ./generateTestPlt.gnu.ex ./InputSamples/generateTestPlt.inp

.. note::

   Only single-level (non-AMR) Cartesian grids are supported. Only ``geometry.coord_sys = 0`` (Cartesian) is accepted; other coordinate systems will cause an abort.

Tool Options
#############
::

   #------------------- Domain -----------------------------------------------------------
   geometry.prob_lo = 0.0 0.0 0.0        # Physical lower bound of the domain (required)
   geometry.prob_hi = 1.0 1.0 1.0        # Physical upper bound of the domain (required)
   geometry.is_periodic = 1 1 1          # Periodicity flags per dimension (required)
   geometry.coord_sys = 0                # DEF: 0; Coordinate system (only 0=Cartesian supported)

These parameters follow the standard AMReX geometry convention and are all required except ``coord_sys``.
::

   #------------------- Grid ------------------------------------------------------------
   amr.n_cell = 64 64 64                 # Number of cells per dimension (required)
   amr.max_grid_size = 32                # DEF: 32; Maximum box size for domain decomposition
   amr.ngrow = 0                         # DEF: 0; Number of ghost cells

``amr.n_cell`` sets the resolution of the generated plot file and is required. ``max_grid_size`` controls the internal parallel decomposition and does not affect the output. ``ngrow`` sets the number of ghost cell layers; ghost cells are filled with the same field values as interior cells.
::

   #------------------- Output ----------------------------------------------------------
   plotfile_name = pltTestFile            # DEF: pltTestFile; Name of the output plot file

The output is a standard single-level AMReX plot file that can be read by any tool in the PeleAnalysis suite or visualised with VisIt or ParaView.
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
   myField.value_inside  = 1.0            # DEF: 1.0; Value inside the region
   myField.value_outside = 0.0            # DEF: 0.0; Value outside the region

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
   * - ``sine``
     - Sinusoidal field: offset + amplitude * sin(2π f_x x) * cos(2π f_y y) * sin(2π f_z z)

.. note::

   ``circle_*`` and ``sphere_*`` are identical in behaviour — both compute a radial distance from a centre point using all spatial dimensions. The same applies to ``ring_*`` and ``spherical_shell_*``.

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

**sine** ::

   myField.frequency = 1.0 1.0 1.0       # DEF: 1.0 1.0 1.0; Spatial frequency per dimension
   myField.phase     = 0.0 0.0 0.0       # DEF: 0.0 0.0 0.0; Phase offset per dimension (radians)
   myField.amplitude = 1.0               # DEF: 1.0; Amplitude
   myField.offset    = 0.0               # DEF: 0.0; Additive offset

   # Result: offset + amplitude * sin(2π f_x x + φ_x) * cos(2π f_y y + φ_y) * sin(2π f_z z + φ_z)

