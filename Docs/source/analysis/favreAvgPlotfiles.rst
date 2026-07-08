favreAvgPlotfiles
==================

Overview
########

Utility to average AMReX plotfiles on the same physical domain, including cases
where the input files have non-matching AMR refinement structures. The tool
computes density-weighted ensemble or time-averaged statistics from multiple
plotfiles. Depending on the selected options, the output can contain mean fields,
Favre-averaged mean fields, and variance-like quantities.

The tool always uses the density field from the input plotfiles and forms
density-weighted moments internally.

Usage
#####

::

   ./favreAvgPlotfiles3d.gnu.ex infiles=<file1 file2 ...> [options]

Example
#######

::

   ./favreAvgPlotfiles3d.gnu.ex infiles="plt29000 plt30000 plt31000" outfile="avg_result"

   ./favreAvgPlotfiles3d.gnu.ex infiles="plt29000 plt30000" do_variance=0 do_divide=1

   ./favreAvgPlotfiles3d.gnu.ex infiles="plt29000 plt30000" variables="x_velocity y_velocity z_velocity temp" outfile="plt_favre_avg"

Tool Options
############

::

   #------------------- IO CONTROL -----------------------------------------------------------
   infiles = plt29000 plt30000 plt31000     # List of AMReX plotfiles to average [REQUIRED]
   outfile = plt_averaged                   # DEF: plt_averaged; Output averaged plotfile name

``infiles`` specifies the list of input AMReX plotfiles to be averaged. All input
files must have the same physical domain geometry. Multiple files are provided as
a space-separated list.

``outfile`` specifies the name of the output plotfile containing the averaged
statistics. If not provided, the default is ``plt_averaged``.

::

   #------------------- VARIABLE SELECTION --------------------------------------------------
   variables = x_velocity y_velocity z_velocity density    # Variables to process (DEF: all)

``variables`` allows selective extraction of variables from the input plotfiles.
If not specified, all variables present in the input files are processed.
Variable names must match exactly the names in the input plotfiles.

The field ``density`` must be present in the input plotfiles, even if it is not
listed in ``variables``. It is used internally to form density-weighted moments.

::

   #------------------- AVERAGING OPTIONS ---------------------------------------------------
   do_average = 1                           # DEF: 1; Output first-moment averages
   do_variance = 1                          # DEF: 1; Output variance-like quantities
   do_divide = 1                            # DEF: 1; Divide density-weighted moments by rho_mean

These flags control which statistics are written to the output plotfile.

``do_average`` enables output of the first-moment fields. When set to 1, the tool
outputs either Favre means or undivided density-weighted means, depending on
``do_divide``.

``do_variance`` enables output of variance-like quantities. At least one of
``do_average`` or ``do_variance`` must be set to 1.

``do_divide`` controls whether the density-weighted moments are divided by the
mean density.

When ``do_divide=1``, the tool outputs Favre-averaged quantities:

.. math::

   \widetilde{\phi}
   =
   \frac{\langle \rho \phi \rangle}{\langle \rho \rangle}

and the Favre variance-like quantity:

.. math::

   \widetilde{\phi''^2}
   =
   \frac{\langle \rho \phi^2 \rangle}{\langle \rho \rangle}
   -
   \left(
   \frac{\langle \rho \phi \rangle}{\langle \rho \rangle}
   \right)^2

When ``do_divide=0``, the tool does not divide the accumulated moments by
``rho_mean``. In this case, the output mean fields are undivided
density-weighted moments:

.. math::

   \langle \rho \phi \rangle

The corresponding variance output is computed from the undivided moments as:

.. math::

   \langle \rho \phi^2 \rangle
   -
   \langle \rho \phi \rangle^2

This quantity is not a standard Reynolds variance and should be interpreted as
an undivided density-weighted diagnostic.

::

   #------------------- AMR CONTROL ---------------------------------------------------------
   output_max_level = 1000                  # DEF: 1000; Maximum refinement level to keep
   output_max_grid_size = 32                # DEF: 32; Maximum grid size in output
   interp_type = 1                          # DEF: 1; 0=piecewise constant, 1=linear interpolation

``output_max_level`` specifies the finest AMR level to include in the output.
The value is zero-indexed in the input option. Internally, the tool adds one to
this value to account for the base level. The default value is large enough to
include all available levels in most cases.

``output_max_grid_size`` sets the maximum grid size for output boxes when the
input files have different grid structures. This option is only applied on levels
where the input BoxArrays are not identical. Smaller values increase the number
of boxes but may improve load balancing.

``interp_type`` controls the interpolation scheme used by
``PltFileManager::fillPatchFromPlt`` when data are filled onto the combined grid.
Use 0 for piecewise constant interpolation or 1 for linear interpolation.

Output Variables
################

The output variable names depend on the selected averaging options.

The first output component is always:

- ``rho_mean`` — mean density :math:`\langle \rho \rangle`

**With do_divide=1 (Favre-averaged, default):**

When Favre averaging is enabled, the output contains:

- ``rho_mean`` — mean density :math:`\langle \rho \rangle`
- ``<variable>_favre_mean`` — Favre-averaged mean, if ``do_average=1``:

  .. math::

     \widetilde{\phi}
     =
     \frac{\langle \rho \phi \rangle}{\langle \rho \rangle}

- ``<variable>_favre_variance`` — Favre variance-like quantity, if ``do_variance=1``:

  .. math::

     \widetilde{\phi''^2}
     =
     \frac{\langle \rho \phi^2 \rangle}{\langle \rho \rangle}
     -
     \widetilde{\phi}^2

**With do_divide=0 (undivided density-weighted output):**

When division by ``rho_mean`` is disabled, the output contains:

- ``rho_mean`` — mean density :math:`\langle \rho \rangle`
- ``rho_<variable>_mean`` — undivided density-weighted mean, if ``do_average=1``:

  .. math::

     \langle \rho \phi \rangle

- ``rho_<variable>_variance`` — variance-like quantity from undivided moments, if ``do_variance=1``:

  .. math::

     \langle \rho \phi^2 \rangle
     -
     \langle \rho \phi \rangle^2

This ``do_divide=0`` variance output is not a standard Reynolds variance. It is
computed directly from the internally accumulated undivided density-weighted
moments.

**Example output with default settings**

For:

::

   do_average = 1
   do_variance = 1
   do_divide = 1

and selected variables:

::

   variables = x_velocity y_velocity z_velocity density rhoh temp

the output variable list is:

::

   0   rho_mean
   1   x_velocity_favre_mean
   2   y_velocity_favre_mean
   3   z_velocity_favre_mean
   4   density_favre_mean
   5   rhoh_favre_mean
   6   temp_favre_mean
   7   x_velocity_favre_variance
   8   y_velocity_favre_variance
   9   z_velocity_favre_variance
   10  density_favre_variance
   11  rhoh_favre_variance
   12  temp_favre_variance

If ``density`` is included in ``variables``, the field
``density_favre_mean`` corresponds to:

.. math::

   \frac{\langle \rho^2 \rangle}{\langle \rho \rangle}

and should not be confused with ``rho_mean``.

Constraints and Notes
#####################

- All input plotfiles must have the same physical domain geometry and extent.
- The input files may have non-matching AMR refinement structures. The tool
  constructs a combined grid and fills data from each plotfile onto that grid.
- All input files must contain a variable named ``density``.
- If ``variables`` is not specified, all variables in the input files are
  processed.
- If ``variables`` is specified, only the selected variables are processed, but
  ``density`` is still required internally.
- At least one of ``do_average`` or ``do_variance`` must be set to 1.
- ``rho_mean`` is always written as the first output component.
- For compressible turbulent flows, ``do_divide=1`` gives Favre-averaged
  quantities and is typically the physically relevant option.
- With ``do_divide=0``, the tool outputs undivided density-weighted moments and
  corresponding variance-like diagnostics.
- Output file size depends on the number of selected variables, selected
  statistics, and AMR refinement levels. Use ``variables`` and
  ``output_max_level`` to reduce output size when needed.
