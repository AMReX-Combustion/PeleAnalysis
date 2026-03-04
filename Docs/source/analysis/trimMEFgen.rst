.. highlight:: bash

trimMEFgen
**********

Trim and post-process 3D isosurface meshes by removing nodes and triangular elements based on field-value thresholds or radial distance criteria. The tool reads a binary isosurface file, applies any requested trimming operations, and writes the processed surface to a new file. Surface area is reported before and after trimming.

Usage: ::

   ./trimMEFgen2d.gnu.ex infile=<s> outfile=<s> [options]

Example: ::

   ./trimMEFgen3d.gnu.ex ./InputSamples/isoTrim.inp

Tool Options
#############

::

   #------------------- IO CONTROL -----------------------------------------------------------
   infile  = surface.dat                      # Input binary isosurface file to be processed
   outfile = surface_trimmed.dat              # Output file for the processed isosurface

`infile` specifies the binary isosurface file to read. This file is expected to contain a header with a label and variable names, followed by node data and triangular element connectivity. `outfile` specifies where the trimmed surface will be written in the same binary format.

::

   #------------------- Field-Value Trimming -------------------------------------------------
   comps   = 3 4                              # Component indices to use for trimming
   signs   = gt lt                            # Comparison operators for each component (lt,le,gt,ge,eq)
   vals    = 0.5 1.0                          # Threshold values; nodes satisfying condition are removed

`comps` takes a space-separated list of zero-based component indices identifying which field variables to use for trimming. For each component, a corresponding entry must be provided in `signs` and `vals`. `signs` specifies the comparison operator applied to each component, and must be one of ``lt`` (less than), ``le`` (less than or equal), ``gt`` (greater than), ``ge`` (greater than or equal), or ``eq`` (equal). `vals` provides the corresponding threshold value. A node is removed if **any** of the specified conditions evaluates to true. Any triangular element with at least one removed node is also discarded.

::

   #------------------- Radial Trimming ------------------------------------------------------
   RXY      = 0.05                            # Radial threshold: trim based on r = sqrt(x^2 + y^2)
   sign_RXY = lt                              # Comparison operator for radial trim (lt,le,gt,ge,eq)

`RXY` enables trimming based on the cylindrical radius :math:`r = \sqrt{x^2 + y^2}` of each node. This is useful for axisymmetric or cylindrical geometries where a clean inner or outer radial boundary is desired. `sign_RXY` specifies the comparison operator and accepts the same values as `signs`. For example, ``sign_RXY = lt`` with ``RXY = 0.05`` removes all nodes closer than 0.05 to the z-axis. Field-value trimming (via `comps`) and radial trimming can be applied simultaneously; both are applied independently before the surface area is recomputed.

::

   #------------------- Output Control -------------------------------------------------------
   remComps = 5 6                             # Component indices to drop from the output file

`remComps` takes a space-separated list of zero-based component indices to remove from the output data. This is useful for stripping intermediate or auxiliary variables (such as the isosurface construction field itself) that are not needed in downstream analysis. The remaining components are written to `outfile` in their original order.

::

   #------------------- Diagnostics ----------------------------------------------------------
   do_area_stats = true                       # Print min/max triangle areas across the surface

`do_area_stats` enables reporting of per-triangle area statistics. When set to ``true``, the tool computes the area of every triangle in the trimmed surface and prints the minimum and maximum values. This is useful for diagnosing mesh quality — a large spread between minimum and maximum triangle areas may indicate poorly resolved regions of the surface that could affect downstream calculations.
