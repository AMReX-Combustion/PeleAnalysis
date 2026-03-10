.. highlight:: bash

trimMEFgen
**********

Trim and post-process 2D & 3D isosurface meshes by removing nodes and elements based on field-value thresholds or radial distance criteria. The tool reads a binary isosurface file, applies any requested trimming operations, and writes the processed surface to a new file. Surface area is reported before and after trimming.

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
   comps   = X Y(H2) RXZ                      # Component names to use for trimming. X, Y, Z are also available.
                                              # RXY, RXZ, RYZ, RXYZ refer to the radial/spherical distance.
   signs   = gt lt gt                         # Comparison operators for each component (lt,le,gt,ge,eq)
   vals    = 0.5 1.0 0.1                      # Threshold values; nodes satisfying condition are removed
   rmBound = 0                                # [0, 1], DEF: 1; If 1, elements with at least one nodes satisfying condition are removed.
                                              # If 0, only elements with all nodes satisfying condition are removed.

`comps` takes a space-separated list of component names identifying which field variables to use for trimming. For each component, a corresponding entry must be provided in `signs` and `vals`. `signs` specifies the comparison operator applied to each component, and must be one of ``lt`` (less than), ``le`` (less than or equal), ``gt`` (greater than), ``ge`` (greater than or equal), or ``eq`` (equal). `vals` provides the corresponding threshold value. `comps` can also be the coordinate comps `X`, `Y`, and `Z`. Additionally, `comps` can be `RXY`, `RXZ`, `RYZ`, and `RXYZ`, representing a condition for the radial (`RXY`, `RXZ`, `RYZ`) or spherical (`RXYZ`) coordinate. A node is flaged for removal if **any** of the specified conditions evaluates to true. If `rmBound = 1`, any triangular element with at least one node flagged for removal is also discarded. If `rmBound = 0`, only triangular elements with all corresponding nodes flagged for removal are also discarded. With this option, surfaces can be split lossless. The following bash script can be used to split a surface in `nPiece`.

::
   
   # === Parameters ===
   nPiece=20
   z_min=0.0000
   z_max=0.02880
   dim=Z
   plt_dir=plt43000
   exe=/home/tl277147/Software/Pele/PeleAnalysis/Src/trimMEFgen3d.gnu.MPI.ex

   # === Auto-computed ===
   z_step=$(echo "scale=10; (${z_max} - ${z_min}) / ${nPiece}" | bc)
   base_scale=$(echo "${z_step}" | awk -F'.' '{print length($2)}')
   scale=$(( base_scale + 2 ))

   # === Loop ===
   for (( i=1; i<=nPiece; i++ )); do
      val_low=$(echo  "scale=${scale}; ${z_min} + ${z_step} * ($i - 1)" | bc)
      val_high=$(echo "scale=${scale}; ${z_min} + ${z_step} * $i"       | bc)
      rmBound=$(( (i)  % 2 ))

      echo "Step ${i}/${nPiece} | Z: [${val_low}, ${val_high}] | rmBound=${rmBound}"

      srun ${exe} \
         infile=${plt_dir}_surf.mef \
         outfile=${plt_dir}_surf_p${i}.mef \
         comps = ${dim} ${dim} \
         signs = le gt \
         vals = ${val_low} ${val_high} \
         rmBound = ${rmBound}
   done

::

   #------------------- Output Control -------------------------------------------------------
   remComps = 5 6                             # Component indices to drop from the output file

`remComps` takes a space-separated list of zero-based component indices to remove from the output data. This is useful for stripping intermediate or auxiliary variables (such as the isosurface construction field itself) that are not needed in downstream analysis. The remaining components are written to `outfile` in their original order.

::

   #------------------- Diagnostics ----------------------------------------------------------
   do_area_stats = true                       # Print min/max triangle areas across the surface

`do_area_stats` enables reporting of per-triangle area statistics. When set to ``true``, the tool computes the area of every triangle in the trimmed surface and prints the minimum and maximum values. This is useful for diagnosing mesh quality — a large spread between minimum and maximum triangle areas may indicate poorly resolved regions of the surface that could affect downstream calculations.

.. warning::

   `do_area_stats` only implemented for 3D.

