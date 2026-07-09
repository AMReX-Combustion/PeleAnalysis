curvature
=========

Description
-----------

``curvature`` computes curvature-related quantities from AMReX plotfiles
based on a user-defined progress variable. The tool reads an input
plotfile, evaluates geometric and kinematic quantities, and writes a
new plotfile containing the computed fields.

Usage
-----

.. code-block:: bash

   curvature curvature.inp


Input File
----------

Example input file:

.. code-block:: none

   infile = plt00000
   outfile = plt00000_curvature
   finestLevel = 0

   Aux_Variables = 2 3 4
   progressName = Y_H2
   progMin = 0.0
   progMax = 1.0
   useFileMinMax = false
   threshold_prog = false
   threshold_value = 0.01

   do_gaussCurv = true
   do_smooth = true
   smoothing_time = 1.0e-7

   do_strain = true
   getStrainTensor = true
   do_velnormal = true


Parameters
----------

``infile``
   Input AMReX plotfile.

``outfile``
   Output plotfile containing computed quantities.

``finestLevel``
   Finest AMR level to be processed.

``Aux_Variables``
   Names of variables copied unchanged from input to output plotfile.

``progressName``
   Name of the progress variable.

``progMin``, ``progMax``
   Limits applied to the progress variable.

``useFileMinMax``
   Use minimum and maximum values from the plotfile.

``threshold_prog``
   Enable clipping outside the flame front.

``threshold_value``
   Threshold used for clipping.

``do_gaussCurv``
   Compute Gaussian curvature.

``do_smooth``
   Apply smoothing to the progress variable.

``smoothing_time``
   Smoothing time scale.

``do_strain``
   Compute strain-related quantities.

``getStrainTensor``
   Output strain tensor components.

``do_velnormal``
   Compute normal velocity.

``n_files``
   Maximum number of binary files used to write the output plotfile data
   (AMReX ``VisMF::SetNOutFiles``). Lower this to reduce the number of files
   created for large parallel post-processing runs. AMReX clamps the value to
   the number of MPI ranks, so a serial run always writes a single data file.
   Default: the AMReX default.
