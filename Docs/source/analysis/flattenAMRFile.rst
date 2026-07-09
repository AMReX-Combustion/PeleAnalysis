flattenAMRFile
==============

Description
-----------

``flattenAMRFile`` converts an AMReX multi-level (AMR) plotfile into a
single-level uniform plotfile. Data from the selected AMR level are
interpolated onto a uniform grid covering the full computational domain.
The resulting plotfile contains all variables on one refinement level.

Usage
-----

.. code-block:: bash

   flattenAMRFile infile=<pltfile> [options]


Input Parameters
----------------

Example usage:

.. code-block:: none

   infile = plt00000
   output_file = plt00000_flatten
   output_level = 0
   output_max_grid_size = 64


Parameters
----------

``infile``
   Input AMReX plotfile.

``output_file``
   Name of the generated single-level plotfile.
   Default: ``<infile>_flatten``.

``output_level``
   AMR level taken from the input plotfile and written to the
   flattened output file.
   Default: ``0``.

``output_max_grid_size``
   Maximum grid size used when generating the uniform output grid.
   Default: ``64``.

``verbose``
   Enable additional runtime output.

``n_files``
   Maximum number of binary files used to write the output plotfile data
   (AMReX ``VisMF::SetNOutFiles``). Lower this to reduce the number of files
   created for large parallel post-processing runs. AMReX clamps the value to
   the number of MPI ranks, so a serial run always writes a single data file.
   Default: the AMReX default.


Output
------

A standard AMReX single-level plotfile containing all variables
interpolated onto a uniform grid.


Notes
-----

Flattening removes the AMR hierarchy but preserves the spatial
dimension and physical fields of the original dataset.
