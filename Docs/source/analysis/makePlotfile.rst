makePlotfile
============

Description
-----------

``makePlotfile`` builds a 3D AMReX plotfile from a sequence of 2D FAB files
(e.g. *z-slices*). The tool reads a list of FABs whose filenames are
constructed from a common prefix/suffix and an integer counter, stacks
them along the z-direction, and writes a single-level plotfile containing
the requested variables.

Usage
-----

.. code-block:: bash

   makePlotfile makePlotfile.inp


Input File
----------

The program is controlled via an input file with the following structure:

.. code-block:: none

   infile_head = plane
   infile_tail = .zslice.0.Level_0.fab
   outfile = output
   start = 0
   end = 10
   interval = 1
   probLo = -0.05 -0.05 -0.05
   probHi =  0.05  0.05  0.05
   time = 0.0
   names = x_velocity y_velocity z_velocity


Parameters
~~~~~~~~~~

``infile_head``
   Prefix of the input FAB filenames.

``infile_tail``
   Suffix of the input FAB filenames.

``outfile``
   Output plotfile name.

``start``
   Start index inserted between ``infile_head`` and ``infile_tail``.
   The index is formatted as a 6-digit integer (``%06d``).

``end``
   End index inserted between ``infile_head`` and ``infile_tail``.

``interval``
   Index increment between successive FAB files.

``names``
   Variable names contained in each FAB (and written to the plotfile).

``probLo``
   Physical lower corner of the output domain. Provide either 3 values
   (x y z) or 2 values (x y). If 2 values are provided, the z-extent is
   inferred from the number of planes and a spacing based on the x-spacing.

``probHi``
   Physical upper corner of the output domain (x y z). Must provide 3 values.

``time``
   Physical time assigned to the output plotfile.
   Default: ``0.0``.

``verbose`` (or ``v``)
   Enable verbose output.

``n_files``
   Maximum number of binary files used to write the output plotfile data
   (AMReX ``VisMF::SetNOutFiles``). Lower this to reduce the number of files
   created for large parallel post-processing runs. AMReX clamps the value to
   the number of MPI ranks, so a serial run always writes a single data file.
   Default: the AMReX default.

Output
------

A standard single-level AMReX plotfile containing the stacked 3D fields.

Typical Applications
--------------------

- Reconstruct a 3D dataset from stored 2D slices
- Convert slice-based outputs into a plotfile for downstream AMReX/yt tools
- Post-processing pipelines requiring plotfile input

Notes
-----

All input FABs must have identical 2D dimensions and contain the same number
of components as listed in ``names``. The output uses a slab decomposition
(one z-plane per box).

