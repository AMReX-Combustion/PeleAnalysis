scaleMEF
========

Description
-----------

``scaleMEF`` scales selected variables stored on an iso-surface (MEF)
file by user-defined constant factors. The tool reads an existing MEF
surface file, multiplies specified data components at all surface nodes,
and writes a modified MEF file preserving the original geometry and
connectivity.

Optional renaming of scaled variables is supported.

Usage
-----

.. code-block:: bash

   scaleMEF scaleMEF.inp


Input File
----------

The program is controlled via an input file with the following structure:

.. code-block:: none

   infile = surface.mef
   outfile = surface_scaled.mef

   comps = 3 4
   vals  = 10.0 0.5

   newNames = curvature_scaled strain_scaled
   newComps = 3 4


Parameters
~~~~~~~~~~

``infile``
   Input MEF / iso-surface file.

``outfile``
   Output MEF file containing scaled variables.

``comps``
   Indices of components to be scaled.

``vals``
   Scaling factors applied to the corresponding components in
   ``comps``.

Alternative component selection:

``sComp``
   Starting component index.

``nComp``
   Number of consecutive components to scale beginning at
   ``sComp``.

``newNames``
   Optional new variable names written to the output file.

``newComps``
   Component indices whose names are replaced by
   ``newNames``.


Output
------

A modified MEF file containing:

- unchanged surface geometry
- unchanged element connectivity
- scaled node-based variables
- optionally renamed variables


Typical Applications
--------------------

- Non-dimensionalization of surface quantities
- Curvature or strain rescaling
- Unit conversion of iso-surface data
- Preparation of surface data for further analysis


Notes
-----

The number of entries in ``vals`` must match the number of selected
components. Geometry and topology of the iso-surface remain unchanged;
only stored data values are modified.
