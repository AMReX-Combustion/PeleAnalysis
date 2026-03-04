qCriterion
==========

Description
-----------

``qCriterion`` computes the Q-criterion from AMReX plotfiles based on the
velocity field. The tool reads a plotfile, evaluates velocity gradients,
computes ``Q`` and a normalized variant ``Q_norm``, and writes a new
plotfile containing these derived fields. Optional auxiliary variables
can be copied from the input plotfile to the output without modification.

Usage
-----

.. code-block:: bash

   qCriterion qCriterion.inp


Input File
----------

The program is controlled via an input file with the following structure:

.. code-block:: none

   infile = plt000000
   finestLevel = 2
   is_per = 0 1 0
   sym_dir = 0 0 0
   outfile = plt000000_qCriterion

   Aux_Variables = Y_H2


Parameters
~~~~~~~~~~

``infile``
   Path to the AMReX plotfile.

``finestLevel``
   Finest AMR level up to which the Q-criterion is computed.
   Default: plotfile finest level.

``is_per``
   Periodicity flags in each spatial direction (0: non-periodic, 1: periodic).
   Default: ``0 0 0``.

``sym_dir``
   Symmetry flags in each spatial direction (0: non-symmetric, 1: symmetric).
   Default: ``0 0 0``.

``outfile``
   Output plotfile name.
   Default: ``<infile>_qCriterion``.

``Aux_Variables``
   Names of variables copied unchanged from input to output plotfile.
   Default: none.


Output
------

A new AMReX plotfile containing:

- copied input variables (velocity components and optional ``Aux_Variables``)
- velocity-gradient components (e.g. ``x_velocity_gx``, ``x_velocity_gy``, ...)
- ``Q`` and ``Q_norm``


Typical Applications
--------------------

- Vortex identification via Q-criterion
- Flow-structure visualization
- Turbulence diagnostics and post-processing

Notes
-----

The input plotfile must contain the velocity components named
``x_velocity``, ``y_velocity``, and ``z_velocity``. The tool is intended
for 3D datasets (compiled with ``DIM=3``).
