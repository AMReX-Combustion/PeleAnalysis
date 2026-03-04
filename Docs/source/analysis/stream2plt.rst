stream2plt
==========

Description
-----------

``stream2plt`` extracts streamline data from AMReX plotfiles and writes
the resulting streamlines to a new plotfile. Streamlines are generated
from selected vector components and may optionally be filtered using
user-defined selection criteria.

Additional variables can be copied from the input plotfile to the output
without modification.


Usage
-----

.. code-block:: bash

   stream2plt stream2plt.inp


Input File
----------

The program is controlled via an input file with the following structure:

.. code-block:: none

   infile = plt000000
   outfile = plt000000_stream
   finestLevel = 2

   is_per = 0 1 0
   sym_dir = 0 0 0

   comps = 0 1 2
   no_filter = true
   nLines = 100

   Aux_Variables = Y_H2


Parameters
~~~~~~~~~~

``infile``
   Input AMReX plotfile.

``outfile``
   Output plotfile containing extracted streamlines.
   Default: ``<infile>_stream``.

``finestLevel``
   Finest AMR level used for streamline extraction.
   Default: plotfile finest level.

``is_per``
   Periodicity flags in each spatial direction
   (0: non-periodic, 1: periodic).

``sym_dir``
   Symmetry flags in each spatial direction
   (0: non-symmetric, 1: symmetric).


Streamline Variables
~~~~~~~~~~~~~~~~~~~~

``comps``
   Component indices defining the vector field used to compute
   streamlines (typically velocity components).

Alternative specification:

``sComp``
   Starting component index.

``nComp``
   Number of consecutive components forming the vector field.


Filtering Options
~~~~~~~~~~~~~~~~~

``no_filter``
   Disable all streamline filtering.
   Default: ``false``.

The following optional filters may be used when filtering is enabled:

``distComp``, ``distVal``
   Distance-based selection criterion.

``maxComps``, ``maxVals``, ``maxSgns``
   Maximum-value filtering.

``minComps``, ``minVals``, ``minSgns``
   Minimum-value filtering.

``compAt``, ``atComps``, ``valAt``, ``atVal``, ``atSgns``
   Conditional selection filters.

``RXY``, ``RXYsgn``
   Radial filtering criterion.


Additional Settings
~~~~~~~~~~~~~~~~~~~

``nLines``
   Number of streamlines extracted from the dataset.

``Aux_Variables``
   Variables copied unchanged from input to output plotfile.


Output
------

A plotfile containing sampled streamline data constructed from the
selected vector field and optional auxiliary variables.


Typical Applications
--------------------

- Flow visualization using streamlines
- Analysis of transport paths
- Flame or vortex trajectory inspection
- Reduced streamline datasets for post-processing


Notes
-----

The selected components must represent a valid vector field
(e.g. velocity). When filtering is disabled, streamlines are sampled
uniformly from available seed locations.
