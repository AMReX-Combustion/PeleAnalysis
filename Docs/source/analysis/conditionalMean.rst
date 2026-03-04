
conditionalMean
===============

Description
-----------

``conditionalMean`` computes conditional averages of selected variables
from AMReX plotfiles based on a specified conditioning variable.
The tool is intended for post‑processing turbulent combustion datasets
and enables analysis such as conditional statistics with respect to
mixture fraction, temperature, or progress variables.

The executable reads AMReX plotfiles and bins data according to the
conditioning variable, computing mean values within each bin.

Usage
-----

.. code-block:: bash

   conditionalMean input.inp

Input File
----------

The program is controlled via an input file with the following structure:

.. code-block:: none

   plotfile = pltXXXXX
   condition_variable = Z
   variables = T Y(H2) Y(O2)
   nbins = 100
   output = conditionalMean.dat

Parameters
~~~~~~~~~~

``plotfile``
   Path to the AMReX plotfile.

``condition_variable``
   Variable used for conditioning (e.g. mixture fraction).

``variables``
   List of variables for which conditional means are computed.

``nbins``
   Number of bins used for the conditioning variable.

``output``
   Output filename.

Output
------

The output file contains bin centers and corresponding conditional
mean values for each requested variable.

Typical Applications
--------------------

- Conditional temperature statistics
- Mixture fraction conditioned quantities
- Turbulence–chemistry interaction analysis
- Flame structure diagnostics

Notes
-----

The conditioning variable must exist in the plotfile.
Large plotfiles may require significant memory depending
on bin resolution.

