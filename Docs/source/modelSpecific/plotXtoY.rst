plotXtoY
========

Description
-----------

``plotXtoY`` converts species mole fractions stored in an AMReX plotfile
to mass fractions using the PelePhysics equation of state. The tool reads
species mole fractions ``X(species)`` together with temperature and writes
a new plotfile containing the corresponding mass fractions ``Y(species)``.

The spatial structure and AMR hierarchy of the original dataset are
preserved.


Usage
-----

.. code-block:: bash

   plotXtoY plotXtoY.inp


Input File
----------

The program is controlled via an input file with the following structure:

.. code-block:: none

   infile = plt000000
   finestLevel = 2


Parameters
~~~~~~~~~~

``infile``
   Input AMReX plotfile containing species mole fractions
   ``X(species)`` and temperature ``Temp``.

``finestLevel``
   Finest AMR level processed by the conversion.
   Default: finest level available in the plotfile.


Output
------

A new AMReX plotfile containing:

- mass fractions ``Y(species)``
- temperature ``Temp``

The output plotfile name is automatically generated as

.. code-block:: none

   <infile>_plotXtoY


Typical Applications
--------------------

- Conversion of PeleLMeX output from mole to mass fractions
- Preparation of datasets for combustion diagnostics
- Post-processing requiring conserved species quantities


Notes
-----

The input plotfile must contain species mole fractions named
``X(species)`` as well as the temperature field ``Temp``.
The output filename cannot be specified manually and is
generated automatically from the input filename.
