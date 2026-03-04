plotYtoX
========

Description
-----------

``plotYtoX`` converts species mass fractions stored in an AMReX plotfile
to mole fractions using the PelePhysics equation of state. The tool reads
species mass fractions ``Y(species)`` together with temperature and writes
a new plotfile containing the corresponding mole fractions ``X(species)``.

The spatial structure and AMR hierarchy of the original dataset are
preserved.


Usage
-----

.. code-block:: bash

   plotYtoX plotYtoX.inp


Input File
----------

The program is controlled via an input file with the following structure:

.. code-block:: none

   infile = plt000000
   finestLevel = 2


Parameters
~~~~~~~~~~

``infile``
   Input AMReX plotfile containing species mass fractions
   ``Y(species)`` and temperature ``Temp``.

``finestLevel``
   Finest AMR level processed by the conversion.
   Default: finest level available in the plotfile.


Output
------

A new AMReX plotfile containing:

- mole fractions ``X(species)``
- temperature ``Temp``

The output plotfile name is automatically generated as

.. code-block:: none

   <infile>_X


Typical Applications
--------------------

- Conversion of PeleLMeX output from mass to mole fractions
- Thermochemical post-processing
- Preparation of datasets requiring mole-based quantities


Notes
-----

The input plotfile must contain species mass fractions named
``Y(species)`` as well as the temperature field ``Temp``.
The output filename cannot be specified manually and is
generated automatically from the input filename.
