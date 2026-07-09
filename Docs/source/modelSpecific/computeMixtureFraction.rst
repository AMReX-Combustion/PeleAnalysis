computeMixtureFraction
======================

Description
-----------

``computeMixtureFraction`` computes the Bilger mixture fraction ``Z`` from an
AMReX plotfile.
The mixture fraction is evaluated from the elemental (C, H, O) composition of
the mixture using the Bilger formulation and the PelePhysics equation of state,
so it is consistent with the compiled chemical mechanism.

For each cell the tool computes

.. math::

   Z = \frac{\beta - \beta_\mathrm{ox}}{\beta_\mathrm{fu} - \beta_\mathrm{ox}},
   \qquad
   \beta = \sum_n \beta_n Y_n

where :math:`\beta_n` are the species Bilger coupling factors built from the
elemental composition and atomic weights, and :math:`\beta_\mathrm{fu}`,
:math:`\beta_\mathrm{ox}` are the values in the pure fuel and oxidizer
streams. The fuel stream is taken as pure fuel; the oxidizer stream defaults
to air (:math:`Y_{O_2}=0.233`, :math:`Y_{N_2}=0.767`) and can be overridden.

For a reaction progress variable, use the :doc:`../analysis/progVar` tool.


Usage
-----

.. code-block:: bash

   computeMixtureFraction infile=<s> [options]


Input File
----------

The program is controlled via command-line arguments or an input file with
the following structure:

.. code-block:: none

   infile        = plt000000
   fuelName      = H2
   finestLevel   = 2
   Aux_Variables = density temp


Parameters
~~~~~~~~~~

``infile``
   Input AMReX plotfile. It must contain the species mass fractions
   ``Y(<species>)`` for all species in the compiled mechanism.

``fuelName``
   Name of the fuel species defining the pure-fuel stream :math:`\beta_\mathrm{fu}`.
   Default: ``H2``.

``YO2ox``
   O\ :sub:`2` mass fraction of the oxidizer stream used for
   :math:`\beta_\mathrm{ox}`. Default: ``0.233`` (air).

``YN2ox``
   N\ :sub:`2` mass fraction of the oxidizer stream used for
   :math:`\beta_\mathrm{ox}`. Default: ``0.767`` (air).

``outsuffix``
   Suffix appended to the input filename to form the output filename.
   Default: ``_ZC``.

``finestLevel``
   Finest AMR level up to which the field is computed.
   Default: plotfile finest level.

``is_per``
   Periodicity flags in each spatial direction (0: non-periodic, 1: periodic).
   Default: ``1 1 1``.

``Aux_Variables``
   Names of variables copied unchanged from the input plotfile to the output
   plotfile, for carrying through fields the tool does not otherwise write.
   Default: none.

``n_files``
   Maximum number of binary files used to write the output plotfile data
   (AMReX ``VisMF::SetNOutFiles``). Lower this to reduce the number of files
   created for large parallel post-processing runs. AMReX clamps the value to
   the number of MPI ranks, so a serial run always writes a single data file.
   Default: the AMReX default.


Output
------

A new AMReX plotfile named ``<infile>outsuffix`` (by default
``<infile>_ZC``) containing:

- ``Z`` — Bilger mixture fraction
- any ``Aux_Variables`` requested, copied unchanged from the input

The output inherits the domain geometry, coordinate system, and box
structure from the input plotfile.


Typical Applications
--------------------

- Mixture-fraction conditioning of reacting-flow data
- Flamelet and manifold post-processing
- Diagnostics of partially premixed and non-premixed combustion


Notes
-----

This tool requires a PelePhysics-enabled build. The chemical mechanism is
compiled in at build time and determines the species list and the elemental
composition used to build the Bilger factors; the tool cannot be used with a
generic AMReX build.
