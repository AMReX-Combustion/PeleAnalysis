computeMixtureFraction
======================

Description
-----------

``computeMixtureFraction`` computes the Bilger mixture fraction ``Z`` from an
AMReX plotfile, and optionally the mass fraction of every element in the
compiled mechanism, one ``Z_<element>`` column each. The fuel stream may be a
single species or a blend of several.
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
streams. The oxidizer stream defaults to air (:math:`Y_{O_2}=0.233`,
:math:`Y_{N_2}=0.767`) and can be overridden.

The fuel stream is a single species by default, but may also be a blend given
as a list of species with their mole or mass fractions — for example an
H\ :sub:`2`/CH\ :sub:`4` mixture. Mole fractions are converted internally with

.. math::

   Y_k = \frac{X_k W_k}{\sum_j X_j W_j}

using the molecular weights of the compiled mechanism.

With ``elementMassFracs=1`` the tool additionally writes the mass fraction of
each element,

.. math::

   Z_e = \sum_n \frac{n_{e,n} A_e}{W_n} Y_n

where :math:`n_{e,n}` is the number of atoms of element :math:`e` in species
:math:`n`, :math:`A_e` the atomic weight and :math:`W_n` the molecular weight.

The element set, its size, its names and its order all come from the compiled
mechanism, so the tool writes one column per element and nothing is dropped. The
order is the mechanism's own and is *not* CHON in general — drm19, for instance,
lists its elements as ``O H C N Ar``, so the columns come out as ``Z_O``,
``Z_H``, ``Z_C``, ``Z_N``, ``Z_Ar``. Select them by name rather than by
position. Because the decomposition is complete, the columns sum to one in every
cell, to the precision of the mechanism's own atomic and molecular weights.
Unlike :math:`Z` they are unnormalised and independent of the chosen fuel and
oxidizer streams, which makes them useful as additional conditioning variables
in their own right.

.. note::

   The Bilger coupling function is built on C, H and O alone, so :math:`Z`
   itself accounts for no other element. That is deliberate for N — the
   formulation excludes it even where the mechanism carries NOx chemistry — and
   harmless for diluents such as Ar and He. Where a mechanism bonds some other
   element to C, H or O, that element is part of the combustion chemistry and
   :math:`Z` cannot see it; the tool prints a warning naming the element and an
   example species. The ``alzeta`` mechanism, which contains fluorine, triggers
   it. The elemental mass fractions are unaffected either way.

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

   infile           = plt000000
   fuelNames        = H2 CH4
   fuelMoleFracs    = 0.95 0.05
   elementMassFracs = 1
   finestLevel      = 2
   Aux_Variables    = density temp


Parameters
~~~~~~~~~~

``infile``
   Input AMReX plotfile. It must contain the species mass fractions
   ``Y(<species>)`` for all species in the compiled mechanism.

``fuelName``
   Name of the fuel species defining the pure-fuel stream
   :math:`\beta_\mathrm{fu}`, when the fuel is a single species.
   Default: ``H2``.

``fuelNames``
   Species making up a blended fuel stream, e.g. ``H2 CH4``. Supply the
   composition with either ``fuelMoleFracs`` or ``fuelMassFracs``. Overrides
   ``fuelName`` when given.

``fuelMoleFracs``
   Mole fractions of the ``fuelNames`` species. Normalised internally, so
   they need not sum exactly to one.

``fuelMassFracs``
   Mass fractions of the ``fuelNames`` species, as an alternative to
   ``fuelMoleFracs``. Giving both is an error.

``elementMassFracs``
   Set to ``1`` to also write one ``Z_<element>`` column for every element in
   the compiled mechanism. The number of columns and their names therefore
   depend on the mechanism. Default: ``0``.

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
- ``Z_<element>``, one per element in the compiled mechanism and in the
  mechanism's own order — elemental mass fractions, only when
  ``elementMassFracs=1``
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
