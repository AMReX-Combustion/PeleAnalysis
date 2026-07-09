progVar
=======

Description
-----------

``progVar`` computes a reaction progress variable and its chemical source
term from an AMReX plotfile. Given a user-selected set of species, the
progress variable is defined as the normalized sum of their mass fractions

.. math::

   C = \frac{\sum_i Y_i - C_\mathrm{unburnt}}{C_\mathrm{burnt} - C_\mathrm{unburnt}}

where :math:`C_\mathrm{unburnt}` and :math:`C_\mathrm{burnt}` are the sum of
the selected species mass fractions in the unburnt and fully burnt states,
respectively. The corresponding source term is obtained by summing the
production rates ``I_R(<species>)`` of the selected species and normalizing
by the same denominator:

.. math::

   \dot{\omega}_C = \frac{\sum_i I_R(Y_i)}{C_\mathrm{burnt} - C_\mathrm{unburnt}}

The tool reads a plotfile, evaluates these quantities at every cell on every
AMR level, and writes a new plotfile containing the raw species sum, the
progress variable, and its source term. The spatial structure and AMR
hierarchy of the original dataset are preserved.


Usage
-----

.. code-block:: bash

   progVar infile=<s> speciesNames=<s> unburntVal=<f> burntVal=<f> [options]


Input File
----------

The program is controlled via command-line arguments or an input file with
the following structure:

.. code-block:: none

   infile        = plt00000
   speciesNames  = Y(H2) Y(H2O)
   unburntVal    = 0.0
   burntVal      = 1.0
   outsuffix     = _prog
   outname       = progVar
   Aux_Variables = density temp


Parameters
~~~~~~~~~~

``infile``
   Path to the input AMReX plotfile. It must contain the mass fractions of
   all selected species, and their production rates ``I_R(<species>)`` when
   the source term is written (``printSource=1``, the default).

``speciesNames``
   List of species mass-fraction fields to include in the progress
   variable, e.g. ``Y(H2) Y(H2O)``. The matching production-rate fields
   ``I_R(<species>)`` are located automatically (a leading ``Y(`` and
   trailing ``)`` are stripped when forming the ``I_R(...)`` name).

``unburntVal``
   Sum of the selected species mass fractions in the unburnt mixture.
   Default: ``0``.

``burntVal``
   Sum of the selected species mass fractions in the fully burnt mixture.
   Default: ``1``.

``outsuffix``
   Suffix appended to the input filename to form the output filename.
   Default: ``_prog``.

``outname``
   Name of the progress-variable field in the output plotfile. The source
   term is written as ``I_R(<outname>)``.
   Default: ``progVar``.

``printSource``
   Whether to compute and write the progress-variable source term
   ``I_R(<outname>)``. Set to ``0`` to output only ``specSum`` and the
   progress variable, in which case the input plotfile does not need the
   ``I_R(<species>)`` fields. Default: ``1`` (source term written).

``Aux_Variables``
   Names of variables copied unchanged from the input plotfile to the output
   plotfile. Useful for carrying through fields (e.g. ``density``, ``temp``)
   that ``progVar`` does not otherwise write. Default: none.

``finestLevel``
   Finest AMR level up to which the progress variable is computed.
   Default: plotfile finest level.

``is_per``
   Periodicity flags in each spatial direction (0: non-periodic, 1: periodic).
   Default: ``1 1 1``.


Output
------

A new AMReX plotfile named ``<infile>outsuffix`` containing:

- ``specSum`` — the sum of the selected species mass fractions
- ``<outname>`` — the normalized progress variable :math:`C`
- ``I_R(<outname>)`` — the progress-variable source term (only when
  ``printSource=1``, the default)
- any ``Aux_Variables`` requested, copied unchanged from the input

The output inherits the domain geometry, coordinate system, and box
structure from the input plotfile.


Typical Applications
--------------------

- Definition of a progress variable for flamelet/manifold analysis
- Conditional averaging in progress-variable space
- Post-processing of reacting-flow datasets from PeleLMeX or PeleC


Notes
-----

The input plotfile must contain the selected species mass fractions, and —
when ``printSource=1`` (the default) — their production rates
``I_R(<species>)``; the tool aborts if any requested field is missing. The
progress variable and its source term are normalized by
``burntVal - unburntVal``, so these two values must not be equal (the tool
aborts if they are).
