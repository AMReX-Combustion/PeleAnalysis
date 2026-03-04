.. highlight:: bash

testQPDtools
************
Print a comprehensive structural summary of the compiled chemical mechanism. Unlike other tools in this suite, ``testQPDtools`` requires no plot file input — all output is derived directly from the mechanism compiled into the executable at build time. It is primarily used to verify that the correct mechanism has been compiled, to inspect reaction connectivity and stoichiometry, and to preview the QPD edge graph that will be used by :doc:`plotQPD`.

The following information is printed to stdout in order:

.. list-table::
   :widths: 30 70
   :header-rows: 1

   * - Section
     - Description
   * - Element list
     - All elements present in the mechanism (e.g. ``C H O N``)
   * - Species list
     - All species in mechanism index order
   * - Species composition
     - Elemental composition of each species (e.g. ``CH4 = { C:1 H:4 }``)
   * - Reaction count
     - Total number of reactions (``NUM_REACTIONS``)
   * - Species in reactions
     - For each species, the indices of reactions where it appears on the left-hand side and right-hand side
   * - RMAP
     - Forward reaction index map used internally by PelePhysics
   * - RRMAP
     - Inverse reaction index map (maps from PelePhysics indices back to mechanism file indices)
   * - Stoichiometric coefficients
     - For each reaction in mapped order, the species involved and their stoichiometric coefficients
   * - Edge list
     - QPD atom-transfer edges for the chosen element (see ``trElem`` below)

Usage: ::

   ./testQPDtools.gnu.MPI.ex [options]

Example: ::

   ./testQPDtools.gnu.MPI.ex ./InputSamples/testQPDtools.inp

.. note::

   This tool requires no plot file. All output is determined entirely by the mechanism compiled into the executable. To inspect a different mechanism, recompile with the desired mechanism.

Tool Options
#############
::

   #------------------- Mechanism Inspection Options -----------------------------------------
   trElem = C                                 # DEF: C; Element for QPD edge graph construction

`trElem` is the only runtime option. It specifies which element to use when constructing the atom-transfer edge graph printed in the final section. The element must be present in the compiled mechanism; the tool will abort with an assertion error if an unrecognised element name is provided. Defaults to ``C`` (carbon), consistent with the default in :doc:`plotQPD`. Common alternatives include ``H``, ``O``, and ``N``.

All other sections (element list, species list, composition, reaction count, reaction participation, RMAP, RRMAP, and stoichiometric coefficients) are always printed regardless of any options.

Output
######
All output is written to stdout. A typical run produces output structured as follows:

.. code-block:: none

   ==> Element list
   { C H O N ... }

   ==> Species list
   { H2 H O O2 OH H2O CH4 ... }

   ==> Species composition
   H2 = { H:2 }
   O2 = { O:2 }
   CH4 = { C:1 H:4 }
   ...

   ==> NumReactions: 325

   ==> Species in reactions
   Spec, rnsL, rnsR: H2: { 0 3 7 ... } { 1 2 5 ... }
   ...

   ==> RMAP
   { 0 1 2 3 ... }

   ==> RRMAP
   { 0 1 2 3 ... }

   ==> Stoich coeffs
   rn[0] (0) = { H2:1 O:1 H:1 OH:1 }
   ...

   ==> Edges
   CH4 -- CH3 [reactions: ...]
   ...

The RMAP and RRMAP sections are particularly useful when cross-referencing reaction indices between this tool and :doc:`plotQPD`, since ``plotQPD`` accumulates fluxes using the mapped reaction indices.

Relation to Other Tools
########################
- The edge list printed by this tool is exactly the edge graph used by :doc:`plotQPD` for the same value of ``trElem`` / ``QPDatom``. Running ``mechInfo`` first is a good way to verify the edge structure before running a full QPD computation on a large plot file.
- The RMAP printed here is the same mapping applied internally by :doc:`plotQPD` when accumulating ``Qfsum`` and ``Qrsum``. Use RRMAP to trace a mapped reaction index back to its position in the original mechanism file.

Dependencies
############
This tool requires a PelePhysics-enabled build. The entire output is determined at compile time by the choice of mechanism. The tool cannot be used with a generic AMReX build.
