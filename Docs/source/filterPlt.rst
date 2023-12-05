.. highlight:: bash


filterPlt - apply a filter to a Plot file
*****************************************

This tool utilizes the PelePhysics PltFileManager utility to read in plot files
and the PelePhysics Filter utility to apply different types of filters. To compile,
it is necessary to define the AMREX_HOME and PELE_PHYSICS_HOME
variables in the GNUmakefile. For multilevel plot files, filtering can either
be done with a constant absolute filter width, or a constant filter to grid
ratio across levels. The filter to grid ratio on the base level must always be
even. Consult `PelePhysics <https://amrex-combustion.github.io/PelePhysics/Utility.html#filter>`_
to see different filter types offered.

Usage: ::

   ./filterPlt3d.gnu.MPI.ex infile=plt00000 [options]

Help: ::

   ./filterPlt3d.gnu.MPI.ex help=true

Example:
