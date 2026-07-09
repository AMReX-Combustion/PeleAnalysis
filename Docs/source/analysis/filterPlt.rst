.. highlight:: bash


filterPlt
*********

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

Example: ::

   ./filterPlt3d.gnu.MPI.ex ./InputSamples/filterPlt.inp

Example Input File ``filterPlt.inp``::

        #------------------- IO CONTROL -----------------------------------------------------------
        infiles = plt00000 plt00001 plt00002      # pltfiles to average
        outfile = plt_averaged			  # DEF: plt_averaged, Name of output file
        n_files = 64                              # DEF: AMReX default, cap on the number of plotfile data files (VisMF), clamped to the MPI rank count

        #------------------- Operation control ----------------------------------------------------
        variables = temp HeatRelease              # DEF: all possible, list of variable names to average
        max_filter_level = 4			  # DEF: 1000, max level to consider for filtering
        filter_type = 1				  # DEF: 1, filter type as defined in PeleC (1->box, 2->Gaussian, etc)
        base_fgr = 2				  # DEF: 2,is the desired filter to grid ratio on the base level, must be even
        same_fgr_all_levels = false		  # DEF: false, if true the same filter to grid ratio is kept on all levels (rather than absolute filter width)

        max_grid_size = 32			  # DEF: 32, output max_grid_size. 
        interp_type = 1				  # DEF: 1, determines the type of interpolation when FillPatching: 0 -> piecewise constant, 1 -> cell cons linear

