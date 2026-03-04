.. highlight:: bash


avgPlotfiles
****************

The tool `avgPlotfiles` averages the given plotfiles. It is flexible with respect to the underlaying BoxArrays: the output file will be refined anywhere *any* of the input files are refined (coarse data is interpolated to the finer levels in each file as needed before averaging to obtain this result). Both tools require that all files have the same domain and base grid. The user can select a specifc list of variables, in which case that list must be present in all input files, but otherwise the input files may contain different sets of variables.

Usage: ::

  ./avgPlotfiles2d.gnu.MPI.ex infiles=$(ls -d plt*) [options]

Example: ::

  ./avgPlotfiles2d.gnu.MPI.ex ./InputSamples/avgPlotfiles.inp

Example Input File ``avgPlotfiles.inp``::

    #------------------- IO CONTROL -----------------------------------------------------------
    infiles = plt00000 plt00001 plt00002      # pltfiles to average
    outfile = plt_averaged                   # Name of output file

    #------------------- Operation control ----------------------------------------------------
    variables = temp HeatRelease             # DEF: all possible, list of variable names to average
    output_max_level = 4                     # DEF: 1000, max level to consider for average
    output_max_grid_size = 32                # DEF: 32, output max_grid_size. If all BoxArrays are the same, this is ignored
    interp_type = 1                          # DEF: 1, determines the type of interpolation when FillPatching: 0 -> piecewise constant, 1 -> cell cons linear

