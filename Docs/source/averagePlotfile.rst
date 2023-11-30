.. highlight:: bash


averagePlotfile - Create a pltfile by averaging existing pltfiles
*****************************************************************

Two tools exist for the purpose of averaging plotfiles. `averagePlotfile`
assumes all plotfiles to be averaged have the same underlying BoxArrays, or
the output is only given at the base level. `averagePlotfileFlexible` relaxes
this assumption: the output file will be refined anywhere *any* of the input
files are refined (coarse data is interpolated to the finer levels in each file
as needed before averaging to obtain this result). Both tools require that all
files have the same domain and base grid, and by default the same variables.
`averagePlotfileFlexible` adds the capability to optionally select a specifc
list of variables, in which case that list must be present in all input files
but otherwise the input files may contain different sets of variables.

Usage: ::

  ./avgPlotfilesFlexible2d.gnu.MPI.ex infiles=$(ls -d plt*) [options]

Example:
