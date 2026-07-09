.. highlight:: bash


combinePlts
***********

Create a new plotfile that is composed of a set of components taken from each of the infiles, which are existing plotfiles. The input plotfiles *must* have the same AMR hierarchy,
up to the finest level requested for the output.

Usage: ::

  ./combinePlts2d.gnu.MPI.ex infiles=<s1 s2 s3> outfile=<s> vars=<s1 s2 s3> [options] 

Example: ::

  ./combinePlts2d.gnu.MPI.ex ./InputSamples/combinePlts.inp

Example Input File ``combinePlts.inp``::

        # Example input file for combinePlts
        
        #------------------- IO CONTROL -----------------------------------------------------------
        infiles = plt00000 plt00001 plt00002      # pltfiles to combine
        outfile = plt_combined			  # Name of output file
        n_files = 64                              # DEF: AMReX default, cap on the number of plotfile data files (VisMF), clamped to the MPI rank count

        #------------------- Operation control ----------------------------------------------------
        vars = temp HeatRelease       	          # DEF: all possible, list of variables to average COMP or NAME
        finestLevel = 4				  # DEF: -1, gets finest level by defaults, or max level user sets to consider for combine 
        is_per = 1 1 0                            # DEF: 0 0 0, sets case periodicity 

