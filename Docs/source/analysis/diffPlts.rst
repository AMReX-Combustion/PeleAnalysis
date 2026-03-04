.. highlight:: bash


diffPlts
********

Create a new plotfile that includes the difference of the reference and the change plotfile for the specified variables. The input plotfiles *must* have the same AMR hierarchy,
up to the finest level requested for the output.

Usage: ::

  ./diffPlts2d.gnu.MPI.ex infiles=<s1 s2 s3> outfile=<s> vars=<s1 s2 s3> [options] 

Example: ::

  ./diffPlts2d.gnu.MPI.ex ./InputSamples/diffPlts.inp

Example Input File ``diffPlts.inp``::

        # Tool to calculate the difference of the same field from two plotfiles
        
        #------------------- IO CONTROL -----------------------------------------------------------
        infiles = plt00000 plt00000_new                       # Plt files to calculate diff. First file is the reference file (for relative difference)
        outfile = plt00000_diff                               # DEF: infiles[0]+"_diff"; Name of output plt file.
        finestLevel = 1                                       # DEF: joint finest level of both plt files.
        is_per = 0 0 0                                        # DEF: 0 0 0, Assumed periodicity
        
        #------------------- Operation control ----------------------------------------------------
        vars = Y(H2) HeatRelease                              # List of variables for diff
        diff_type = relative                                  # [absolute, relative], DEF absolute; Calculate absolute or relative difference.


