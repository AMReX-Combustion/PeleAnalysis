.. highlight:: bash


grad
****

Given a plotfile that contains a scalar quantity, compute the components of the gradient, and its
magnitude, of that scalar at all cells in the solution. It creates a new plotfile
containing only these computed quantities, plus the specified Aux_Variables.

Usage: ::

  ./grad2d.gnu.MPI.ex infile=<s> gradVar=<s> [options]

Example: ::

  ./grad2d.gnu.MPI.ex ./InputSamples/grad.inp

Example Input File ``grad.inp``::

        #------------------- IO CONTROL -----------------------------------------------------------
        infile = plt00000                         # input pltfile to calcuate grad on
        outfile = plt_grad			  # DEF: <infile>_gt, Name of output file
        n_files = 64                              # DEF: AMReX default, cap on the number of plotfile data files (VisMF), clamped to the MPI rank count

        #------------------- Operation control ----------------------------------------------------
        gradVar = temp		      	          # variable to calculate gradient on
        finestLevel = 4				  # DEF: 1000, max level to consider for gradient 
        Aux_Variables = density x_velocity	  # DEF: None, variables that are just carried through the script and written to the output file   
        sym_dir = 0 0 0				  # DEF: 0 0 0, sets case symmetry 
        is_per = 1 1 0                            # DEF: 0 0 0, sets case periodicity
