.. highlight:: bash


dumpFABslice
************

This tool reads a plotfile and outputs a FAB slice with a specified normal and index position.

Usage: ::

   ./dumpFABslice3d.gnu.MPI.ex infile=plt00000 [options]

Help: ::

   ./dumpFABslice3d.gnu.MPI.ex help=true

Example: ::

   ./dumpFABslice3d.gnu.MPI.ex ./InputSamples/dumpFABslice.inp

Example Input File ``dumpFABslice.inp``::

        #------------------- IO CONTROL -----------------------------------------------------------
        infile = plt00000                  # AMReX plotfile to slice

        #------------------- Slice control --------------------------------------------------------
        dir = 2                            # Direction to slice (0=x, 1=y, 2=z), (3D only, REQ)
        num = 32                           # Plane index along dir to extract,    (3D only, REQ)
        
        #------------------- Variable control -----------------------------------------------------
        varNames = temp x_velocity         # DEF: all variables, list of variables to slice
        finestLevel = 2                    # DEF: finest available, finest AMR level to use
