.. highlight:: bash


checkIso
********

The tool `checkIso` checks the correctness of an isofile, i.e., the consistent numbering of edges.

Usage: ::

  ./checkIso2d.gnu.MPI.ex isoFile=<s>

Example: ::

  ./checkIso2d.gnu.MPI.ex ./InputSamples/checkIso.inp

Example Input File ``checkIso.inp``::

        # Example input file for checkIso
        #
        # This tool checks the correctness of an isofile, i.e., the consistent numbering of edges.
        
        #------------------- IO CONTROL -----------------------------------------------------------
        isoFile = plt00000_surf                    # Isosurface file to check
