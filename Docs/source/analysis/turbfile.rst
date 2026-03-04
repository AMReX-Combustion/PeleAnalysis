.. highlight:: bash


turbfile
*********

This tool reads a plotfile and outputs a directory ofile including /HDR and /DAT containing slices with the normal direction dir. 

Usage: ::

   ./turbfile3d.gnu.MPI.ex infile=plt00000 [options]

Help: ::

   ./turbfile3d.gnu.MPI.ex help=true

Example: ::

   ./turbfile3d.gnu.MPI.ex ./InputSamples/turbfile.inp

Example Input File ``turbfile.inp``::


        #------------------- IO CONTROL -----------------------------------------------------------
        ifile = plt00000		        # pltfile to read for creating the inflow file
        ofile = turbInflow			# Name of output directory including /HDR and /DAT
        
        #------------------- Operation control ----------------------------------------------------
        dir = 2			                # DEF: AMREX_SPACEDIM -1, sets normal direction of inflow plane
