.. highlight:: bash


integral
********

Calculate integrals of 2D and 3D plot files in a given number of directions. Details on controling the integral directions can be found below.

Usage: ::

   ./integral3d.gnu.MPI.ex infile=<s> vars=<s, list<s>> integralDimension=<i> dir=<i> [options]

Example: ::

   ./integral3d.gnu.MPI.ex ./InputSamples/integral.inp


Tool Options
#############

::

   #------------------- IO CONTROL -----------------------------------------------------------
   infile = plt00000                          # Plot file for surface construction
   finestLevel = 0                            # DEF: finest level of plot file; Sets the finest level to read.

   #------------------- Variables ------------------------------------------------------------
   vars = I_R(H2)                             # Variables to integrate
   cVar = Y(H2)                               # Set variable to restrict integrated region.
   cMin = 1e-5                                # Cell below this value (for cVar) are skipped for integral
   cMax = 0.0111                              # Cell above this value (for cVar) are skipped for integral

`vars` takes a space-separated list of variable names. Some varnames need to be wrapped in `""`. The integration space can be restricted to a specific progress variable space based on `cVar`. Cells, where the value of `cVar` is lower than `cMin` or `cMax` are excluded from the integral.

::

   #------------------- Integral options -----------------------------------------------------
   integralDimension = 3                      # [1, 2, 3]; Integral dimension up to AMREX_SPACEDIM. For integralDimension<AMREX_SPACEDIM additional flags for the direction must be provided.
   avg = 0                                    # [0, 1], DEF: 0; If 1, divides the integral by the area.

   #---------- Only relevant for integralDimension = 1 -------------------------
   dir = 0                                    # [0, 1, 2], up to AMREX_SPACEDIM-1. Dimansion along which to integrate.

   #---------- Only relevant for integralDimension = 2 and AMREX_SPACEDIM=3 ----
   dir1 = 0                                   # [0, 1, 2]. First dimesion along which to integrate.
   dir2 = 1                                   # [0, 1, 2]. First dimesion along which to integrate.

The `integralDimension` defines in how many dimensions the integral is performed. `integralDimension = 1` calculates the integral along lines. The axis of the lines is defined via `dir`. E.g. `dir = 0` calculates line integrals along the x axis. `integralDimension = 2` calculates the integral on planes / slices. For `AMREX_SPACEDIM=2`, this is the full integral and needs no specification of the direction. For `AMREX_SPACEDIM=3`, the orientation of the planes is provided via `dir1` and `dir2`. The option `avg` results in returning the average insteag of the integral.

::

   #----------------- Additional options ----------------------------------------------------
   format = dat                               # [dat, ppm], DEF: dat; Option to create a ppm (portable pixmap) image file. Only available for AMREX_SPACEDIM=3 and integralDimension=1
   useminmax1 = -1e5 1e5                      # Minimum and maximum for normalization in ppm. Need to provide 2 values for each integrated var.
   #useminmax2 = -2e5 4e5                      # Minimum and maximum for normalization in ppm. Need to provide 2 values for each integrated var.
   goPastMax = 0                              # [0, 1], DEF: 1; For format = ppm, specify behaviour over vMax. 0: cap, 1: extra color scale over magenta to white.

The format option allows to output ASCI coded `dat` files or `ppm` (portable pixmap) image files. For the latter, `useminmax1`, `useminmax2`, and so on allow to normalize fariables 1, 2, and so on, respectively. `goPastMax` extends the colorscale to a superunitiy region.

