.. highlight:: bash


rmsVel
******

The tool `rmsVel` calculates the rms velocity flcutuations over the given plotfiles as
.. math::

    u_{\mathrm{rms}} =
    \sqrt{
        \frac{
            (\langle u_x^2 \rangle - \langle u_x \rangle^2)
            +
            (\langle u_y^2 \rangle - \langle u_y \rangle^2)
            +
            (\langle u_z^2 \rangle - \langle u_z \rangle^2)
        }{3}
    }



Usage: ::

  ./rmsVel2d.gnu.MPI.ex infiles=$(ls -d plt*) [options]

Example: ::

  ./rmsVel2d.gnu.MPI.ex ./InputSamples/rmsVel.inp

Example Input File ``rmsVel.inp``::

        #------------------- IO CONTROL -----------------------------------------------------------
        infiles = plt00000 plt00001 plt00002      # pltfiles to compute fluctuations over

        #------------------- Operation control ----------------------------------------------------
        finestLevel = 4				  # DEF: -1, gets finest level by defaults, or max level user sets to consider for rms calculation
