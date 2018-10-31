
PeleAnalysis
============

this repository contains a collection of standalone routines for processing plotfiles created with the AMReX software framework for block-structure adaptive mesh refinement simulations.  AMReX is available at ::

    https://github.com/AMReX-Codes/amrex

In order to build these processing tools, you should clone or fork the amrex repository, and set the environment variable `AMREX_HOME` to point to the local folder where that is placed.  Then clone this repository, `cd Src` and edit the `GNUmakefile` to select which tool to build.  If AMReX is configured properly, a stand-alone executable will be built locally, based on the selected options, including spatial dimension (2 or 3), compiler choices, whether to build with MPI and/or OpenMP enabled, and whether to build a debugging or optimized version.  Note that some of the tools require building a companion `f90` source file - you must manually set the flag in the `GNUmakefile` accordingly.  More extensive documentation is available (see building instructions below).

Contributions
-------------

To add a new feature to PeleAnalysis, the procedure is:

1. Create a branch for the new feature (locally) ::

    git checkout -b AmazingNewFeature

2. Develop the feature, merging changes often from the development branch into your AmazingNewFeature branch ::
   
    git commit -m "Developed AmazingNewFeature"
    git checkout development
    git pull                     [fix any identified conflicts between local and remote branches of "development"]
    git checkout AmazingNewFeature
    git merge development        [fix any identified conflicts between "development" and "AmazingNewFeature"]

3. Push feature branch to PeleAnalysis repository (if you have write access, otherwise fork the repo and
push the new branch to your fork)::

    git push -u origin AmazingNewFeature [Note: -u option required only for the first push of new branch]

4.  Submit a merge request through the github project page - be sure you are requesting to merge your branch to the development branch.




Documentation
-------------
Documentation for the analysis routines exists in the Docs directory. To build the documentation::

    cd Docs
    make html


Acknowledgment
--------------
This research was supported by the Exascale Computing Project (ECP), Project
Number: 17-SC-20-SC, a collaborative effort of two DOE organizations -- the
Office of Science and the National Nuclear Security Administration --
responsible for the planning and preparation of a capable exascale ecosystem --
including software, applications, hardware, advanced system engineering, and
early testbed platforms -- to support the nation's exascale computing
imperative.
