## Contributions Model

New contributions to *PeleAnalysis* are welcome !

To add a new feature to PeleAnalysis, the procedure is:

1. Create a branch for the new feature (locally) ::

    git checkout -b AmazingNewFeature

2. Develop the feature, merging changes often from the ``development`` branch into your ``AmazingNewFeature`` branch ::
   
    git commit -m "Developed AmazingNewFeature"
    git checkout development
    git pull                     [fix any identified conflicts between local and remote branches of "development"]
    git checkout AmazingNewFeature
    git merge development        [fix any identified conflicts between "development" and "AmazingNewFeature"]

3. Push feature branch to PeleAnalysis repository (if you have write access, otherwise fork the repo and
push the new branch to your fork)::

    git push -u origin AmazingNewFeature [Note: -u option required only for the first push of new branch]

4.  Submit a merge request through the github project page - be sure you are requesting to merge your branch to the ``development`` branch of the ``ITV-RWTH/PeleAnalysis`` and not to the main repository ``AMReX-Combustion/PeleAnalysis``.


## PeleAnalysis Coding Style Guide

Source code files can be automatically formatted to adhere to the appropriate formatting rules using ``clang-format``. To format all files, use the command:

    find Src \( -name "*.cpp" -o -name "*.H" \) -exec clang-format -i {} +

from within the PeleAnalysis base directory. You can also format files individually using ``clang-format -i /path/to/file``. Adherence to this format is checked for all PRs.

Beyond that, as much as possible, `PeleAnalysis` adheres to [AMReX Coding Style](https://github.com/AMReX-Codes/amrex/blob/development/CONTRIBUTING.md#amrex-coding-style-guide)
and we are encouraging contributors to follow those guidelines.
