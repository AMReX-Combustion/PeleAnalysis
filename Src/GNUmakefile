# PELE_ANALYSIS_HOME should be pointing to the home directory of your PeleAnalysis clone. 
# If the relative path does not compile, try using the absolute path.
PELE_ANALYSIS_HOME = ..

# If you recursively cloned the repository, you should not need to set the paths for 
# PELE_PHYSICS_HOME and AMREX_HOME. If you chose to not recursively clone, or you want 
# to link PeleAnalysis to your versions of PELE_PHYSICS_HOME and AMREX_HOME used by 
# PeleLMeX or PeleC, make the following lines to point to the respective directories.
#
# PELE_PHYSICS_HOME = ../../PeleLMeX/Submodules/PelePhysics
# AMREX_HOME = ../../PeleLMeX/Submodules/PelePhysics/Submodules/amrex

# AMReX
DEBUG	      = FALSE
DIM	      = 3
PRECISION     = DOUBLE

# Compilation
COMP          = gnu
USE_MPI       = FALSE
USE_OMP       = FALSE
USE_CUDA      = FALSE

# Executable
# Uncomment the desired exec
EBASE		= template
#EBASE		= convert2hdf5
#EBASE         = grad
#EBASE		= jpdf
#EBASE         = combinePlts
#EBASE         = isosurface
#EBASE         = flattenAMRFile
#EBASE         = filterPlt
#EBASE         = curvature
#EBASE         = surfMEFtoDAT
#EBASE         = surfDATtoMEF
#EBASE          = amrToFE
#EBASE         = conditionalMean
#EBASE		= isosurface
#EBASE		= decimateMEF
#EBASE         = subPlt
#EBASE         = partStream
#EBASE		= stream2plt
#EBASE		= makePlotfile
#EBASE		= dumpFABslice
#EBASE		= turbfile

# tools that need review/combining/fixing/work
#EBASE         = avgToPlane


CEXE_sources += $(EBASE).cpp

ifeq ($(EBASE),partStream)
  USE_PARTICLES = TRUE
  CEXE_headers += StreamPC.H
  CEXE_sources += StreamPC.cpp
endif

ifeq ($(EBASE),stream2plt)
  CEXE_headers += StreamData.H
  CEXE_sources += StreamData.cpp
endif

ifeq ($(EBASE),$(filter $(EBASE), isosurface buildDistance))
   NEED_DISTANCE = TRUE
endif

ifeq ($(EBASE),convert2hdf5)
   USE_HDF5 = TRUE
   USE_HDF5_ZFP = TRUE
   HDF5_HOME = /path/to/hdf5
   ZFP_HOME = /path/to/zfp
   H5Z_HOME = /path/to/h5z-zfp
endif

ifeq ($(EBASE),decimateMEF) #make qslim first!
   NEED_QSLIM = TRUE
endif

Blocs := .

include $(PELE_ANALYSIS_HOME)/Tools/GNUmake/Make.Analysis
