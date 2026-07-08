# PELE_ANALYSIS_HOME should be pointing to the home directory of your PeleAnalysis clone. 
# If the relative path does not compile, try using the absolute path.
PELE_ANALYSIS_HOME = ..

# If you recursively cloned the repository, you should not need to set the paths for 
# PELE_PHYSICS_HOME and AMREX_HOME. If you chose to not recursively clone, or you want 
# to link PeleAnalysis to your versions of PELE_PHYSICS_HOME and AMREX_HOME used by 
# PeleLMeX or PeleC, make the following lines to point to the respective directories.
#
#PELE_PHYSICS_HOME = ../../PelePhysics
#AMREX_HOME = ../../PelePhysics/Submodules/amrex


# AMReX
DEBUG	      = FALSE
DIM	      	  = 2
PRECISION     = DOUBLE

# Compilation
COMP          = gnu
USE_MPI       = FALSE
USE_OMP       = FALSE
USE_CUDA      = FALSE

# Executable
# Uncomment the desired exec
#EBASE   = template
#EBASE   = amrToFE
#EBASE   = arithmetics
#EBASE   = avgPlotfiles
#EBASE   = avgToPlane
#EBASE   = checkIso
#EBASE   = ciao2plt
#EBASE   = combinePlts
#EBASE   = conditionalMean
#EBASE   = convert2hdf5
#EBASE   = curvature
#EBASE   = diffPlts
#EBASE   = dumpFABslice
#EBASE   = favreAvgPlotfiles
#EBASE   = filterPlt
#EBASE   = flattenAMRFile
#EBASE   = generateTestPlt
#EBASE   = grad
#EBASE   = integral
#EBASE   = isoMEF
#EBASE   = isosurface
#EBASE   = jpdf
#EBASE   = makePlotfile
#EBASE   = partStream
#EBASE   = progVar
#EBASE   = qCriterion
#EBASE   = regridPlt
#EBASE   = rmsVel
#EBASE   = scaleMEF
#EBASE   = stream2plt
#EBASE   = streamBinTubeStats
#EBASE   = subPlt
#EBASE   = surfMEFtoDAT
#EBASE   = surfDATtoMEF
#EBASE   = trimMEFgen
#EBASE   = turbfile

# tools that need review/combining/fixing/work.  Comment with double ## t escape from CI
##EBASE   = combineMEF
##EBASE   = decimateMEF
##EBASE   = mergeMEF
##EBASE   = multMEF
##EBASE   = sliceMEF
##EBASE   = smoothMEF



CEXE_sources += $(EBASE).cpp

ifeq ($(EBASE),$(filter $(EBASE), template favreAvgPlotfiles)) # add targets that need utils as space-separated list
   USE_UTILS = TRUE
endif

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
   USE_HDF5_ZFP = FALSE
   HDF5_HOME = $(HDF5_DIR) # /path/to/hdf5
   ZFP_HOME = $(ZFP_DIR)  # /path/to/zfp
   H5Z_HOME = $(H5Z_DIR) # /path/to/h5z-zfp
endif

ifeq ($(EBASE),ciao2plt)
  # Serial build: use HDF5_DIR.
  # MPI build: use HDF5_PAR_DIR (parallel HDF5, compiled with MPI support).
  # If HDF5_PAR_DIR is not set, HDF5_DIR is used for both modes.
  # Library names are the same (-lhdf5 -lhdf5_cpp); only the installation path differs.
  # On some systems the parallel build uses distinct names (e.g. -lhdf5_par); override
  # LIBRARIES manually in that case.
  ifeq ($(USE_MPI),TRUE)
    HDF5_CUR_DIR ?= $(HDF5_PAR_DIR)
    HDF5_CUR_DIR ?= $(HDF5_DIR)
  else
    HDF5_CUR_DIR = $(HDF5_DIR)
  endif
  LIBRARIES += -lhdf5 -lhdf5_hl -lz -lhdf5_cpp
  INCLUDE_LOCATIONS += $(HDF5_CUR_DIR)/include
  LIBRARY_LOCATIONS += $(HDF5_CUR_DIR)/lib
endif

ifeq ($(EBASE),decimateMEF) #make qslim first!
   NEED_QSLIM = TRUE
endif

Blocs := .

include $(PELE_ANALYSIS_HOME)/Tools/GNUmake/Make.Analysis
