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
DIM	      = 2
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
#EBASE   = avgPlotfiles
#EBASE   = checkIso
#EBASE   = combinePlts
#EBASE   = conditionalMean
#EBASE   = curvature
#EBASE   = filterPlt
#EBASE   = flattenAMRFile
#EBASE   = grad
#EBASE   = integral
#EBASE   = isoMEF
#EBASE   = isosurface
#EBASE   = jpdf
#EBASE   = makePlotfile
#EBASE   = partStream
#EBASE   = qCriterion
#EBASE   = rmsVel
#EBASE   = scaleMEF
#EBASE   = stream2plt
#EBASE   = streamBinTubeStats
#EBASE   = subPlt
#EBASE   = surfMEFtoDAT
#EBASE   = surfDATtoMEF
#EBASE   = trimMEFgen

# tools that need review/combining/fixing/work.  Comment with double ## t escape from CI
##EBASE   = avgToPlane
##EBASE   = convert2hdf5
##EBASE   = combineMEF
##EBASE   = decimateMEF
##EBASE   = dumpFABslice
##EBASE   = interp
##EBASE   = mergeMEF
##EBASE   = multMEF
##EBASE   = regridPlt
##EBASE   = sliceMEF
##EBASE   = smoothMEF
##EBASE   = turbfile



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
