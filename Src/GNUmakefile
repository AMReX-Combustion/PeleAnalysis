PELE_ANALYSIS_HOME = ..

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

#BROKEN
#EBASE         = avgToPlane

#STREAM TOOLS - need combining/fixing/work
#EBASE         = stream
#EBASE         = partStream
#EBASE         = stream2plt
#EBASE         = sampleStreamlines
#EBASE         = streamScatter
#EBASE         = streamSub
#EBASE         = streamTubeStats



CEXE_sources += $(EBASE).cpp

ifeq ($(EBASE),partStream)
  USE_PARTICLES = TRUE
  CEXE_headers += StreamData.H    StreamPC.H
  CEXE_sources += StreamData.cpp  StreamPC.cpp
endif

ifeq ($(EBASE),$(filter $(EBASE), isosurface buildDistance))
   NEED_DISTANCE = TRUE
endif

ifeq ($(EBASE),decimateMEF) #make qslim first!
   NEED_QSLIM = TRUE
endif

Blocs := .

include $(PELE_ANALYSIS_HOME)/Tools/GNUmake/Make.Analysis
