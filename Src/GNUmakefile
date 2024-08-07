PELE_HOME ?= ../../PeleLMeX
PELE_PHYSICS_HOME ?= $(PELE_HOME)/Submodules/PelePhysics
AMREX_HOME ?= $(PELE_PHYSICS_HOME)/Submodules/amrex

# AMReX
DEBUG	      = FALSE
DIM	      = 3
PRECISION     = DOUBLE

# Compilation
COMP          = gnu
USE_MPI       = TRUE
USE_OMP       = FALSE
USE_CUDA      = FALSE

# Executable
# Uncomment the desired exec
EBASE		= template
#EBASE         = grad
#EBASE         = combinePlts
#EBASE         = isosurface
#EBASE         = flattenAMRFile
#EBASE         = filterPlt
#EBASE         = stream
#EBASE         = partStream
#EBASE         = stream2plt
#EBASE         = sampleStreamlines
#EBASE         = streamScatter
#EBASE         = streamSub
#EBASE         = streamTubeStats
#EBASE         = curvature
#EBASE         = surfMEFtoDAT
#EBASE         = surfDATtoMEF
#EBASE         = buildPMF
#EBASE          = amrToFE
#EBASE         = conditionalMean

#NEEDS FIXING:
#EBASE         = subPlt
#EBASE         = avgToPlane

NEEDS_f90_SRC = FALSE
#NEEDS_f90_SRC = TRUE

Pdirs   := Base Boundary AmrCore Extern/amrdata LinearSolvers/MLMG
Ppack   += $(foreach dir, $(Pdirs), $(AMREX_HOME)/Src/$(dir)/Make.package)

PPdirs  := Utility/PltFileManager Utility/Filter
Ppack   += $(foreach dir, $(PPdirs), $(PELE_PHYSICS_HOME)/Source/$(dir)/Make.package)

CEXE_sources += $(EBASE).cpp
ifeq ($(NEEDS_f90_SRC),TRUE)
  f90EXE_sources += $(EBASE)_nd.f90
endif

ifeq ($(EBASE),partStream)
  USE_PARTICLES = TRUE
  Pdirs += Particle
  CEXE_headers += StreamData.H    StreamPC.H
  CEXE_sources += StreamData.cpp  StreamPC.cpp
endif

ifeq ($(EBASE),isosurface)
  SDF_LOC = ../Tools/SDFGen
  include $(SDF_LOC)/Make.package
  INCLUDE_LOCATIONS += $(SDF_LOC)
  VPATH_LOCATIONS += $(SDF_LOC)
endif

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

INCLUDE_LOCATIONS += .
VPATH_LOCATIONS   += .

include $(Ppack)
INCLUDE_LOCATIONS += $(Blocs)
VPATH_LOCATIONS   += $(Blocs)

include $(AMREX_HOME)/Tools/GNUMake/Make.rules
