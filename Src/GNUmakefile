AMREX_HOME ?= ../../..

DEBUG	      = FALSE
DIM	      = 2
COMP          = gcc
PRECISION     = DOUBLE
USE_MPI       = TRUE
USE_OMP       = FALSE
#EBASE         = grad
#EBASE         = combinePlts
#EBASE         = isosurface
EBASE         = stream
#EBASE         = sampleStreamlines
#EBASE         = streamScatter
#EBASE         = streamSub
#EBASE         = streamTubeStats
#EBASE         = curvature
#EBASE         = surfMEFtoDAT
#EBASE         = surfDATtoMEF
#EBASE         = buildPMF
#EBASE         = subPlt
#EBASE         = conditionalMean
#EBASE         = jpdf
#NEEDS_f90_SRC = FALSE
NEEDS_f90_SRC = TRUE

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

CEXE_sources += $(EBASE).cpp
ifeq ($(NEEDS_f90_SRC),TRUE)
  f90EXE_sources += $(EBASE)_nd.f90
endif

INCLUDE_LOCATIONS += .
VPATH_LOCATIONS   += .

Pdirs   := Base Boundary Amr AmrCore Extern/amrdata
Ppack   += $(foreach dir, $(Pdirs), $(AMREX_HOME)/Src/$(dir)/Make.package)

include $(Ppack)
INCLUDE_LOCATIONS += $(Blocs)
VPATH_LOCATIONS   += $(Blocs)

CEXE_sources += AppendToPlotFile.cpp WritePlotFile.cpp
CEXE_headers += AppendToPlotFile.H   WritePlotFile.cpp
INCLUDE_LOCATIONS += $(AMREX_HOME)/Tools/C_util
VPATH_LOCATIONS += $(AMREX_HOME)/Tools/C_util

SDF_LOC = ../Tools/SDFGen
include $(SDF_LOC)/Make.package
INCLUDE_LOCATIONS += $(SDF_LOC)
VPATH_LOCATIONS += $(SDF_LOC)

vpath %.c   : $(VPATH_LOCATIONS)
vpath %.h   : $(VPATH_LOCATIONS)
vpath %.cpp : $(VPATH_LOCATIONS)
vpath %.H   : $(VPATH_LOCATIONS)
vpath %.F   : $(VPATH_LOCATIONS)
vpath %.f   : $(VPATH_LOCATIONS)
vpath %.f90 : $(VPATH_LOCATIONS)

include $(AMREX_HOME)/Tools/GNUMake/Make.rules
