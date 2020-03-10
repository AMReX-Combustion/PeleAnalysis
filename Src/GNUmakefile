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
#EBASE         = stream
EBASE         = stream2plt
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
NEEDS_f90_SRC = FALSE
#NEEDS_f90_SRC = TRUE

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

CEXE_sources += $(EBASE).cpp
ifeq ($(NEEDS_f90_SRC),TRUE)
  f90EXE_sources += $(EBASE)_nd.f90
endif
CEXE_headers += StreamData.H
CEXE_sources += StreamData.cpp

INCLUDE_LOCATIONS += .
VPATH_LOCATIONS   += .

Pdirs   := Base Boundary AmrCore Extern/amrdata
Ppack   += $(foreach dir, $(Pdirs), $(AMREX_HOME)/Src/$(dir)/Make.package)

include $(Ppack)
INCLUDE_LOCATIONS += $(Blocs)
VPATH_LOCATIONS   += $(Blocs)

CEXE_sources += AppendToPlotFile.cpp WritePlotFile.cpp
CEXE_headers += AppendToPlotFile.H   WritePlotFile.cpp
INCLUDE_LOCATIONS += $(AMREX_HOME)/Tools/C_util
VPATH_LOCATIONS += $(AMREX_HOME)/Tools/C_util

CEXE_sources += AMReX_Extrapolater.cpp
CEXE_headers += AMReX_Extrapolater.H
f90EXE_sources += AMReX_extrapolater_$(DIM)d.f90
INCLUDE_LOCATIONS += $(AMREX_HOME)/Src/Amr
VPATH_LOCATIONS += $(AMREX_HOME)/Src/Amr

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
