DEBUG	      = FALSE
DIM	      = 3
COMP          = intel
PRECISION     = DOUBLE
USE_ASSERTION = TRUE
USE_MPI       = TRUE
USE_OMP       = FALSE
USE_EB        = TRUE
EBASE         = main

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

CEXE_sources += $(EBASE).cpp

CEXE_sources += AMReX_EB2_IF_Triangulated.cpp
CEXE_headers += AMReX_EB2_IF_Triangulated.H

INCLUDE_LOCATIONS += .
VPATH_LOCATIONS   += .

Pdirs   := Base Boundary EB AmrCore
Ppack   += $(foreach dir, $(Pdirs), $(AMREX_HOME)/Src/$(dir)/Make.package)

include $(Ppack)
INCLUDE_LOCATIONS += $(Blocs)
VPATH_LOCATIONS   += $(Blocs)

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
