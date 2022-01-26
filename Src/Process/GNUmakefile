PRECISION  = DOUBLE
PROFILE    = FALSE
DEBUG      = FALSE

#DIM        = 1
#DIM        = 2
DIM        = 3

COMP	   = gcc
USE_MPI    = TRUE
USE_OMP    = FALSE

EBASE = process

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

CEXE_sources += $(EBASE).cpp

Pdirs 	:= Base Boundary AmrCore Extern/amrdata

Bpack	+= $(foreach dir, $(Pdirs), $(AMREX_HOME)/Src/$(dir)/Make.package)
Blocs	+= $(foreach dir, $(Pdirs), $(AMREX_HOME)/Src/$(dir))

include $(Bpack)

INCLUDE_LOCATIONS += $(Blocs)
VPATH_LOCATIONS   += $(Blocs)

AMREX_UTILS_HOME = $(AMREX_HOME)/Tools/C_util
include $(AMREX_UTILS_HOME)/Make.package
INCLUDE_LOCATIONS += $(AMREX_UTILS_HOME)
VPATH_LOCATIONS   += $(AMREX_UTILS_HOME)

include $(AMREX_HOME)/Tools/GNUMake/Make.rules

print-%: ; @echo $* is $($*) ; $(SILENT) $(RM) buildInfo.cpp
