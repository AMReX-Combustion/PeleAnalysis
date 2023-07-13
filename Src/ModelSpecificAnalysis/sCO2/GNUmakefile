AMREX_HOME ?= ../../CCSE/amrex

DEBUG	= FALSE

DIM	= 3

COMP    = llvm

USE_CUDA     = FALSE
USE_MPI      = TRUE
USE_OMP      = FALSE
USE_EB       = TRUE
TINY_PROFILE = FALSE
BL_NO_FORT   = FALSE

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

CEXE_sources += main.cpp

include $(AMREX_HOME)/Src/Base/Make.package
include $(AMREX_HOME)/Src/Extern/amrdata/Make.package
include $(AMREX_HOME)/Src/EB/Make.package
include $(AMREX_HOME)/Src/AmrCore/Make.package
include $(AMREX_HOME)/Src/Boundary/Make.package
include $(AMREX_HOME)/Src/F_Interfaces/Base/Make.package

include $(AMREX_HOME)/Tools/GNUMake/Make.rules
