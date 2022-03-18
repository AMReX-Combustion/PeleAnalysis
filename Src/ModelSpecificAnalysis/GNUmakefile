PELE_ANALYSIS_HOME = ../../

# AMReX
DIM        = 3
PRECISION  = DOUBLE
PROFILE    = FALSE
VERBOSE    = TRUE
DEBUG      = FALSE

# Compiler
COMP	   = gnu
USE_MPI    = FALSE
USE_OMP    = FALSE
USE_CUDA   = FALSE
USE_HIP    = FALSE
TINY_PROFILE = FALSE

EBASE=bruno
EBASE=plotXtoY
EBASE=plotYtoX
EBASE=plotTransportCoeff
EBASE=plotQPD
EBASE=testQPDtools
NEEDS_f90_SRC = FALSE

CEXE_sources += $(EBASE).cpp
ifeq ($(NEEDS_f90_SRC),TRUE)
  f90EXE_sources += $(EBASE)_nd.f90
endif

Eos_Model       = Fuego
Chemistry_Model = drm19
Transport_Model = Simple

Blocs   := .

include $(PELE_ANALYSIS_HOME)/Tools/GNUmake/Make.Analysis
