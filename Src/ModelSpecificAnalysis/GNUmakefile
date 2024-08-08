PELE_ANALYSIS_HOME = ../..

# AMReX
DIM        = 3
PRECISION  = DOUBLE
PROFILE    = FALSE
VERBOSE    = TRUE
DEBUG      = FALSE

# Compiler
COMP	   = gnu
USE_MPI    = TRUE
USE_OMP    = FALSE
USE_CUDA   = FALSE
USE_HIP    = FALSE
TINY_PROFILE = FALSE

EBASE=bruno
EBASE=plotXtoY
EBASE=plotYtoX
EBASE=testQPDtools
EBASE=plotTransportCoeff

CEXE_sources += $(EBASE).cpp

Eos_Model       = Fuego
Chemistry_Model = drm19
Transport_Model = Simple

Blocs   := .

include $(PELE_ANALYSIS_HOME)/Tools/GNUmake/Make.ModelSpecific
