PELE_ANALYSIS_HOME ?= $(HOME)/src/CCSE/Pele/PeleAnalysis
PELE_PHYSICS_HOME ?= $(HOME)/src/CCSE/Pele/PelePhysics

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
FUEGO_GAS  = TRUE

EBASE=bruno
EBASE=plotXtoY
EBASE=plotYtoX
EBASE=testQPDtools
NEEDS_f90_SRC = FALSE

CEXE_sources += $(EBASE).cpp
ifeq ($(NEEDS_f90_SRC),TRUE)
  f90EXE_sources += $(EBASE)_nd.f90
endif

#######################
# ODE solver OPTIONS: DVODE (default) / SUNDIALS / RK explicit
#######################
# Activates use of SUNDIALS: CVODE (default) / ARKODE
USE_SUNDIALS_PP = TRUE
USE_RK64_PP = FALSE
ifeq ($(USE_SUNDIALS_PP), TRUE)
  USE_ARKODE_PP = FALSE
  ifeq ($(USE_CUDA), TRUE)
      USE_KLU_PP = FALSE
  else
      ifeq ($(USE_HIP), TRUE)
          USE_KLU_PP = FALSE
      else
          USE_KLU_PP = FALSE
      endif
  endif
else
  USE_RK64_PP = TRUE
  ifeq ($(USE_RK64_PP), TRUE)
  # Activates use of Hari explicit RK
  # will only work if USE_SUNDIALS_PP = FALSE
  TRANSPORT_CPP = TRUE
  endif
endif

ifeq ($(FUEGO_GAS), TRUE)
  Eos_Model       = Fuego
  Chemistry_Model = LiDryer
  Reactor_dir     = cvode
  Transport_Model = Simple
else
  Eos_Model       = GammaLaw
  Chemistry_Model = Null
  Reactor_dir     = null
  Transport_Model = Constant
endif

Blocs   := .

include $(PELE_ANALYSIS_HOME)/Tools/GNUmake/Make.Analysis
