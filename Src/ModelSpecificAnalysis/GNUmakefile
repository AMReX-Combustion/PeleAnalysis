PELE_ANALYSIS_HOME ?= /scratch/bsoriano/PeleAnalysis
PELE_PHYSICS_HOME ?= /scratch/bsoriano/PeleC/Submodules/PelePhysics/
AMREX_HOME ?= /scratch/bsoriano/PeleC/Submodules/AMReX/

PRECISION  = DOUBLE
PROFILE    = FALSE
DEBUG      = FALSE
DIM        = 3
COMP	   = intel
USE_MPI    = TRUE
USE_OMP    = FALSE
NEEDS_f90_SRC = FALSE

# EBASE = bruno
# EBASE = newbruno
# EBASE = regridPlt_new
EBASE = reduce_size

USE_SUNDIALS_PP = TRUE

CEXE_sources += $(EBASE).cpp
include $(AMREX_HOME)/Tools/GNUMake/Make.defs

ifeq ($(NEEDS_f90_SRC),TRUE)
  f90EXE_sources += $(EBASE)_nd.f90
endif


# This sets the EOS directory in $(PELE_PHYSICS_HOME)/Eos
#Eos_dir     := Soave-Redlich-Kwong
Eos_dir     := Fuego

# This sets the network directory in $(PELE_PHYSICS_HOME)/Reactions
Eos_Model := Fuego
#Reactions_dir := Null

# This sets the model directory in $(PELE_PHYSICS_HOME)/Support/Fuego/Mechanism/Models
Chemistry_Model := TPRF_98species_gpu
# Chemistry_Model := grimech30_nasa9
# Chemistry_Model := LiDryer

# This sets the transport directory in $(PELE_PHYSICS_HOME)/Transport
Transport_dir := Simple

include $(PELE_ANALYSIS_HOME)/Tools/GNUmake/Make.Analysis
