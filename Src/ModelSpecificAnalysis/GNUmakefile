# PELE_ANALYSIS_HOME should be pointing to the home directory of your PeleAnalysis clone. 
# If the relative path does not compile, try using the absolute path.
PELE_ANALYSIS_HOME = ../..

# If you recursively cloned the repository, you should not need to set the paths for 
# PELE_PHYSICS_HOME and AMREX_HOME. If you chose to not recursively clone, or you want 
# to link PeleAnalysis to your versions of PELE_PHYSICS_HOME and AMREX_HOME used by 
# PeleLMeX or PeleC, make the following lines to point to the respective directories.
#
# PELE_PHYSICS_HOME = $(abspath ../../../PelePhysics)
# AMREX_HOME = $(abspath ../../../PelePhysics/Submodules/amrex)

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

#EBASE=plotQPD
#EBASE=plotTYtoLe
#EBASE=plotTransportCoeff
#EBASE=plotXtoY
#EBASE=plotYtoX
#EBASE=computeMixtureFraction
#EBASE=testQPDtools
#EBASE=testTsolve


CEXE_sources += $(EBASE).cpp

Eos_Model       = Fuego
Chemistry_Model = drm19
Transport_Model = Simple

Blocs   := .

include $(PELE_ANALYSIS_HOME)/Tools/GNUmake/Make.ModelSpecific
