PELE_ANALYSIS_HOME ?= $(HOME)/src/CCSE/Pele/PeleAnalysis
PELE_PHYSICS_HOME ?= $(HOME)/src/CCSE/Pele/PelePhysics

PRECISION  = DOUBLE
PROFILE    = FALSE
DEBUG      = FALSE

#DIM        = 1
DIM        = 2
#DIM        = 3

COMP	   = gcc
USE_MPI    = FALSE

EBASE = plotYtoX
EBASE = plotXtoY

# This sets the EOS directory in $(PELE_PHYSICS_HOME)/Eos
#Eos_dir     := Soave-Redlich-Kwong
Eos_dir     := Fuego

# This sets the network directory in $(PELE_PHYSICS_HOME)/Reactions
Reactions_dir := Fuego

# This sets the model directory in $(PELE_PHYSICS_HOME)/Support/Fuego/Mechanism/Models
#Chemistry_Model := Davis
#Chemistry_Model := LiDryer
Chemistry_Model := drm19

# This sets the transport directory in $(PELE_PHYSICS_HOME)/Transport
Transport_dir := EGLib
#Transport_dir := Simple

include $(PELE_ANALYSIS_HOME)/Tools/GNUmake/Make.Analysis
