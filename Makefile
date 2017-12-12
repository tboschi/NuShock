#
# Makefile
#
# Author: Tommaso Boschi
#

.PHONY: clean

#SHELL = /bin/sh
#NAME = all
#MAKEFILE = Makefile
SRC_DIR = src

# Include machine specific flags and locations (inc. files & libs)
#
include $(GENIE)/src/make/Make.include


#Main executable to be compiled
NEW =	Exclusion	\
	ExcluMix	\
	MegaExcl	\
	MegaPlot	\
	Timing		\
	Simulation	\
	Probability	\
	DecayPlot	\
	#Plotter	\
	MakeFlux	\
	MakeEfficiency	\
	CrossExcl	\
	EnergyBack	\
	#RegionSearch	\
	SearchRegion	\
	FeldCous	\
	Plotter		\
	Coupler		\
	CosmoBounds
	#EFTgamma	\
	EFTexcl		\
	GammaRequired
	#XSec		\
	#PSscatter	\
	Rate		\
	GenBack		\
	Eps2Dat		\
	Kine		\
	CLs		\
	Width		\
	#Eps2Root	\
	#FakeElectron	\
	PionMuonFlux	\
	Kine		\
	DUNE_FGT

TGT :=  $(NEW:%=$(SRC_DIR)/Apps/%)
BIN :=  $(NEW:%=bin/%)

#Dependencies of the Main
DEP =	Tools/Tools		\
	DecayRates/DecayRates	\
	DecayRates/3Body	\
	SterileFlux/Flux	\
	SterileFlux/FluxDriver	\
	EventGenerator/EventGenerator	\
	Detector/Detector	\
	Detector/Efficiency	\
	Particle/Particle	\
	Background/Background	\
	Scattering/Nucleon	\

INC_DIR := $(patsubst %,-I$(SRC_DIR)/%,$(subst /, ,$(DEP)))
DEP :=  $(DEP:%=$(SRC_DIR)/%.o)

GENIE_LIBS  = $(shell $(GENIE)/src/scripts/setup/genie-config --libs)
LDFLAGS  := $(LDFLAGS) $(GENIE_LIBS) $(LIBRARIES) $(CERN_LIBRARIES) -L$(CUBA)/lib
LDLIBS   := $(LDLIBS) -lcuba
CXXFLAGS := $(CXXFLAGS) -I$(CUBA)/include $(INCLUDES) $(INC_DIR) 

#Using implicit rules for G++ and linker
#Then moves exec into main folder
all: $(TGT)
	mv $(TGT) ./bin

#xsec: $(DEP)
#	$(CC) $(CXXFLAGS) $(INCLUDES) $(SRC_DIR)/Apps/XSec.cpp -o XSec $(LDFLAGS)
#	mv XSec ./bin

$(TGT): $(DEP)

wow:
	@echo $(GENIE_LIBS)
	@echo $(LIBRARIES)
	@echo $(CERN_LIBRARIES)
	@echo $(LDFLAGS)
	@echo $(LDLIBS)
	@echo $(CXXFLAGS)

#################### CLEANING

clean: 
	find $(SRC_DIR) -name "*.o" -delete
	find $(SRC_DIR) -name "*~" -delete
	find $(SRC_DIR) -name "core" -delete
	$(RM) $(BIN)
