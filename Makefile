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
NEW =	DecayPlot
TGT :=  $(NEW:%=$(SRC_DIR)/Apps/%)

#Dependencies of the Main
DEP =	Tools/Tools		\
	DecayRates/DecayRates

INC_DIR := $(patsubst %,-I$(SRC_DIR)/%,$(subst /, ,$(DEP)))
DEP :=  $(DEP:%=$(SRC_DIR)/%.o)

GENIE_LIBS  = $(shell $(GENIE)/src/scripts/setup/genie-config --libs)
LDFLAGS  := $(LDFLAGS) $(GENIE_LIBS) $(LIBRARIES) $(CERN_LIBRARIES)
CXXFLAGS := $(CXXFLAGS) $(INCLUDES) $(INC_DIR)

#Using implicit rules for G++ and linker
#Then moves exec into main folder
all: $(TGT)
	mv $(TGT) ./

$(TGT): $(DEP)


#################### CLEANING

clean: 
	$(RM) $(TGT)
	$(RM) $(DEP)
	$(RM) *.o *~ core 
	$(RM) $(NEW)
