.PHONY: clean

INCDIR =	include
APPDIR =	app
BINDIR =	bin
LIBDIR =	lib

ROOTLIB		= $(shell root-config --glibs)
ROOTCXX		= $(shell root-config --cflags)
#GENIELIB	= $(shell genie-config --libs)
CUBALIB		= -L$(CUBA)/lib# -lcuba
CUBACXX		= -I$(CUBA)/include 
#LHAPDFLIB	= -L$(LHAPDF)/lib -lLHAPDF
#LHAPDFCXX	= -I$(LHAPDF)/include

LDFLAGS  := -Wl,--no-as-needed $(LDFLAGS) $(ROOTLIB) $(GENIELIB) $(CUBALIB) $(LHAPDFLIB) -L$(LIBDIR)
LDLIBS   := -lcuba
CXXFLAGS := $(CXXFLAGS) -fPIC -std=c++11 -O3 -mavx $(ROOTCXX) $(CUBACXX) $(LHAPDFCXX) -I$(INCDIR)

#apps and exctuables
CPP =	ProductionScale	\
	DecayBranch	\
	MakeFlux
	#MakeFlux	\
	Eps2Dat		\
	Eps2Root	

#header folders
HPP =	Tools/Const		\
	Tools/asa047		\
	Tools/Integration	\
	Tools/Particle		\
	Physics/Amplitude	\
	Physics/DecayRates	\
	Physics/Production	\
	Physics/PhaseSpace	\
	Physics/Neutrino	\
	Flux/Flux		\
	Flux/FluxDriver		\
	Detector/Detector	\
	Detector/Efficiency	\
	Event/EventGenerator	\

#main target
TGT :=	$(CPP:%=$(APPDIR)/%)

#dependencies of target
#INCDEP := $(HPP:%=$(INCDIR)/%/*.cpp)
INCDEP := $(HPP:%=$(INCDIR)/%.cpp)
DEP :=	$(patsubst %.cpp,%.o,$(wildcard $(INCDEP)))

all: $(TGT)
	@mkdir -p $(BINDIR)
	@mkdir -p $(LIBDIR)
	@echo "Moving stuff..."
	@mv $(TGT) $(BINDIR)
	@cp $(DEP) $(LIBDIR)
	@echo "Done!"

$(TGT): $(DEP)

test: app/TestGenPS
	@mv app/TestGenPS $(BINDIR)
	
dep:
	echo $(INCDEP)

clean:
	find $(INCDIR) -mindepth 1 -name "*.o" -delete
	find $(INCDIR) -mindepth 1 -name "*~"  -delete
	find $(APPDIR) -mindepth 1 -name "*~"  -delete
	find $(BINDIR) -mindepth 1 -name "*"   -delete
