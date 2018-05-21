.PHONY: clean

INCDIR =	include
APPDIR =	app
BINDIR =	bin
LIBDIR =	lib

ROOTLIB		= $(shell root-config --glibs)
ROOTCXX		= $(shell root-config --cflags)
GENIELIB	= $(shell genie-config --libs)
#CUBALIB		= -L$(CUBA)/lib -lcuba
CUBACXX		= -I$(CUBA)/include 
#LHAPDFLIB	= -L$(LHAPDF)/lib -lLHAPDF
LHAPDFCXX	= -I$(LHAPDF)/include

LDFLAGS  := -Wl,--no-as-needed $(LDFLAGS) $(ROOTLIB) $(GENIELIB) $(CUBALIB) $(LHAPDFLIB) -L$(LIBDIR)
CXXFLAGS := $(CXXFLAGS) -std=c++11 -O3 -mavx $(ROOTCXX) $(CUBACXX) $(LHAPDFCXX) -I$(INCDIR)

#apps and exctuables
CPP =	MakeFlux	\
	Eps2Dat		\
	Eps2Root	

#header folders
HPP =	Tools		\
	Flux		\
	Pysics		\
	Detector	\
	Event		\

#main target
TGT :=	$(CPP:%=$(APPDIR)/%)

#dependencies of target
INCDEP := $(HPP:%=$(INCDIR)/%/*.cpp)
DEP :=	$(patsubst %.cpp,%.o,$(wildcard $(INCDEP)))

all: $(TGT)
	@mkdir -p $(BINDIR)
	@mkdir -p $(LIBDIR)
	@echo "Moving stuff..."
	@mv $(TGT) $(BINDIR)
	@cp $(DEP) $(LIBDIR)
	@echo "Done!"

$(TGT): $(DEP)

clean:
	find $(INCDIR) -mindepth 1 -name "*.o" -delete
	find $(INCDIR) -mindepth 1 -name "*~"  -delete
	find $(APPDIR) -mindepth 1 -name "*~"  -delete
	find $(BINDIR) -mindepth 1 -name "*"   -delete
