.PHONY: clean include

INCDIR =	include
SRCDIR =	$(INCDIR)/src
APPDIR =	app
BINDIR =	bin
LIBDIR =	lib

ROOTLIB		= $(shell root-config --glibs)
ROOTCXX		= $(shell root-config --cflags)
#GENIELIB	= $(shell genie-config --libs)
#CUBALIB	= -L$(CUBA)/lib# -lcuba
#CUBACXX	= -I$(CUBA)/include 
#LHAPDFLIB	= -L$(LHAPDF)/lib -lLHAPDF
#LHAPDFCXX	= -I$(LHAPDF)/include

LDFLAGS  := -Wl,--no-as-needed $(LDFLAGS) $(ROOTLIB) $(GENIELIB) $(CUBALIB) $(LHAPDFLIB) -L$(LIBDIR)
#LDLIBS   := -lcuba
CXXFLAGS := $(CXXFLAGS) -fPIC -std=c++11 -O3 -mavx $(ROOTCXX) $(CUBACXX) $(LHAPDFCXX) -I$(INCDIR)

#apps and exctuables
CPP := $(shell find $(APPDIR) -maxdepth 1 -name '*.cpp')
SRC := $(shell find $(SRCDIR) -maxdepth 2 -name '*.cpp')

#main target
TGT := $(CPP:.cpp=)
DEP := $(SRC:.cpp=.o)

#dependencies of target
#INCDEP := $(HPP:%=$(INCDIR)/%/*.cpp)
#INCDEP := $(HPP:%=$(INCDIR)/%.cpp)
#DEP :=	$(patsubst %.cpp,%.o,$(wildcard $(INCDEP)))

all: $(TGT)
	@mkdir -p $(LIBDIR)
	@mkdir -p $(BINDIR)
	@echo "Cleaning up..."
	@cp $(DEP) $(LIBDIR)
	@mv $(TGT) $(BINDIR)
	@echo "Done!"

$(TGT): $(DEP)

include:
	$(eval DEP := $(shell find $(LIBDIR) -maxdepth 1 -name '*.o'))
	echo "dep " $(DEP)

clean:
	-find $(SRCDIR) -name "*.o" -delete
	-find $(INCDIR) -name "*~"  -delete
	-find $(LIBDIR) -mindepth 1 -name "*"   -delete
	-find $(APPDIR) -mindepth 1 -name "*~"  -delete
	-find $(BINDIR) -mindepth 1 -name "*"   -delete
