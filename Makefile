INCDIR = include
APPDIR = app
BINDIR = bin
DOCDIR = doc

## root
ROOTLIB	= $(shell root-config --glibs)
ROOTCXX	= $(shell root-config --cflags)

ROOTLIB		= $(shell root-config --glibs)	#libs for ROOT
ROOTCXX		= $(shell root-config --cflags)	#libs for ROOT
CUBALIB		= -L$(CUBA)/lib		#libs for CUBA
CUBACXX		= -I$(CUBA)/include 	#libs for CUBA
LHAPDFLIB	= -L$(LHAPDF)/lib	#includes for LHAPDF
LHAPDFCXX	= -I$(LHAPDF)/include	#includes for LHAPDF

# no need to link genie anymore
#GENIELIB	= $(shell genie-config --libs)	#libs for GENIE
#GENIEINC	= $(shell genie-config --flags)	#libs for GENIE


#optimization
ARCH ?= -march=native

LDFLAGS  := -Wl,--no-as-needed $(LDFLAGS) $(ROOTLIB) $(GENIELIB) $(CUBALIB) $(LHAPDFLIB)
LDLIBS   := -lcuba -lLHAPDF
CXXFLAGS := $(CXXFLAGS) -fPIC -fopenmp -std=c++11 -O3 -mavx
CXXFLAGS := $(DEBUG) $(WARNING) -fPIC -std=c++11 -O3 $(ARCH)  $(ROOTCXX) $(CUBACXX) $(LHAPDFCXX) -I$(INCDIR)


#apps and exctuables
TARGETS := $(shell find $(APPDIR) -maxdepth 1 -name '*.cpp')
SOURCES := $(shell find $(INCDIR) -maxdepth 2 -name '*.cpp')
HEADERS := $(shell find $(INCDIR) -maxdepth 2 -name '*.h')

OBJECTS := $(SOURCES:.cpp=.o)
DEPENDS := $(SOURCES:.cpp=.d)
TARGETS := $(if $(APP), $(APPDIR)/$(APP), $(TARGETS:.cpp=))


all: $(TARGETS)
	@mkdir -p $(BINDIR)
	@cp $^ $(BINDIR)
	@echo "Done!"

help:
	@echo Targets found are $(TARGETS)
	@echo Sources found are $(SOURCES)
	@echo Headers found are $(HEADERS)
	@echo "If you need to build just one file, do make APP=name"
	@echo "or if you need to specify an architecture, do make ARCH=arch"
	@echo "To build documentation, make doc"
	@echo "Enjoy your compilation"

$(TARGETS): $(OBJECTS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -MMD -MP -c $< -o $@

-include $(DEPENDS)

#$(OBJECT): $(HEADERS)


clean:
	$(RM) $(TARGETS)
	$(RM) $(OBJECTS)
	$(RM) $(DEPENDS)
	$(RM) -r $(BINDIR)

	$(MAKE) -C $(DOCDIR) clean


.PHONY: all doc help clean
