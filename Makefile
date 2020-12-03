APPDIR = app
INCDIR = include
SRCDIR = src
BINDIR = bin

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
TARGETS := $(wildcard $(APPDIR)/*.cpp)
SOURCES := $(wildcard $(SRCDIR)/*.cpp)
HEADERS := $(wildcard $(INCDIR)/*/*.h*)

DEPENDS := $(SOURCES:.cpp=.d)
OBJECTS := $(SOURCES:.cpp=.o)
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

.PHONY: all help clean
