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
CXXFLAGS := $(CXXFLAGS) -fPIC -fopenmp -std=c++11 -O3 -mavx $(ROOTCXX) $(CUBACXX) $(LHAPDFCXX) -I$(INCDIR)

#apps and exctuables
CPP := $(shell find $(APPDIR) -maxdepth 1 -name '*.cpp')
SRC := $(shell find $(SRCDIR) -maxdepth 2 -name '*.cpp')

#main target
TARGET := $(CPP:.cpp=)
DEPEND := $(SRC:.cpp=.o)

all: $(TARGET)
	@mkdir -p $(LIBDIR)
	@mkdir -p $(BINDIR)
	@echo "Cleaning up..."
	@cp $(DEPEND) $(LIBDIR)
	@cp $(TARGET) $(BINDIR)
	@echo "Done!"

$(TARGET): $(DEPEND)

include:
	$(eval DEPEND := $(shell find $(LIBDIR) -maxdepth 1 -name '*.o'))
	echo "dep " $(DEPEND)

clean:
	-find $(SRCDIR) -name "*.o" -delete
	-find $(INCDIR) -name "*~"  -delete
	-find $(LIBDIR) -mindepth 1 -name "*"   -delete
	-find $(APPDIR) -mindepth 1 -name "*~"  -delete
	-find $(BINDIR) -mindepth 1 -name "*"   -delete
