.PHONY: clean include

INCDIR =	include
SRCDIR =	$(INCDIR)/src
APPDIR =	app
BINDIR =	bin
LIBDIR =	lib

ROOTLIB		= $(shell root-config --glibs)	#libs for ROOT
ROOTINC		= $(shell root-config --cflags)	#libs for ROOT
#GENIELIB	= $(shell genie-config --libs)	#libs for GENIE
#GENIEINC	= $(shell genie-config --flags)	#libs for GENIE
CUBALIB		= -L$(CUBA)/lib		#libs for CUBA
CUBAINC		= -I$(CUBA)/include 	#libs for CUBA
LHAPDFLIB	= -L$(LHAPDF)/lib	#includes for LHAPDF
LHAPDFINC	= -I$(LHAPDF)/include	#includes for LHAPDF

LDFLAGS  := -Wl,--no-as-needed $(LDFLAGS) $(ROOTLIB) $(GENIELIB) $(CUBALIB) $(LHAPDFLIB) -L$(LIBDIR)
LDLIBS   := -lcuba -lLHAPDF
CXXFLAGS := $(CXXFLAGS) -fPIC -fopenmp -std=c++11 -O3 -mavx $(ROOTINC) $(CUBAINC) $(LHAPDFINC) -I$(INCDIR)

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
