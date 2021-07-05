APPDIR = app
TESTDIR = test
INCDIR = include
SRCDIR = src
BINDIR = bin

CUBA ?= cuba
LHAPDF ?= lhapdf


ROOTLIB		= $(shell root-config --libdir)	#libs for ROOT
ROOTINC		= $(shell root-config --cflags)	#libs for ROOT
ROOTLD		= $(shell root-config --glibs) #libs for ROOT
CUBALIB		= -L$(CUBA)/lib	#libs for CUBA
CUBAINC		= -I$(CUBA)/include 	#libs for CUBA
LHAPDFLIB	= -L$(LHAPDF)/lib	#includes for LHAPDF
LHAPDFINC	= -I$(LHAPDF)/include	#includes for LHAPDF
FCCLINC		= -Ifccl/include	# include for fccl


#optimization
ARCH ?= -march=native

LDFLAGS  := $(LDFLAGS) $(CUBALIB) $(LHAPDFLIB)
LDLIBS   := -Wl,--as-needed -lcuba -lLHAPDF $(ROOTLD) -lTMVA -lTMVAGui
CXXFLAGS := -Wall -fPIC -std=c++11 -O2 $(ARCH) $(ROOTINC) $(CUBAINC) $(LHAPDFINC) $(FCCLINC) -I$(INCDIR)


#apps and exctuables
TARGETS := $(wildcard $(APPDIR)/*.cpp)
TESTS   := $(wildcard $(TESTDIR)/*.cpp)
SOURCES := $(wildcard $(SRCDIR)/*.cpp)
HEADERS := $(wildcard $(INCDIR)/*/*.h*)
#SOURCES := $(SRCDIR)/Particle.cpp $(SRCDIR)/Track.cpp $(SRCDIR)/OpenQQ.cpp
#HEADERS := $(INCDIR)/physics/Particle.h $(INCDIR)/physics/Track.h $(INCDIR)/physics/OpenQQ.h

DEPENDS := $(SOURCES:.cpp=.d)
OBJECTS := $(SOURCES:.cpp=.o)
TARGETS := $(if $(APP), $(APPDIR)/$(APP), $(TARGETS:.cpp=))
TESTS   := $(if $(TEST), $(TESTDIR)/$(TEST), $(TESTS:.cpp=))


all: $(TARGETS)
	@mkdir -p $(BINDIR)
	@cp $^ $(BINDIR)
	@echo "Done!"

test: $(TESTS)
	@echo "Done!"

help:
	@echo Targets found are $(TARGETS)
	@echo Sources found are $(SOURCES)
	@echo Headers found are $(HEADERS)
	@echo "If you need to build just one file, do make APP=name"
	@echo "or if you need to specify an architecture, do make ARCH=arch"
	@echo "Enjoy your compilation"

$(TARGETS): $(OBJECTS)

$(TESTS): $(OBJECTS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -MMD -MP -c $< -o $@

-include $(DEPENDS)

#$(OBJECT): $(HEADERS)


clean:
	@$(RM) $(TARGETS)
	@$(RM) $(TESTS)
	@$(RM) $(OBJECTS)
	@$(RM) $(DEPENDS)
	@$(RM) -r $(BINDIR)

.PHONY: all test help clean
