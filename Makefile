CXX = g++




SRC = src
LIB = lib

CPPFLAGS = 
CXXFLAGS = `root-config --cflags --glibs` -fPIC


BUILD = $(SRC)/Apps/DecayPlot.cpp

all:
	$(CXX) $(SRC)/Tools/Tools.cpp $(SRC)/Apps/DecayPlot.cpp -I $(SRC)/Tools -I$(GENIE)/src $(CXXFLAGS) -o test
