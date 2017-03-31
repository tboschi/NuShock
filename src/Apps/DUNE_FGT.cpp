#include <iostream>
#include <fstream>
#include <vector>
#include <cstring>
#include <getopt.h>

#include "Detector.h"

int main(int argc, char** argv)
{

	const struct option longopts[] = 
	{
		{"config", 	required_argument,	0, 'c'},
		{"output", 	required_argument,	0, 'o'},
		{"efficiency", 	required_argument,	0, 'e'},
		{"help",	no_argument,		0, 'h'},
		{0,	0, 	0,	0},
	};

	int index; 
	int iarg = 0;
	opterr = 1;
	
	std::ofstream OutFile;
	
	//Initialize variables for getopt
	std::string InFile;
	double EnEff = 1.0;

	while((iarg = getopt_long(argc,argv, "c:o:e:h", longopts, &index)) != -1)	
	{
		switch(iarg)
		{
			case 'c':
				InFile.assign(optarg);
				break;
			case 'o':
				OutFile.open(optarg);
				break;
			case 'e':
				EnEff = strtod(optarg, NULL);
				break;
			case 'h':
				std::cout << "Test for detector class: prints properties and efficiency at given energy (GeV)" << std::endl;
				std::cout << "Usage : " << std::endl;
				std::cout << "DUNE_FGT [OPTIONS]" << std::endl;
				std::cout << "\n  -c,  --config" << std::endl;
				std::cout << "\t\tSpecify the configuration file for the detector" << std::endl;
				std::cout << "\n  -o,  --output" << std::endl;
				std::cout << "\t\tSave output tu file instead of standard output" << std::endl;
				std::cout << "\n  -e,  --efficiency" << std::endl;
				std::cout << "\t\tEnergy at which the efficiency is calculated, default 1.0 GeV" << std::endl;
				std::cout << "\n  -h,  --help" << std::endl;
				std::cout << "\t\tPrint this message and exit" << std::endl;
				return 1;
			default:
				std::cout << "Type DUNE_FGT --help or -h to see usage" << std::endl; 
				return 1;
		}
	
	}

	//To have multiple output, handled by usage
	std::ostream &Out = (OutFile.is_open()) ? OutFile : std::cout;

	//Create object from classes
	Detector * DUNE = new Detector(InFile);
	
	//Main body
	std::vector<std::string> ListOfKey = DUNE->ListKey();
	for (int i = 0; i < ListOfKey.size(); ++i)
	{
		Out << i << "\t" << ListOfKey.at(i) << "\t";
		Out << DUNE->GetElement(ListOfKey.at(i)) << std::endl;
	}

	Out << "Efficiency at " << EnEff << ": " << DUNE->Efficiency("MUPI", EnEff) << std::endl;

	//Garbage collection
	delete DUNE;

	return 0;
}
	
