#include <iostream>
#include <fstream>
#include <vector>
#include <cstring>
#include <getopt.h>

#include "Tools.h"
#include "EventGenerator.h"
#include "3Body.h"

int main(int argc, char** argv)
{
	const struct option longopts[] = 
	{
		{"output", 	required_argument, 	0, 'o'},
		{"sterile", 	required_argument, 	0, 's'},
		{"ue4", 	required_argument, 	0, 'e'},
		{"um4", 	required_argument, 	0, 'm'},
		{"ut4", 	required_argument, 	0, 't'},
		{"help", 	no_argument,		0, 'h'},
		{0,	0, 	0,	0},
	};

	int index; 
	int iarg = 0;
	opterr = 1;
	
	std::ofstream OutFile;
	std::string SMConfig;

	while((iarg = getopt_long(argc,argv, "o:s:h", longopts, &index)) != -1)	
	{
		switch(iarg)
		{
			case 'o':
				OutFile.open(optarg);
				if (!OutFile.is_open())
				{
					std::cout << "Wrong input!" << std::endl;
					return 1;
				}
				break;
			case 's':
				SMConfig.assign(optarg);
				break;
			case 'h':
				std::cout << "Compute decay spectrum from 0 MeV to 500 MeV" << std::endl;
				std::cout << "Usage : " << std::endl;
				std::cout << "decayplot [OPTIONS]" << std::endl;
				std::cout <<"\n  -o,  --output" << std::endl;
				std::cout << "\t\tOutput file to save plot; if not specified the standard output is used instead" << std::endl;
				std::cout <<"\n  -s,  --smconfig" << std::endl;
				std::cout << "\t\tSMconfig" << std::endl;
				std::cout <<"\n  -h,  --help" << std::endl;
				std::cout << "\t\tPrint this message and exit" << std::endl;
				return 1;
			default:
				break;
		}
	}

	double Ms, Ue, Um, Ut;
	
	std::string Line, Key;
	std::stringstream ssL;
	double Element;

	std::ifstream ConfigFile(SMConfig.c_str());
	while (std::getline(ConfigFile, Line))
	{
		if (Line[0] == '#') continue;

		ssL.str("");
		ssL.clear();
		ssL << Line;

		ssL >> Key >> Element;
		if (Key == "M_Sterile") Ms = Element;
		if (Key == "U_e") Ue = Element;
		if (Key == "U_m") Um = Element;
		if (Key == "U_t") Ut = Element;
	}
	ConfigFile.close();

	std::ostream &Out = (OutFile.is_open()) ? OutFile : std::cout;

	Decay * SuperGamma = new Decay(Ms, Ue, Um, Ut);

	std::vector<std::string> vChannel(SuperGamma->ListChannels());
	
	Out << "#Mass\tTotal\t";
	for (int i = 1; i < vChannel.size(); ++i)
		Out << vChannel.at(i) << "\t";
	Out << std::endl;

	for (Ms = 0; Ms < 0.5; Ms += 0.001)
	{
		SuperGamma->SetMass(Ms);

		Out << Ms << "\t" << SuperGamma->Total() << "\t";
		for (int i = 1; i < vChannel.size(); ++i)
			Out << SuperGamma->Branch(vChannel.at(i)) << "\t";
		Out << std::endl;
	}
	

	return 0;
}
