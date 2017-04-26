#include <iostream>
#include <fstream>
#include <vector>
#include <cstring>
#include <getopt.h>

#include "Tools.h"
#include "DecayRates.h"

int main(int argc, char** argv)
{
	const struct option longopts[] = 
	{
		{"output", 	required_argument, 	0, 'f'},
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
	double M_Sterile = 0.350;	//350 MeV
	double U_e = 1.0/sqrt(3.0);	//All mixing enabled to be maximal
	double U_m = 1.0/sqrt(3.0);
	double U_t = 1.0/sqrt(3.0);
	
	while((iarg = getopt_long(argc,argv, "f:s:e:m:t:h", longopts, &index)) != -1)	
	{
		switch(iarg)
		{
			case 'f':
				OutFile.open(optarg);
				if (!OutFile.is_open())
				{
					std::cout << "Wrong input!" << std::endl;
					return 1;
				}
				break;
			case 'e':
				U_e = strtod(optarg,NULL);
				break;
			case 'm':
				U_m = strtod(optarg,NULL);
				break;
			case 't':
				U_t = strtod(optarg,NULL);
				break;
			case 'h':
				std::cout << "Compute decay spectrum from 0 MeV to 500 MeV" << std::endl;
				std::cout << "Usage : " << std::endl;
				std::cout << "decayplot [OPTIONS]" << std::endl;
				std::cout <<"\n  -o,  --output" << std::endl;
				std::cout << "\t\tOutput file to save plot; if not specified the standard output is used instead" << std::endl;
				std::cout <<"\n  -e,  --ue4" << std::endl;
				std::cout << "\t\tSpecifiy electron-heavy mixing, default 1/sqrt(3)" << std::endl;
				std::cout <<"\n  -m,  --um4" << std::endl;
				std::cout << "\t\tSpecifiy muon-heavy mixing, default 1/sqrt(3)" << std::endl;
				std::cout <<"\n  -t,  --ut4" << std::endl;
				std::cout << "\t\tSpecifiy tau-heavy mixing, default 1/sqrt(3)" << std::endl;
				std::cout <<"\n  -h,  --help" << std::endl;
				std::cout << "\t\tPrint this message and exit" << std::endl;
				return 1;
			default:
				break;
		}
	}

	std::ostream &Out = (OutFile.is_open()) ? OutFile : std::cout;

	Decay * SuperGamma = new Decay(M_Sterile, U_e, U_m, U_t);
	std::vector<std::string> vChannel = SuperGamma->ListChannels();

	Out << "#1Mass";
	for (int i = 0; i < vChannel.size(); ++i)
		Out << i+2 << vChannel.at(i) << "\t";
	Out << std::endl;

	
	for (M_Sterile = 0; M_Sterile < 0.500; M_Sterile += 0.005)
	{
		double Sum = 0;
		SuperGamma->SetMass(M_Sterile);
		for (int i = 1; i < vChannel.size(); ++i)
		{
			Sum += SuperGamma->Branch(vChannel.at(i));
			Out << Sum << "\t";
		}
		Out << std::endl;
	}

	delete SuperGamma;

	return 0;
}
