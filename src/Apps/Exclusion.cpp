#include <iostream>
#include <fstream>
#include <cstring>
#include <getopt.h>

#include "EventGenerator.h"
#include "FluxDriver.h"
#include "DecayRates.h"
#include "Detector.h"

void Usage(char* argv0);

int main(int argc, char** argv)
{
	const struct option longopts[] = 
	{
		{"smconfig", 	required_argument, 	0, 's'},
		{"detconfig", 	required_argument,	0, 'd'},
		{"fluxconfig", 	required_argument,	0, 'f'},
		{"output", 	required_argument,	0, 'o'},
		{"help", 	no_argument,	 	0, 'h'},
		{0,	0, 	0,	0},
	};

	int index; 
	int iarg = 0;
	opterr = 1;
	
	//Initialize variables
	std::string SMConfig, DetConfig;
	std::string FluxConfig;
	std::ofstream OutFile;
	//TFile *OutFile;
	bool UeFlag = false;
	bool UmFlag = false;
	bool UtFlag = false;
	
	while((iarg = getopt_long(argc,argv, "s:d:f:o:EMTh", longopts, &index)) != -1)
	{
		switch(iarg)
		{
			case 's':
				SMConfig.assign(optarg);
				break;
			case 'd':
				DetConfig.assign(optarg);
				break;
			case 'f':
				FluxConfig.assign(optarg);
				break;
			case 'o':
				OutFile.open(optarg);
				//OutFile = new TFile(optarg, "RECREATE");
				break;
			case 'E':
				UeFlag = true;
				break;
			case 'M':
				UmFlag = true;
				break;
			case 'T':
				UtFlag = true;
				break;
			case 'h':
				Usage(argv[0]);
				return 1;
				return 1;
			default:
				break;
		}
	}

	//To have multiple output, handled by usage
//	std::ostream &Out = (OutFile.is_open()) ? OutFile : std::cout;

	std::cout << "M0" << std::endl;
	EventGenerator * EvGen = new EventGenerator(SMConfig, DetConfig, FluxConfig);
	std::cout << "M1" << std::endl;
	
	EvGen->SetChannel("EPI");

	EvGen->SetMass(0);
	
	EvGen->SetUe(1e-10);
	EvGen->SetUm(1e-10);
	EvGen->SetUt(1e-10);
	
	for (double Mass = 0.05; Mass < 0.5; Mass += 0.05)	//increase mass
	{
		EvGen->SetMass(Mass);

		for (double logUu = -5.0; logUu < -0.0; logUu += 0.1)	//increase Uu linearly
		{
			double Uu = pow(10.0, logUu);
			std::cout << "Uu " << Uu << std::endl << std::endl;
			if (UeFlag)
				EvGen->SetUe(Uu);
			if (UmFlag)
				EvGen->SetUm(Uu);
			if (UtFlag)
				EvGen->SetUt(Uu);

			EvGen->MakeSterileFlux();

			int Nevent = 0;
			
			/*
			for (int i = 0; i < 1000; ++i)	//Generate 1000 events
			{
				EvGen->SampleEnergy();	
				if (EvGen->EventInDetector()) ++Nevent;
			}*/

			EvGen->SampleEnergy();	
			OutFile << Mass << "\t" << Uu << "\t";
			OutFile	<< EvGen->EventProbability() << "\t";
			OutFile << EvGen->GetDecayPtr()->Total() << std::endl;
		}
	}

	return 0;
}
	
void Usage(char* argv0)
{
	std::cout << "Description" << std::endl;
	std::cout << "Usage : " << std::endl;
	std::cout << argv0 << " [OPTIONS]" << std::endl;
	std::cout <<"\n  -s,  --smconfig" << std::endl;
	std::cout << "\t\tStandard Model configuration file" << std::endl;
	std::cout <<"\n  -d,  --detconfig" << std::endl;
	std::cout << "\t\tDetector configuration file" << std::endl;
	std::cout <<"\n  -f,  --fluxconfig" << std::endl;
	std::cout << "\t\tFlux configuration file" << std::endl;
	std::cout <<"\n  -o,  --output" << std::endl;
	std::cout << "\t\tOutput file" << std::endl;
	std::cout <<"\n  -E,  -M,  -T" << std::endl;
	std::cout << "\t\tSelect which mixing element" << std::endl;
	std::cout <<"\n  -h,  --help" << std::endl;
	std::cout << "\t\tPrint this message and exit" << std::endl;
}
