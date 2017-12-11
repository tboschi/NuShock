#include <iostream>
#include <fstream>
#include <cstring>
#include <getopt.h>

#include "EventGenerator.h"
#include "FluxDriver.h"
#include "DecayRates.h"
#include "Detector.h"

#include "TH1D.h"

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
	std::string SMConfig, DetConfig, FluxConfig;
	std::string Channel;
	std::ofstream OutFile;
	char Flag;
	
	while((iarg = getopt_long(argc,argv, "s:d:f:c:o:U:h", longopts, &index)) != -1)
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
			case 'c':
				Channel.assign(optarg);
				break;
			case 'o':
				OutFile.open(optarg);
				break;
			case 'U':
				Flag = *optarg;
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
	std::ostream &Out = (OutFile.is_open()) ? OutFile : std::cout;

	EventGenerator * EvGenLow = new EventGenerator(SMConfig, DetConfig, FluxConfig);
	EventGenerator * EvGenHig = new EventGenerator(SMConfig, DetConfig, FluxConfig);
	
	EvGenLow->SetChannel(Channel);
	EvGenHig->SetChannel(Channel);

	EvGenLow->SetEnergy(1);
	EvGenHig->SetEnergy(19);

	double Uu;
	if (Flag == 'E')
		Uu = EvGenLow->GetUe();
	if (Flag == 'M')
		Uu = EvGenLow->GetUm();
	if (Flag == 'T')
		Uu = EvGenLow->GetUt();
	
	for (double Mass = 0.0; Mass < 0.5; Mass += 0.001)	//increase mass
	{
		EvGenLow->SetMass(Mass);
		EvGenHig->SetMass(Mass);

		Out << Mass << "\t" << EvGenLow->DecayProb()/Uu/Uu << "\t" << EvGenHig->DecayProb()/Uu/Uu << std::endl;
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
	std::cout <<"\n  -h,  --help" << std::endl;
	std::cout << "\t\tPrint this message and exit" << std::endl;
}
