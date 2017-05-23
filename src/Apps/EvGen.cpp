#include <iostream>
#include <fstream>
#include <vector>
#include <cstring>
#include <getopt.h>

#include "EventGenerator.h"
#include "FluxDriver.h"
#include "DecayRates.h"

void Usage(char* argv0);

int main(int argc, char** argv)
{
	const struct option longopts[] = 
	{
		{"smconfig", 	required_argument, 	0, 's'},
		{"detconfig", 	required_argument,	0, 'd'},
		{"fluxconfig", 	required_argument,	0, 'f'},
		{"root", 	required_argument,	0, 'r'},
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
	TFile *OutFile;
	std::ofstream Out;
	
	while((iarg = getopt_long(argc,argv, "s:d:f:r:o:h", longopts, &index)) != -1)
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
			case 'r':
				OutFile = new TFile(optarg, "RECREATE");
				break;
			case 'o':
				Out.open(optarg);
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

	EventGenerator * EvGen = new EventGenerator(SMConfig, DetConfig, FluxConfig);

	
	//TheFlux->SetBaseline(500);
	
	std::string Base = "Mass";
	std::string Name = Base;
	std::stringstream ssL;
	for (double Mass = 0; Mass < 0.5; Mass += 0.002)
	{
		ssL.str("");
		ssL.clear();
		ssL << Base << Mass;
		OutFile->mkdir(ssL.str().c_str());
		OutFile->cd(ssL.str().c_str());

		EvGen->SetMass(Mass);
		EvGen->MakeSterileFlux();

		EvGen->GetFluxDriverPtr()->GetTotal()->Write();
		Out << Mass << "\t";
		Out << EvGen->GetFluxDriverPtr()->GetTotal()->Integral("WIDTH") << "\t";
		Out << EvGen->GetFluxDriverPtr()->GetPion()->Integral("WIDTH") << "\t";
		Out << EvGen->GetFluxDriverPtr()->GetKaon()->Integral("WIDTH") << "\t";
		Out << EvGen->GetFluxDriverPtr()->GetKaon0()->Integral("WIDTH") << "\t";
		Out << EvGen->GetFluxDriverPtr()->GetMuon()->Integral("WIDTH") << std::endl;
	}

	OutFile->Close();

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
	std::cout <<"\n  -r,  --root" << std::endl;
	std::cout << "\t\tRoot output file" << std::endl;
	std::cout <<"\n  -h,  --help" << std::endl;
	std::cout << "\t\tPrint this message and exit" << std::endl;
}
