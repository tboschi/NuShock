#include <iostream>
#include <fstream>
#include <vector>
#include <cstring>
#include <getopt.h>
#include <iomanip>

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

	std::string base;
	
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
				base.assign(optarg);
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
	//std::ostream &Out = (OutFile.is_open()) ? OutFile : std::cout;

	EventGenerator * EvGen = new EventGenerator(SMConfig, DetConfig, FluxConfig);

	
	//TheFlux->SetBaseline(500);
	

	for (int j = 0; j < 500; ++j)
	{
		double Mass = j*0.001;
		std::stringstream ssL;
		ssL << base << std::setw(3) << std::setfill('0') << j;

		std::cout << ssL.str() << std::endl;

		Out.open(ssL.str().c_str());

		EvGen->SetMass(Mass);
		EvGen->MakeSterileFlux();

		double E;
		for (int bin = 0; bin <= 100; ++bin)
		{
			E = bin*0.2;
			Out << E << "\t";
			if (EvGen->GetFluxDriverPtr()->GetTotal()->GetBinContent(bin) == 0)
			       Out << 1e-20 << "\t";
			else Out << EvGen->GetFluxDriverPtr()->GetTotal()->GetBinContent(bin) << "\t";
			//Out << EvGen->GetFluxDriverPtr()->GetPion()->GetBinContent(bin) << "\t";
			//Out << EvGen->GetFluxDriverPtr()->GetKaon()->GetBinContent(bin) << "\t";
			//Out << EvGen->GetFluxDriverPtr()->GetKaon0()->GetBinContent(bin) << "\t";
			//Out << EvGen->GetFluxDriverPtr()->GetMuon()->GetBinContent(bin) << std::endl;
			Out << std::endl;
		}

		Out.close();
	}

	//OutFile->Close();

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
