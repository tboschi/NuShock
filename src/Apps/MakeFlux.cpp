#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <cstring>
#include <getopt.h>

#include "EventGenerator.h"
#include "FluxDriver.h"
#include "DecayRates.h"

#include "TFile.h"
#include "TH1D.h"

void Usage(char* argv0);

int main(int argc, char** argv)
{
	const struct option longopts[] = 
	{
		{"detconfig", 	required_argument,	0, 'd'},
		{"fluxconfig", 	required_argument,	0, 'f'},
		{"basename", 	required_argument,	0, 'b'},
		{"help", 	no_argument,	 	0, 'h'},
		{0,	0, 	0,	0},
	};

	int index; 
	int iarg = 0;
	opterr = 1;
	
	//Initialize variables
	std::string DetConfig, FluxConfig;
	std::string BaseName;
	
	while((iarg = getopt_long(argc,argv, "d:f:b:h", longopts, &index)) != -1)
	{
		switch(iarg)
		{
			case 'd':
				DetConfig.assign(optarg);
				break;
			case 'f':
				FluxConfig.assign(optarg);
				break;
			case 'b':
				BaseName.assign(optarg);
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

	FluxDriver * TheFlux = new FluxDriver(FluxConfig);
	Detector * TheBox = new Detector(DetConfig);
	
	std::stringstream ssL;
	double Mass;
	for (unsigned int iMass = 0; iMass < 500; ++iMass)
	{
		Mass = iMass/1000.0;
		ssL.str("");
		ssL.clear();
		ssL << BaseName << std::setfill('0') << std::setw(3) << iMass;

		std::cout << "m " << Mass << "\t in " << ssL.str() << std::endl;
		TheFlux->MakeFlux(Mass);
		TheFlux->SetBaseline(TheBox->GetElement("Baseline"));
		TheFlux->SetPOT(TheBox->GetElement("POT/s"));

		std::ofstream Out(ssL.str().c_str());

		double EnStep = (TheFlux->GetRangeEnd() - TheFlux->GetRangeStart()) / (10.0*TheFlux->GetBinNumber());
		for (double E = TheFlux->GetRangeStart(); E < TheFlux->GetRangeEnd(); E += EnStep)
			Out << E + EnStep/2.0 << "\t" << TheFlux->GetIntensityNeut(E) + TheFlux->GetIntensityAnti(E) << std::endl;

		Out.close();
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
	std::cout <<"\n  -r,  --root" << std::endl;
	std::cout << "\t\tRoot output file" << std::endl;
	std::cout <<"\n  -h,  --help" << std::endl;
	std::cout << "\t\tPrint this message and exit" << std::endl;
}
