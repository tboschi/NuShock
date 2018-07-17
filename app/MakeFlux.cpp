#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <cstring>
#include <getopt.h>

#include "Tools.h"
#include "Flux.h"
#include "Detector.h"
#include "Physics.h"

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
			default:
				break;
		}
	}

	//To have multiple output, handled by usage
	//std::ostream &Out = (OutFile.is_open()) ? OutFile : std::cout;

	Engine *TheEngine = new Engine(FluxConfig, 2, 2);	//creating 2FHC and 2RHC fluxedrivers

	//Detector * TheBox = new Detector(DetConfig);

	//neutrino left has same flux as antineutrino right
	//neutrino right has same flux as antineutrino left
	Neutrino *N_L = new Neutrino(0, Neutrino::Dirac | Neutrino::Left );
	Neutrino *N_R = new Neutrino(0, Neutrino::Dirac | Neutrino::Right);
	N_L->SetMixings(1, 1, 1);
	N_R->SetMixings(1, 1, 1);
	
	TheEngine->BindNeutrino(N_L, Engine::FHC, 0);
	TheEngine->BindNeutrino(N_R, Engine::FHC, 1);

	TheEngine->BindNeutrino(N_L, Engine::RHC, 0);	//Antinu R
	TheEngine->BindNeutrino(N_R, Engine::RHC, 1);	//Antinu L

	std::stringstream ssL;
	double Mass;
	for (double Mass = 0.0; Mass < 2.0; Mass += 0.05)
	{
		N_L->SetMass(Mass);
		N_R->SetMass(Mass);

		TheEngine->MakeFlux();

		TheEngine->ScaleBaseline(574);
		TheEngine->ScalePOT(1e20);

		ssL.str("");
		ssL.clear();
		ssL << BaseName << std::setfill('0') << std::setw(4) << Mass*1000 << ".dat";
		std::cout << "Mass " << Mass << "\t in " << ssL.str() << std::endl;

		std::ofstream Out(ssL.str().c_str());

		double Start, End;
		double EnStep = TheEngine->RangeWidth(Start, End);
		for (double Energy = Start; Energy < End; Energy += EnStep)
		{
			N_L->SetEnergyKin(Energy);
			N_R->SetEnergyKin(Energy);

			Out << Energy << "\t"; 
			Out << TheEngine->Intensity(Engine::FHC, 0) << "\t"; 
			Out << TheEngine->Intensity(Engine::FHC, 1) << "\t"; 
			Out << TheEngine->Intensity(Engine::RHC, 0) << "\t"; 
			Out << TheEngine->Intensity(Engine::RHC, 1) << "\t"; 
			Out << std::endl;
		}

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
