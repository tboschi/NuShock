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
	std::string DetConfig, FluxConfig_A, FluxConfig_P;
	std::string BaseName;
	
	while((iarg = getopt_long(argc,argv, "d:a:p:b:h", longopts, &index)) != -1)
	{
		switch(iarg)
		{
			case 'd':
				DetConfig.assign(optarg);
				break;
			case 'a':
				FluxConfig_A.assign(optarg);
				break;
			case 'p':
				FluxConfig_P.assign(optarg);
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

	FluxDriver * TheFlux_P = new FluxDriver(FluxConfig_P);
	FluxDriver * TheFlux_A = new FluxDriver(FluxConfig_A);
	Detector * TheBox = new Detector(DetConfig);

	Neutrino *N_up = new Neutrino(0, Neutrino::Dirac | Neutrino::Particle | Neutrino::Left );
	Neutrino *N_do = new Neutrino(0, Neutrino::Dirac | Neutrino::Particle | Neutrino::Right);
	
	std::stringstream ssL;
	double Mass;
	for (double Mass = 0.0; Mass < 2.0; Mass += 0.02)
	{
		N_up->SetMass(Mass);
		N_do->SetMass(Mass);

		ssL.str("");
		ssL.clear();
		ssL << BaseName << std::setfill('0') << std::setw(4) << Mass*1000;

		std::cout << "mass " << Mass << "\t in " << ssL.str() << std::endl;

		TheFlux_P->MakeFlux(N_up);
		TheFlux_A->MakeFlux(N_do);

		TheFlux_P->SetBaseline(574);
		TheFlux_A->SetBaseline(574);
		TheFlux_P->SetPOT(1e20);
		TheFlux_A->SetPOT(1e20);

		std::ofstream Out(ssL.str().c_str());

		double Start, End;
		double EnStep = TheFlux_P->RangeBin(Start, End);
		for (double Energy = Start+EnStep/2.0; Energy < End+EnStep/2.0; Energy += EnStep)
		{
			N_up->SetEnergyKin(Energy);
			N_do->SetEnergyKin(Energy);

			Out << Energy << "\t"; 
			Out << TheFlux_P->Intensity(N_up) << "\t"; 
			Out << TheFlux_A->Intensity(N_do) << "\t"; 
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
