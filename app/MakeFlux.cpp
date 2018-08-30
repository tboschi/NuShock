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
		{"dirac", 	no_argument,		0, 'D'},
		{"majorana", 	no_argument,		0, 'M'},
		{"help", 	no_argument,	 	0, 'h'},
		{0,	0, 	0,	0},
	};

	int index; 
	int iarg = 0;
	opterr = 1;
	
	//Initialize variables
	std::string DetConfig, FluxConfig, Channel = "ALL", BaseName;
	bool Dirac = true;	//default unpolarised
	
	while((iarg = getopt_long(argc,argv, "c:d:f:b:DMh", longopts, &index)) != -1)
	{
		switch(iarg)
		{
			case 'c':
				Channel.assign(optarg);
				break;
			case 'd':
				DetConfig.assign(optarg);
				break;
			case 'f':
				FluxConfig.assign(optarg);
				break;
			case 'b':
				BaseName.assign(optarg);
				break;
			case 'D':
				Dirac = true;
				break;
			case 'M':
				Dirac = false;
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

	Tracker *TheBox = new Tracker(DetConfig);
	Engine *TheEngine = new Engine(FluxConfig, 2, 2);	//creating 2FHC and 2RHC fluxedrivers

	//Detector * TheBox = new Detector(DetConfig);

	//neutrino left has same flux as antineutrino right
	//neutrino right has same flux as antineutrino left
	//Neutrino *N_L = new Neutrino(0, Neutrino::Dirac | Neutrino::Left );
	//Neutrino *N_R = new Neutrino(0, Neutrino::Dirac | Neutrino::Right);
	Neutrino *TheNu0_L, *TheNu0_R, *TheNuB_L, *TheNuB_R;
	if (Dirac)
	{
		TheNu0_L = new Neutrino(0, Neutrino::Left  | Neutrino::Dirac);
		TheNuB_L = new Neutrino(0, Neutrino::Left  | Neutrino::Dirac | Neutrino::Antiparticle);
		TheNu0_R = new Neutrino(0, Neutrino::Right | Neutrino::Dirac);
		TheNuB_R = new Neutrino(0, Neutrino::Right | Neutrino::Dirac | Neutrino::Antiparticle);
		BaseName += "_dirac_";
	}
	else
	{
		TheNu0_L = TheNuB_L = new Neutrino(0, Neutrino::Left  | Neutrino::Majorana);
		TheNu0_R = TheNuB_R = new Neutrino(0, Neutrino::Right | Neutrino::Majorana);
		BaseName += "_major_";
	}

	TheNu0_L->SetMixings(5.0e-3, 5.0e-3, 5.0e-3);
	TheNuB_L->SetMixings(5.0e-3, 5.0e-3, 5.0e-3);
	TheNu0_R->SetMixings(5.0e-3, 5.0e-3, 5.0e-3);
	TheNuB_R->SetMixings(5.0e-3, 5.0e-3, 5.0e-3);
	TheNu0_L->SetDecayChannel(Channel);
	TheNuB_L->SetDecayChannel(Channel);
	TheNu0_R->SetDecayChannel(Channel);
	TheNuB_R->SetDecayChannel(Channel);

	TheEngine->BindNeutrino(TheNu0_L, Engine::FHC, 0);
	TheEngine->BindNeutrino(TheNu0_R, Engine::FHC, 1);
	TheEngine->BindNeutrino(TheNuB_L, Engine::RHC, 0);
	TheEngine->BindNeutrino(TheNuB_R, Engine::RHC, 1);

	std::stringstream ssL;
	double Mass;
	for (double Mass = 0.45; Mass < 0.46; Mass += 0.05)
	{
		TheNu0_L->SetMass(Mass);
		TheNuB_L->SetMass(Mass);
		TheNu0_R->SetMass(Mass);
		TheNuB_R->SetMass(Mass);

		TheEngine->MakeFlux();
		TheEngine->ScaleArea(TheBox);
		TheEngine->ScaleBaseline(TheBox);
		TheEngine->ScalePOT(1e20);

		std::vector<double> vWeight;
		double Total = TheEngine->MakeSampler(TheBox, vWeight);
		for (unsigned int i = 0; i < vWeight.size(); ++i)
			std::cout << vWeight.at(i) << std::endl;

		ssL.str("");
		ssL.clear();
		ssL << BaseName << std::setfill('0') << std::setw(4) << Mass*1000 << ".dat";
		std::cout << "Mass " << Mass << "\t in " << ssL.str() << std::endl;

		std::ofstream Out(ssL.str().c_str());

		double Start, End;
		double EnStep = TheEngine->RangeWidth(Start, End);
		for (double Energy = Start + EnStep/2.0; Energy < End; Energy += EnStep)
		{
			TheNu0_L->SetEnergyKin(Energy);
			TheNuB_L->SetEnergyKin(Energy);
			TheNu0_R->SetEnergyKin(Energy);
			TheNuB_R->SetEnergyKin(Energy);

			Out << Energy << "\t"; 
			if (Dirac)
			{
				Out << TheEngine->IntensitySample(Engine::FHC, 0) << "\t"; 
				Out << TheEngine->IntensitySample(Engine::FHC, 1) << "\t"; 
				Out << TheEngine->IntensitySample(Engine::RHC, 0) << "\t"; 
				Out << TheEngine->IntensitySample(Engine::RHC, 1) << "\t"; 
				Out << std::endl;
			}
			else
			{
				Out << TheEngine->IntensitySample(Engine::FHC, 0) +
				       TheEngine->IntensitySample(Engine::RHC, 0) << "\t"; 
				Out << TheEngine->IntensitySample(Engine::FHC, 1) +
				       TheEngine->IntensitySample(Engine::RHC, 1) << "\t"; 
				Out << std::endl;
			}
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
