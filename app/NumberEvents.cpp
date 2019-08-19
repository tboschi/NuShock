#include <iostream>
#include <fstream>
#include <vector>
#include <cstring>
#include <getopt.h>

#include "tools.h"
#include "flux.h"
#include "detector.h"
#include "physics.h"
#include "analysis.h"

void Usage(char* argv0);
int main(int argc, char** argv)
{
	const struct option longopts[] = 
	{
		{"detconfig", 	required_argument,	0, 'd'},
		{"fluxconfig", 	required_argument,	0, 'f'},
		{"channel", 	required_argument,	0, 'c'},
		{"output", 	required_argument,	0, 'o'},
		{"threshold", 	required_argument,	0, 't'},
		{"massdepend", 	required_argument,	0, 'q'},
		{"efficiency", 	no_argument,		0, 'W'},
		{"majorana", 	no_argument,		0, 'j'},
		{"help", 	no_argument,	 	0, 'h'},
		{0,	0, 	0,	0},
	};

	int index; 
	int iarg = 0;
	opterr = 1;
	
	//Initialize variables
	std::string detConfig, fluxConfig, fileName;
	std::string channel = "ALL";
	double ue = 0, um = 0, ut = 0;
	double mass = 0;
	
	bool left = false, right = false;	//default unpolarised

	while((iarg = getopt_long(argc,argv, "d:f:c:m:E:M:T:LRrjAPh", longopts, &index)) != -1)
	{
		switch(iarg)
		{
			case 'd':
				detConfig.assign(optarg);
				break;
			case 'f':
				fluxConfig.assign(optarg);
				break;
			case 'c':
				channel.assign(optarg);
				break;
			case 'm':
				mass = std::strtod(optarg, NULL);
				break;
			case 'E':
				ue = std::strtod(optarg, NULL);
				break;
			case 'M':
				um = std::strtod(optarg, NULL);
				break;
			case 'T':
				ut = std::strtod(optarg, NULL);
				break;
			case 'L':
				left = true;
				right = false;
				break;
			case 'R':
				left = false;
				right = true;
				break;
			case 'h':
				Usage(argv[0]);
				return 1;
			default:
				break;
		}
	}

	ue = std::sqrt(ue);
	um = std::sqrt(um);
	ut = std::sqrt(ut);

	//constructing the detector
	unsigned int opt0, optB, optHel;

	double mod = 1.0;

	if (left)
		optHel = Neutrino::Left;	//-1 helicity
	else if (right)
		optHel = Neutrino::Right;	//+1 helicity
	else
		optHel = Neutrino::Unpolarised;

	Detector *theBox = new Detector(detConfig);
	Engine *theEngine = new Engine(fluxConfig);

	Neutrino dNu0(mass, optHel | Neutrino::Dirac);
	Neutrino dNuB(mass, optHel | Neutrino::Dirac | Neutrino::Antiparticle);
	Neutrino mNu0(mass, optHel | Neutrino::Majorana);
	//Neutrino majorNuB(0, optHel | Neutrino::Majorana | Neutrino::Antiparticle);

	dNu0.SetDecayChannel(channel);
	dNuB.SetDecayChannel(channel);
	mNu0.SetDecayChannel(channel);
	//majorNuB.SetDecayChannel(channel);

	theEngine->BindNeutrino("dirac_nu0", dNu0, Engine::FHC);
	theEngine->BindNeutrino("dirac_nuB", dNuB, Engine::RHC);
	theEngine->BindNeutrino("major_nu0", mNu0, Engine::FHC);
	theEngine->BindNeutrino("major_nuB", mNu0, Engine::RHC);

	theEngine->MakeFlux();
	theEngine->ScaleToDetector(theBox);

	std::cout << "computing for " << mass << ", ("
		  << ue << ",\t" << um << ",\t" << ut << ")" << std::endl;
	std::cout << "and channel " << channel << std::endl;

	std::map<std::string, double> numevts;
	double tot = theEngine->MakeSampler(theBox, numevts, ue, um, ut);

	double vDirac0 = numevts["dirac_nu0"];
	double vDiracB = numevts["dirac_nuB"];
	double vMajor  = (numevts["major_nu0"]+numevts["major_nuB"]) / 2.0;

	std::cout << "Number of events " << vDirac0 << "\t" << vDiracB << "\t" << vMajor << std::endl;


	return 0;
}
	
void Usage(char* argv0)
{
	std::cout << "Description" << std::endl;
	std::cout << "Usage : " << std::endl;
	std::cout << argv0 << " [OPTIONS]" << std::endl;
	std::cout <<"\n  -d,  --detconfig [CONFIG]" << std::endl;
	std::cout << "\t\tDetector configuration file" << std::endl;
	std::cout <<"\n  -f,  --fluxconfig [CONFIG]" << std::endl;
	std::cout << "\t\tFlux configuration file" << std::endl;
	std::cout <<"\n  -c,  --channel [STRING]" << std::endl;
	std::cout << "\t\tDecay channel, defaul ALL" << std::endl;
	std::cout <<"\n  -E,  -M,  -T" << std::endl;
	std::cout << "\t\tSelect which mixing element. Multiple selection allowed" << std::endl;
	std::cout <<"\n  -L,  -R,  -U" << std::endl;
	std::cout << "\t\tSelect neutrino polarisation." << std::endl;
	std::cout <<"\n  -h,  --help" << std::endl;
	std::cout << "\t\tPrint this message and exit" << std::endl;
}

