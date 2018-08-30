#include <iostream>
#include <fstream>
#include <vector>
#include <cstring>
#include <getopt.h>

#include "Tools.h"
#include "Flux.h"
#include "Detector.h"
#include "Physics.h"

#include "Analysis.h"

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
		{"particle", 	no_argument,		0, 'P'},
		{"antipart", 	no_argument,		0, 'A'},
		{"dirac", 	no_argument,		0, 'r'},
		{"majorana", 	no_argument,		0, 'j'},
		{"help", 	no_argument,	 	0, 'h'},
		{0,	0, 	0,	0},
	};

	int index; 
	int iarg = 0;
	opterr = 1;
	
	//Initialize variables
	std::string DetConfig, FluxConfig, FileName;
	std::string Channel = "ALL";
	bool UeFlag = false, UmFlag = false, UtFlag = false;
	bool Efficiency = false;
	double Thr = 2.44, Qct = 0.0;
	//for mass dependency of threshold as in T = Qct * Mass + Thr
	
	bool Left = false, Right = false;		//default unpolarised
	bool Particle = false, Antipart = false;	//default both for Dirac
	bool Dirac = false;				//default majorana neutrino

	while((iarg = getopt_long(argc,argv, "d:f:c:o:t:q:WEMTLRAPjrh", longopts, &index)) != -1)
	{
		switch(iarg)
		{
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
				FileName.assign(optarg);
				break;
			case 't':
				Thr = std::strtod(optarg, NULL);
				break;
			case 'q':
				Qct = std::strtod(optarg, NULL);
				break;
			case 'W':
				Efficiency = true;
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
			case 'L':
				Left = true;
				Right = false;
				break;
			case 'R':
				Left = false;
				Right = true;
				break;
			case 'r':
				Dirac = true;
				break;
			case 'j':
				Dirac = false;
				break;
			case 'P':
				Particle = true;
				Antipart = false;
				break;
			case 'A':
				Particle = false;
				Antipart = true;
				break;
			case 'h':
				Usage(argv[0]);
				return 1;
			default:
				break;
		}
	}

	//constructing the detector
	Detector *TheBox = new Detector(DetConfig);
	std::vector<char> vFlag;

	unsigned int OptHel, Opt0, OptB;

	std::string First;

	if (!UeFlag && !UmFlag && !UtFlag)
	{
		std::cerr << "You have to select at least one mixing, -E -M or -T" << std::endl;
		return 1;
	}

	if (UeFlag)
	{
		vFlag.push_back('E');
		FileName += "_E";
		First += "Ue\t";
	}
	if (UmFlag)
	{
		vFlag.push_back('M');
		FileName += "_M";
		First += "Um\t";
	}
	if (UtFlag)
	{
		vFlag.push_back('T');
		FileName += "_T";
		First += "Ut\t";
	}
	

	if (Dirac)
	{
		Opt0 = Neutrino::Dirac;
		OptB = Neutrino::Dirac | Neutrino::Antiparticle;
		if (Particle)
			FileName += "_dP";
		else if (Antipart)
			FileName += "_dA";
		else
			FileName += "_d";
		if (Efficiency)
			TheBox->SetEfficiency(Channel, 1);
	}
	else	//Majorana
	{
		Opt0 = OptB = Neutrino::Majorana;
		FileName += "_m";
		if (Efficiency)
			TheBox->SetEfficiency(Channel, 0);
	}

	if (Left)
	{
		OptHel = Neutrino::Left;	//-1 helicity
		FileName += "_L";
	}
	else if (Right)
	{
		OptHel = Neutrino::Right;	//+1 helicity
		FileName += "_R";
	}
	else
	{
		OptHel = Neutrino::Unpolarised;
		FileName += "_U";
	}

	FileName += ".dat";
	std::ofstream Out(FileName.c_str());
	Out << "#Mass\t" << First << "Events" << std::endl;

	Neutrino *TheNu0 = new Neutrino(0, OptHel | Opt0);
	Neutrino *TheNuB = new Neutrino(0, OptHel | OptB);

	TheNu0->SetDecayChannel(Channel);
	TheNuB->SetDecayChannel(Channel);

	Engine *TheEngine = new Engine(FluxConfig, 1, 1);	//creating 1FHC and 1RHC fluxedrivers

	TheEngine->BindNeutrino(TheNu0, Engine::FHC, 0);		//left neutrino
	TheEngine->BindNeutrino(TheNuB, Engine::RHC, 0);		//is a right antineutrino

	unsigned int Grid = 500;
	unsigned int nD = vFlag.size();	//number of dimensions
	double Mass, lMStart = log10(0.01), lMEnd = log10(2.0), lMStep = log10(2.0/0.01)/Grid;
	//std::cout << "Scanning over " << nD << " dimensions" << std::endl;

	Exclusion *Solver;
	if (Dirac && Particle)
		Solver = new Exclusion(TheEngine, Engine::FHC, TheBox, vFlag, Thr);
	else if (Dirac && Antipart)
		Solver = new Exclusion(TheEngine, Engine::RHC, TheBox, vFlag, Thr);
	else	//Dirac both or Majorana
		Solver = new Exclusion(TheEngine, Engine::Both, TheBox, vFlag, Thr);

	std::vector<double> vMass, vU2Bot, vU2Top;

	//logscale loop
	for (double logMass = lMStart; logMass < lMEnd; logMass += lMStep)
	{
		Mass = pow(10.0, logMass);
		double Threshold = Thr + Mass * Qct;
		if (Threshold < 2.44)
			Threshold = 2.44;
		Solver->SetThr(Threshold);

		TheNu0->SetMass(Mass);
		TheNuB->SetMass(Mass);

		if (TheNu0->IsDecayAllowed() &&
		    TheNu0->IsProductionAllowed() &&
		    TheNuB->IsDecayAllowed() &&
		    TheNuB->IsProductionAllowed())
		{
			TheEngine->MakeFlux();
			TheEngine->ScaleDetector(TheBox);

			double lU2Bot = -16.0;
			double lU2Top =   0.0;
			double lU2Mid;

			if (Solver->FindInterval(lU2Bot, lU2Mid, lU2Top))
			{
				lU2Bot = Solver->Bisect(lU2Bot, lU2Mid);
				lU2Top = Solver->Bisect(lU2Mid, lU2Top);

				vMass.push_back(Mass);
				vU2Bot.push_back(pow(10, lU2Bot));
				vU2Top.push_back(pow(10, lU2Top));
			}
			else
			{
				for (unsigned int i = 0; i < vMass.size(); ++i)
					Out << vMass.at(i) << "\t" << vU2Bot.at(i) << std::endl;;
				for (unsigned int i = vMass.size(); i > 0; --i)
					Out << vMass.at(i-1) << "\t" << vU2Top.at(i-1) << std::endl;;
				if (vMass.size())
				{
					Out << vMass.front() << "\t" << vU2Bot.front() << std::endl;;
					Out << std::endl << std::endl;
				}

				vMass.clear();
				vU2Bot.clear();
				vU2Top.clear();
			}
		}
	}


	for (unsigned int i = 0; i < vMass.size(); ++i)
		Out << vMass.at(i) << "\t" << vU2Bot.at(i) << std::endl;;
	for (unsigned int i = vMass.size(); i > 0; --i)
		Out << vMass.at(i-1) << "\t" << vU2Top.at(i-1) << std::endl;;
	if (vMass.size())
		Out << vMass.front() << "\t" << vU2Bot.front() << std::endl;;

	delete TheNu0;
	delete TheNuB;
	delete TheEngine;
	delete Solver;

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
	std::cout <<"\n  -o,  --output [path/STRING]" << std::endl;
	std::cout << "\t\tBase name of output file" << std::endl;
	std::cout <<"\n  -c,  --channel [STRING]" << std::endl;
	std::cout << "\t\tDecay channel, defaul ALL" << std::endl;
	std::cout <<"\n  -t,  --threshold [FLOAT]" << std::endl;
	std::cout << "\t\tThreshold intercept. If -q not set, then is mass independent. Defaul 2.44" << std::endl;
	std::cout <<"\n  -q,  --massdepend [FLOAT]" << std::endl;
	std::cout << "\t\tMass dependance of threshold. Needs -t flag. Defaul 0" << std::endl;
	std::cout <<"\n  -E,  -M,  -T" << std::endl;
	std::cout << "\t\tSelect which mixing element. Multiple selection allowed" << std::endl;
	std::cout <<"\n  -L,  -R,  -U" << std::endl;
	std::cout << "\t\tSelect neutrino polarisation." << std::endl;
	std::cout <<"\n  --dirac,  --majorana" << std::endl;
	std::cout << "\t\tSelect fermionic nature of neutrino. Default Majorana" << std::endl;
	std::cout <<"\n  -P, -A {incompatible with --majorana}" << std::endl;
	std::cout << "\t\tSelect either neutrino or antineutrino components for dirac. Default, both are used." << std::endl;
	std::cout <<"\n  -h,  --help" << std::endl;
	std::cout << "\t\tPrint this message and exit" << std::endl;
}
