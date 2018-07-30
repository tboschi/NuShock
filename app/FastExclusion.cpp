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
		{"threshold", 	required_argument,	0, 't'},
		{"efficiency", 	no_argument,		0, 'W'},
		{"output", 	required_argument,	0, 'o'},
		{"particle", 	no_argument,		0, 'P'},
		{"antipart", 	no_argument,		0, 'A'},
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
	bool Particle = false, Antipart = false;	//default majorana

	while((iarg = getopt_long(argc,argv, "d:f:c:o:t:q:WEMTLRAPBmh", longopts, &index)) != -1)
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

	Detector *TheBox = new Detector(DetConfig);
	std::vector<char> vFlag;

	Neutrino *TheNu;
	unsigned int OptHel, OptFerm;

	std::string First;
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

	if (Particle)
	{
		OptFerm = Neutrino::Dirac;
		FileName += "_p";
		if (Efficiency)
			TheBox->SetEfficiency(Channel, 1);
	}
	else if (Antipart)
	{
		OptFerm = Neutrino::Dirac | Neutrino::Antiparticle;
		FileName += "_a";
		if (Efficiency)
			TheBox->SetEfficiency(Channel, 1);
	}
	else
	{
		OptFerm = Neutrino::Majorana;
		FileName += "_m";
		if (Efficiency)
			TheBox->SetEfficiency(Channel, 0);
	}

	if (Left)
	{
		OptHel = Neutrino::Left;
		FileName += "_L";
	}
	else if (Right)
	{
		OptHel = Neutrino::Right;
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

	TheNu = new Neutrino(0, OptHel | OptFerm);
	TheNu->SetDecayChannel(Channel);

	Engine *TheEngine = new Engine(FluxConfig, 1, 1);	//creating 1FHC and 1RHC fluxedrivers
	TheEngine->BindNeutrino(TheNu, Engine::FHC, 0);		//left neutrino
	TheEngine->BindNeutrino(TheNu, Engine::RHC, 0);		//is a right antineutrino

	unsigned int Grid = 250;
	unsigned int nD = vFlag.size();	//number of dimensions
	double Mass;
	//std::cout << "Scanning over " << nD << " dimensions" << std::endl;

	Exclusion *Solver;
	if (TheNu->IsParticle())
		Solver = new Exclusion(TheEngine, Engine::FHC, TheBox, vFlag, Thr);
	else if (TheNu->IsAntiparticle())
		Solver = new Exclusion(TheEngine, Engine::RHC, TheBox, vFlag, Thr);
	else if (TheNu->IsMajorana())
		Solver = new Exclusion(TheEngine, Engine::Both, TheBox, vFlag, Thr);

	std::vector<double> vMass, vU2Bot, vU2Top;

	for (double logMass = -2.0; logMass < 0.3; logMass += 2.3/Grid)	//increase mass log
	{
		Mass = pow(10.0, logMass);
		//std::cout << "Mass " << Mass << std::endl;
		double Threshold = Thr + Mass * Qct;
		if (Threshold < 2.44)
			Threshold = 2.44;
		Solver->SetThr(Threshold);

		TheNu->SetMass(Mass);

		if (TheNu->IsDecayAllowed() &&
		    TheNu->IsProductionAllowed())
		{

			TheEngine->MakeFlux();
			TheEngine->ScaleDetector(TheBox);

			double lU2Bot = -12.0;
			double lU2Top = - 0.0;
			double lU2Mid;

			if (Solver->FindInterval(lU2Bot, lU2Mid, lU2Top))
			{
				lU2Bot = Solver->Bisect(lU2Bot, lU2Mid);
				lU2Top = Solver->Bisect(lU2Mid, lU2Top);

				vMass.push_back(Mass);
				vU2Bot.push_back(pow(10, lU2Bot));
				vU2Top.push_back(pow(10, lU2Top));

				//Out << Mass << "\t" << pow(10, lU2Bot) << "\t" << pow(10, lU2Top) << std::endl;
			}
		}
	}


	for (unsigned int i = 0; i < vMass.size(); ++i)
		Out << vMass.at(i) << "\t" << vU2Bot.at(i) << std::endl;;
	for (unsigned int i = vMass.size(); i > 0; --i)
		Out << vMass.at(i-1) << "\t" << vU2Top.at(i-1) << std::endl;;
	if (vMass.size())
		Out << vMass.front() << "\t" << vU2Bot.front() << std::endl;;

	delete TheNu;
	delete TheEngine;
	delete Solver;

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
	std::cout <<"\n  -c,  --channel" << std::endl;
	std::cout << "\t\tDecay channel, defaul ALL" << std::endl;
	std::cout <<"\n  -E,  -M,  -T" << std::endl;
	std::cout << "\t\tSelect which mixing element" << std::endl;
	std::cout <<"\n  -h,  --help" << std::endl;
	std::cout << "\t\tPrint this message and exit" << std::endl;
}
