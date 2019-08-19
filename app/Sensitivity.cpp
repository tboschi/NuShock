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
	std::string detConfig, fluxConfig, fileName;
	std::string channel = "ALL";
	bool UeFlag = false, UmFlag = false, UtFlag = false;
	
	bool left = false, right = false;	//default unpolarised
	bool particle = false, antipart = false;	//default both for dirac, used in lnv studies
	bool dirac = false;				//default majorana neutrino
	double CL = 0.90;			//confidence level
	double thr = 2.44, qct = 0.0;
	double ue = -1, um = -1, ut = -1;

	while((iarg = getopt_long(argc,argv, "d:f:c:o:C:t:q:EMTLRrjAPh", longopts, &index)) != -1)
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
			case 'o':
				fileName.assign(optarg);
				break;
			case 'C':
				CL = std::strtod(optarg, NULL);
				break;
			case 't':
				thr = std::strtod(optarg, NULL);
				break;
			case 'q':
				qct = std::strtod(optarg, NULL);
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
			case 'r':
				dirac = true;
				break;
			case 'j':
				dirac = false;
				break;
			case 'P':
				particle = true;
				antipart = false;
				break;
			case 'A':
				particle = false;
				antipart = true;
				break;
			case 'h':
				Usage(argv[0]);
				return 1;
			default:
				break;
		}
	}

	//constructing the detector
	unsigned int opt0, optB, optHel;

	if (ue < 0 && um < 0 && ut < 0)
	{
		std::cerr << "You have to select at least one mixing, -E -M or -T" << std::endl;
		return 1;
	}

	std::vector<char> vFlag;

	if (fileName.find(".dat") == std::string::npos)
		fileName += ".dat";

	if (UeFlag)
		fileName.insert(fileName.find(".dat"), "_E");
	if (UmFlag)
		fileName.insert(fileName.find(".dat"), "_M");
	if (UtFlag)
		fileName.insert(fileName.find(".dat"), "_T");
	
	Engine::Current horn;
	double mod = 1.0;


	if (particle)			//for LNV
		horn = Engine::FHC;
	else if (antipart)		//for LNV
		horn = Engine::RHC;
	else
		horn = Engine::both;

	/*
	if (dirac)
	{
		//opt0 = Neutrino::Dirac;
		//optB = Neutrino::Dirac | Neutrino::Antiparticle;

		if (particle)			//for LNV
		{
			horn = Engine::FHC;
			fileName.insert(fileName.find(".dat"), "_dP");
		}
		else if (antipart)		//for LNV
		{
			horn = Engine::RHC;
			fileName.insert(fileName.find(".dat"), "_dA");
		}
		else
		{
			horn = Engine::both;
			fileName.insert(fileName.find(".dat"), "_d");
		}

		//if (Efficiency)
		//	TheBox->SetEfficiency(Channel, 1);
	}
	else	//Majorana
	{
		//opt0 = optB = Neutrino::Majorana;

		if (particle)			//for lnv
		{
			mod = 0.5;
			horn = Engine::FHC;
			fileName.insert(fileName.find(".dat"), "_mP");
		}
		else if (antipart)		//for LNV
		{
			mod = 0.5;
			horn = Engine::RHC;
			fileName.insert(fileName.find(".dat"), "_mA");
		}
		else
		{
			horn = Engine::both;
			fileName.insert(fileName.find(".dat"), "_m");
		}

		//if (Efficiency)
		//	TheBox->SetEfficiency(Channel, 0);
	}
	*/

	if (left)
	{
		optHel = Neutrino::Left;	//-1 helicity
		fileName.insert(fileName.find(".dat"), "_L");
	}
	else if (right)
	{
		optHel = Neutrino::Right;	//+1 helicity
		fileName.insert(fileName.find(".dat"), "_R");
	}
	else
	{
		optHel = Neutrino::Unpolarised;
		fileName.insert(fileName.find(".dat"), "_U");
	}

	std::string dirFile = fileName;
	dirFile.insert(dirFile.find(".dat"), "_d");
	std::string majFile = fileName;
	majFile.insert(majFile.find(".dat"), "_m");

	std::ofstream outdir(dirFile.c_str());
	std::ofstream outmaj(majFile.c_str());

	Detector *theBox = new Detector(detConfig);
	Engine *dirEngine = new Engine(fluxConfig);
	Engine *majEngine = new Engine(fluxConfig);

	Exclusion *dirSolver = new Exclusion(dirEngine, theBox, horn,
						UeFlag, UmFlag, UtFlag, thr, mod);
	Exclusion *majSolver = new Exclusion(majEngine, theBox, horn,
						UeFlag, UmFlag, UtFlag, thr, mod);

	Neutrino diracNu0(0, optHel | Neutrino::Dirac);
	Neutrino diracNuB(0, optHel | Neutrino::Dirac | Neutrino::Antiparticle);
	Neutrino majorNu0(0, optHel | Neutrino::Majorana);
	//Neutrino majorNuB(0, optHel | Neutrino::Majorana | Neutrino::Antiparticle);

	diracNu0.SetDecayChannel(channel);
	diracNuB.SetDecayChannel(channel);
	majorNu0.SetDecayChannel(channel);
	//majorNuB.SetDecayChannel(channel);

	dirEngine->BindNeutrino("dirac_nu0", diracNu0, Engine::FHC);
	dirEngine->BindNeutrino("dirac_nuB", diracNuB, Engine::RHC);
	majEngine->BindNeutrino("major_nu0", majorNu0, Engine::FHC);
	majEngine->BindNeutrino("major_nuB", majorNu0, Engine::RHC);

	int grid = 500;
	double m0 = 0.01, mE = 2.0;
	double U20 = 1e-12, U2E = 1;

	double mass0 = diracNu0.DecayThreshold(channel);
	double massE = diracNu0.ProductionThreshold();
	if (mass0 > massE)
	{
		double tmp = massE;
		massE = mass0;
		mass0 = tmp;
	}

	std::vector<double> dirMass, dirU2bot, dirU2top, dirDirbot, dirDirtop, dirMajbot, dirMajtop;
	std::vector<double> majMass, majU2bot, majU2top, majDirbot, majDirtop, majMajbot, majMajtop;

	for (int i = 0; i < grid; ++i)
	{
		double mass = m0 * pow(mE/m0, i / double(grid));
		std::cout << "mass " << mass << std::endl;

		if (mass < mass0 || mass > massE)
			continue;

		diracNu0.SetMass(mass);
		diracNuB.SetMass(mass);
		majorNu0.SetMass(mass);

		double threshold = thr + mass * qct;
		if (threshold < 2.44)
			threshold = 2.44;

		//solver->SetThr(threshold);

		//creating 1FHC and 1RHC fluxedrivers
		//Engine *fluxDirac = new Engine(fluxConfig);
		//Engine *fluxMajor = new Engine(fluxConfig);

		//theEngine->BindNeutrino("major_nu0", majorana, Engine::FHC);
		//theEngine->BindNeutrino("major_nuB", majorana, Engine::RHC);

		dirEngine->MakeFlux();
		majEngine->MakeFlux();
		dirEngine->ScaleToDetector(theBox);
		majEngine->ScaleToDetector(theBox);

		////////////////////////////////////////

		double lU2bot = -16.0;
		double lU2top =   0.0;
		double lU2mid;

		//dirac solver
		if (dirSolver->FindInterval(lU2bot, lU2mid, lU2top))
		{
			double dirbot, dirtop;
			double majbot, majtop;

			lU2bot = dirSolver->Bisect(lU2bot, lU2mid, dirbot);
			lU2top = dirSolver->Bisect(lU2mid, lU2top, dirtop);

			dirMass.push_back(mass);
			dirU2bot.push_back(pow(10, lU2bot));
			dirU2top.push_back(pow(10, lU2top));
			dirDirbot.push_back(dirbot);
			dirDirtop.push_back(dirtop);

			majbot = majSolver->NumberEvents(lU2bot);
			majtop = majSolver->NumberEvents(lU2top);
			dirMajbot.push_back(majbot);
			dirMajtop.push_back(majtop);
		}

		lU2bot = -16.0;
		lU2top =   0.0;

		//majorana solver
		if (majSolver->FindInterval(lU2bot, lU2mid, lU2top))
		{
			double dirbot, dirtop;
			double majbot, majtop;

			lU2bot = majSolver->Bisect(lU2bot, lU2mid, majbot);
			lU2top = majSolver->Bisect(lU2mid, lU2top, majtop);

			majMass.push_back(mass);
			majU2bot.push_back(pow(10, lU2bot));
			majU2top.push_back(pow(10, lU2top));
			majMajbot.push_back(majbot);
			majMajtop.push_back(majtop);

			dirbot = dirSolver->NumberEvents(lU2bot);
			dirtop = dirSolver->NumberEvents(lU2top);
			majDirbot.push_back(dirbot);
			majDirtop.push_back(dirtop);
		}
	}

	for (int i = 0; i < dirMass.size(); ++i)
		outdir << dirMass[i] << "\t" << dirU2bot[i] << "\t"
		       << dirDirbot[i] << "\t" << dirMajbot[i] << std::endl;

	for (int i = dirMass.size(); i > 0; --i)
		outdir << dirMass[i-1] << "\t" << dirU2top[i-1] << "\t"
		       << dirDirtop[i-1] << "\t" << dirMajtop[i-1] << std::endl;

	outdir << dirMass[0] << "\t" << dirU2bot[0] << "\t"
	       << dirDirbot[0] << "\t" << dirMajbot[0] << std::endl;

	for (int i = 0; i < majMass.size(); ++i)
		outmaj << majMass[i] << "\t" << majU2bot[i] << "\t"
		       << majDirbot[i] << "\t" << majMajbot[i] << std::endl;

	for (int i = majMass.size(); i > 0; --i)
		outmaj << majMass[i-1] << "\t" << majU2top[i-1] << "\t"
		       << majDirtop[i-1] << "\t" << majMajtop[i-1] << std::endl;

	outmaj << majMass[0] << "\t" << majU2bot[0] << "\t"
	       << majDirbot[0] << "\t" << majMajbot[0] << std::endl;



	//vMass.clear();
	//vU2bot.clear();
	//vU2top.clear();

	//delete theBox;
	//delete theEngine;
	//delete solver;

	//out.close();

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
