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

#include "omp.h"

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
		{"efficiency", 	optional_argument,	0, 'w'},
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
	std::string detConfig, module, fluxConfig, fileName, background;
	std::string channel = "ALL";
	bool UeFlag = false, UmFlag = false, UtFlag = false;
	bool effSet = false;
	
	bool left = false, right = false;	//default unpolarised
	double CL = 0.99;			//confidence level

	while((iarg = getopt_long(argc,argv, "d:l:f:c:w::o:C:EMTLRh", longopts, &index)) != -1)
	{
		switch(iarg)
		{
			case 'd':
				detConfig.assign(optarg);
				break;
			case 'l':
				module.assign(optarg);
				break;
			case 'f':
				fluxConfig.assign(optarg);
				break;
			case 'c':
				channel.assign(optarg);
				break;
			case 'w':
				effSet = true;
				if (optarg)
					background.assign(optarg);
				break;
			case 'o':
				fileName.assign(optarg);
				break;
			case 'C':
				CL = std::strtod(optarg, NULL);
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

	//constructing the detector
	unsigned int optHel;

	if (!UeFlag && !UmFlag && !UtFlag)
	{
		std::cerr << "You have to select at least one mixing, -E -M or -T" << std::endl;
		return 1;
	}

	if (fileName.find(".dat") == std::string::npos)
		fileName.append(".dat");

	fileName.insert(fileName.find(".dat"), "_" + channel);

	if (UeFlag)
		fileName.insert(fileName.find(".dat"), "_E");
	if (UmFlag)
		fileName.insert(fileName.find(".dat"), "_M");
	if (UtFlag)
		fileName.insert(fileName.find(".dat"), "_T");
	

	/*
	if (Dirac)
	{
		opt0 = Neutrino::Dirac;
		optB = Neutrino::Dirac | Neutrino::Antiparticle;
		if (Particle)			//for LNV
			FileName += "_dP";
		else if (Antipart)		//for LNV
			FileName += "_dA";
		else
			FileName += "_d";
		if (Efficiency)
			TheBox->SetEfficiency(Channel, 1);
	}
	else	//Majorana
	{
		opt0 = OptB = Neutrino::Majorana;
		if (Particle)			//for LNV
			FileName += "_mP";
		else if (Antipart)		//for LNV
			FileName += "_mA";
		else
			FileName += "_m";
		//FileName += "_m";
		if (Efficiency)
			TheBox->SetEfficiency(Channel, 0);
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

	Detector *box_LNC_dir = new Detector(detConfig, module);
	Detector *box_LNC_maj = new Detector(detConfig, module);
	Detector *box_LNV_dir = new Detector(detConfig, module);
	Detector *box_LNV_maj = new Detector(detConfig, module);

	std::map<std::string, std::vector<double> > mCurve; //a collection of parabolas for each channel
	std::map<std::string, std::vector<double> >::iterator ic;

	std::string pre = "c";
	std::string key = channel + "_";
	if (!background.empty())
	{
		CardDealer threshold(background);
		std::map<std::string, std::vector<double> > mBeta;
		if (threshold.Get(key, mBeta))
		{
			std::string ferm[2] = {"dirac", "major"};
			std::string lnv[2] = {"LNV", "LNC"};

			for (int f = 0; f < 2; ++f)
				for (int l = 0; l < 2; ++l)
			{
				std::vector<double> alpha(3);
				ic = mBeta.begin();
				while (ic != mBeta.end())
				{
					if (ic->first.find(ferm[f]) != std::string::npos &&
					    ic->first.find(lnv[l])  != std::string::npos &&
					    ic->first.find(module)  != std::string::npos)
					{
						std::cout << "adding " << ic->first << std::endl;
						for (int j = 0; j < alpha.size(); ++j)
							alpha[j] += ic->second[j];
					}
					++ic;
				}
				std::cout << "Filling " << ferm[f] + "_" + lnv[l] << std::endl;
				mCurve[ferm[f] + "_" + lnv[l]] = alpha;
			}

			pre = "b";
		}

		fileName.insert(fileName.find(".dat"), "_backgr");
	}

	if (effSet && background.empty())
		fileName.insert(fileName.find(".dat"), "_charge");

	if (effSet && module.empty())
	{
		box_LNC_dir->SetEfficiency(pre + key + "LNC_FGT_dirac", channel);
		box_LNV_dir->SetEfficiency(pre + key + "LNV_FGT_dirac", channel);
		box_LNC_dir->SetEfficiency(pre + key + "LNC_LAr_dirac", channel);
		box_LNV_dir->SetEfficiency(pre + key + "LNV_LAr_dirac", channel);

		box_LNC_maj->SetEfficiency(pre + key + "LNC_FGT_major", channel);
		box_LNV_maj->SetEfficiency(pre + key + "LNV_FGT_major", channel);
		box_LNC_maj->SetEfficiency(pre + key + "LNC_LAr_major", channel);
		box_LNV_maj->SetEfficiency(pre + key + "LNV_LAr_major", channel);
	}
	else if (effSet)
	{
		box_LNC_dir->SetEfficiency(pre + key + "LNC_" + module + "_dirac", channel);
		box_LNV_dir->SetEfficiency(pre + key + "LNV_" + module + "_dirac", channel);

		box_LNC_maj->SetEfficiency(pre + key + "LNC_" + module + "_major", channel);
		box_LNV_maj->SetEfficiency(pre + key + "LNV_" + module + "_major", channel);
	}

	int grid = 500;

	double m0 = 0.01, mE = 2.0;
	double U20 = 1e-10, U2E = 1;

	Neutrino nuthr(0, optHel | Neutrino::Dirac);
	double mass0 = nuthr.DecayThreshold(channel);
	double massE = nuthr.ProductionThreshold();
	if (mass0 > massE)
	{
		double tmp = massE;
		massE = mass0;
		mass0 = tmp;
	}

	std::vector<std::vector<bool> > vGrid;
	std::vector<double> vMass;
	//std::vector<double> vU2top, vU2low;

	Engine *theEngine = new Engine(fluxConfig);

	std::cout << "Saving in " << fileName << std::endl;

	for (int i = 0; i < grid; ++i)
	{
		double mass = m0 * pow(mE/m0, i / double(grid));

		if (mass < mass0 || mass > massE)
			continue;

		std::cout << "mass " << mass << std::endl;

		vMass.push_back(mass);

		Neutrino nu0_dir(mass, optHel | Neutrino::Dirac);
		Neutrino nuB_dir(mass, optHel | Neutrino::Dirac | Neutrino::Antiparticle);
		Neutrino nu__maj(mass, optHel | Neutrino::Majorana);

		nu0_dir.SetDecayChannel(channel);
		nuB_dir.SetDecayChannel(channel);
		nu__maj.SetDecayChannel(channel);

		//creating 1FHC and 1RHC fluxedrivers
		//Engine *fluxDirac = new Engine(fluxConfig);
		//Engine *fluxMajor = new Engine(fluxConfig);

		theEngine->Reset();

		theEngine->BindNeutrino("dirac_nu0", nu0_dir, Engine::FHC);
		theEngine->BindNeutrino("dirac_nuB", nuB_dir, Engine::RHC);
		theEngine->BindNeutrino("major_nu0", nu__maj, Engine::FHC);
		theEngine->BindNeutrino("major_nuB", nu__maj, Engine::RHC);

		theEngine->MakeFlux();
		//each detector has same dimensions, only efficiency is different
		theEngine->ScaleToDetector(box_LNC_dir);

		int bak_LNC_dir = 0, bak_LNV_dir = 0;
		int bak_LNC_maj = 0, bak_LNV_maj = 0;
		if (!background.empty())
		{
			bak_LNC_dir = ceil(mCurve["dirac_LNC"][0] + mass *
					  (mCurve["dirac_LNC"][1] + mass * mCurve["dirac_LNC"][2]));
			bak_LNV_dir = ceil(mCurve["dirac_LNV"][0] + mass *
					  (mCurve["dirac_LNV"][1] + mass * mCurve["dirac_LNV"][2]));

			bak_LNC_maj = ceil(mCurve["major_LNC"][0] + mass *
					  (mCurve["major_LNC"][1] + mass * mCurve["major_LNC"][2]));
			bak_LNV_maj = ceil(mCurve["major_LNV"][0] + mass *
					  (mCurve["major_LNV"][1] + mass * mCurve["major_LNV"][2]));
			std::cout << "Adding backgrounds"
				  << "\n\tDirac LNC : " << bak_LNC_dir << ",\tLNV : " << bak_LNV_dir
				  << "\n\tMajor LNC : " << bak_LNC_maj << ",\tLNV : " << bak_LNV_maj
				  << std::endl;
		}


		//std::vector<bool> vU2top, vU2bot;
		std::vector<bool> vU2;
		for (int j = 0; j < grid; ++j)
		{
			double U2 = U20 * pow(U2E/U20, j / double(grid));

			double ue = sqrt(U2) * (UeFlag ? 1.0 : 1.0/100.0);
			double um = sqrt(U2) * (UmFlag ? 1.0 : 1.0/100.0);
			double ut = sqrt(U2) * (UtFlag ? 1.0 : 1.0/100.0);

			//diracnu0.SetMixings(ue, um, ut);
			//diracnuB.SetMixings(ue, um, ut);
			//majorana.SetMixings(ue, um, ut);

			std::map<std::string, double> diracNum, majorNum;
			double sig_LNC_dir = theEngine->MakeSampler(box_LNC_dir, "dirac_nu0", ue, um, ut);
			double sig_LNV_dir = theEngine->MakeSampler(box_LNV_dir, "dirac_nuB", ue, um, ut);
			double sig_LNC_maj = theEngine->MakeSampler(box_LNC_maj, "major_nu0", ue, um, ut);
			double sig_LNV_maj = theEngine->MakeSampler(box_LNV_maj, "major_nuB", ue, um, ut);

			sig_LNC_maj = sig_LNV_maj = (sig_LNC_maj + sig_LNV_maj) / 2.0;

			//std::cout << "With mixing " << U2
			//	  << "\n\tDirac LNC : " << LNC_dirac << ",\tLNV : " << LNV_dirac
			//	  << "\n\tMajor LNC : " << LNC_major << ",\tLNV : " << LNV_major
			//	  << "\n\tALL : " << all
			//	  << "\n\tration = " << (LNC_major - LNC_dirac) << std::endl;
			//std::cout << U2 << "\t"
			//	  << LNC_dirac << "\t"
			//	  << LNV_dirac << "\t"
			//	  << LNC_major << "\t";
			//std::cout << std::endl;
			//std::cout << "\tEvents: " << vDirac0[k] << "\t" << vDiracB[k] << "\t"
			//	  << vMajor[k] << "(total " << tot << ")" << std::endl;

			vU2.push_back( IsSensitiveToLNV(CL, sig_LNC_dir, sig_LNV_dir,
							    sig_LNC_maj, sig_LNV_maj,
							    bak_LNC_maj, bak_LNV_maj,
							    bak_LNC_maj, bak_LNV_maj ) );
			//std::cout << "\t" << U2 << " : " << std::boolalpha << vU2.back() << std::endl;
		}

		vGrid.push_back(vU2);
	}

	std::ofstream out(fileName.c_str());

	std::vector<double> vU2low(vGrid.size()), vU2top(vGrid.size());
	for (int i = 0; i < vGrid.size(); ++i)
	{
		bool prev = false;
		vU2top[i] = 1;
		vU2low[i] = 1;
		for (int j = 0; j < vGrid[i].size(); ++j)
		{
			double U2 = U20 * pow(U2E/U20, j / double(grid));
			if (vGrid[i][j] != prev)
			{
				if (prev)
					vU2top[i] = U2;
				else 
					vU2low[i] = U2;

				prev = !prev;
			}
		}
	}

	for (int i = 0; i < vMass.size(); ++i)
		out << vMass[i] << "\t" << vU2low[i] << std::endl;
	for (int i = vMass.size(); i > 0; --i)
		out << vMass[i-1] << "\t" << vU2top[i-1] << std::endl;
	out << vMass.front() << "\t" << vU2low.front() << std::endl;

	out.close();

	delete box_LNC_dir;
	delete box_LNV_dir;
	delete box_LNC_maj;
	delete box_LNV_maj;
	delete theEngine;

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
