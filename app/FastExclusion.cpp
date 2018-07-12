#include <iostream>
#include <fstream>
#include <vector>
#include <cstring>
#include <getopt.h>

#include "Tools.h"
#include "Flux.h"
#include "Detector.h"
#include "Physics.h"

#include "TH2D.h"
#include "TFile.h"

void Usage(char* argv0);
int main(int argc, char** argv)
{
	const struct option longopts[] = 
	{
		{"smconfig", 	required_argument, 	0, 's'},
		{"detconfig", 	required_argument,	0, 'd'},
		{"fluxconfig", 	required_argument,	0, 'f'},
		{"channel", 	required_argument,	0, 'c'},
		{"threshold", 	required_argument,	0, 't'},
		{"efficiency", 	no_argument,		0, 'W'},
		{"output", 	required_argument,	0, 'o'},
		{"help", 	no_argument,	 	0, 'h'},
		{0,	0, 	0,	0},
	};

	int index; 
	int iarg = 0;
	opterr = 1;
	
	//Initialize variables
	std::string DetConfig, FluxConfig;
	std::ofstream OutFile;
	//TFile *OutFile;
	std::string Channel = "ALL";
	bool UeFlag = false, UmFlag = false, UtFlag = false;
	bool Left = false, Right = false;	//default unpolarised
	bool Efficiency = false;
	double Threshold = 2.44;
	
	while((iarg = getopt_long(argc,argv, "d:f:c:o:WEMTLRh", longopts, &index)) != -1)
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
				OutFile.open(optarg);
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
			case 't':
				Threshold = std::strtod(optarg, NULL);
				break;
			case 'h':
				Usage(argv[0]);
				return 1;
			default:
				break;
		}
	}

	std::ostream &Out = (OutFile.is_open()) ? OutFile : std::cout;

	Neutrino *TheNu_, *TheNuB;
	if (Left)
	{
		TheNu_ = new Neutrino(0, Neutrino::Dirac | Neutrino::Left );
		TheNuB = new Neutrino(0, Neutrino::Dirac | Neutrino::Left | Neutrino::Antiparticle);
	}
	else if (Right)
	{
		TheNu_ = new Neutrino(0, Neutrino::Dirac | Neutrino::Right );
		TheNuB = new Neutrino(0, Neutrino::Dirac | Neutrino::Right | Neutrino::Antiparticle);
	}
	else
	{
		TheNu_ = new Neutrino(0, Neutrino::Dirac | Neutrino::Unpolarised );
		TheNuB = new Neutrino(0, Neutrino::Dirac | Neutrino::Unpolarised | Neutrino::Antiparticle);
	}
	TheNu_->SetDecayChannel(Channel);
	TheNuB->SetDecayChannel(Channel);

	Detector *TheBox = new Detector(DetConfig);

	std::vector<char> vFlag;
	Out << "#Mass\t";
	if (UeFlag)
	{
		Out << "Ue\t";
		if (Efficiency)
			TheBox->SetEfficiency(Channel, Detector::E);
		vFlag.push_back('E');
	}
	if (UmFlag)
	{
		Out << "Um\t";
		if (Efficiency)
			TheBox->SetEfficiency(Channel, Detector::M);
		vFlag.push_back('M');
	}
	if (UtFlag)
	{
		Out << "Ut\t";
		if (Efficiency)
			TheBox->SetEfficiency(Channel, Detector::T);
		vFlag.push_back('T');
	}
	Out << "Events" << std::endl;

	Engine *TheEngine = new Engine(FluxConfig, 1, 1);	//creating 1FHC and 1RHC fluxedrivers
	TheEngine->BindNeutrino(TheNu_, Engine::FHC, 0);
	TheEngine->BindNeutrino(TheNuB, Engine::RHC, 0);

	unsigned int Grid = 250;
	unsigned int nD = vFlag.size();	//number of dimensions
	double Mass;
	std::cout << "Scanning over " << nD << " dimensions" << std::endl;

	std::vector<double> vSignal;	//summing over energy, array of Uus
	std::vector<std::vector<double> > vGridlU2;

	double lU2_A = -12.0;
	double lU2_B = - 6.0;
	double lU2_C = - 6.0;
	double lU2_D =   0.0;
	double lU2Step = (lU2End - lU2Start) / Grid;
	double Ue, Um, Ut;
	bool UeSet, UmSet, UtSet;

	std::vector<double> vlU2_A(nD), vlU2_B(nD), vlU2_C(nD), vlU2_D(nD);
	for (unsigned int i = 0; i < nD; ++i)
	{
		vlU2_A.at(i) = lU2_A;
		vlU2_B.at(i) = lU2_B;
		vlU2_C.at(i) = lU2_C;
		vlU2_D.at(i) = lU2_D;
	}	

	for (double logMass = -2.0; logMass < 0.3; logMass += 2.3/Grid)	//increase mass log
	{
		Mass = pow(10.0, logMass);
		std::cout << "Mass " << Mass << std::endl;

		TheNu_->SetMass(Mass);
		TheNuB->SetMass(Mass);

		if (!TheNu_->IsDecayAllowed() && !TheNuB->IsDecayAllowed())
			continue;

		TheEngine->MakeFlux();
		TheEngine->ScaleDetector(TheBox);

		bool SetGrid = false;
		vGridlU2.clear();

		vSignal.clear();
		vSignal.resize(pow(Grid, nD));	//number of Uus probing

		unsigned int g = 0;

		double A = 0, B = 0.5, C = 0.5, D = 1;
		while (A-B > 1e-6 && C-D > 1e-6)
		{
			SetMix(TheNu, TheNuB, vFlag, vlU2_A);
			A = TheEngine->DecayNumberIntegrated(TheBox, Efficiency) - Thr;

			SetMix(TheNu, TheNuB, vFlag, vlU2_B);
			B = TheEngine->DecayNumberIntegrated(TheBox, Efficiency) - Thr;

			SetMix(TheNu, TheNuB, vFlag, vlU2_C);
			C = TheEngine->DecayNumberIntegrated(TheBox, Efficiency) - Thr;

			SetMix(TheNu, TheNuB, vFlag, vlU2_D);
			D = TheEngine->DecayNumberIntegrated(TheBox, Efficiency) - Thr;

			if (A * B < 0)		//bot is in between here
			{
				std::vector<double> vlU2_M(nD);
				for (unsigned int i = 0; i < nD; ++i)
					vlU2_M.at(i) = (vlU2_A.at(i) + vlU2_B.at(i)) / 2.0;
				SetMix(TheNu, TheNuB, vFlag, vlU2_M);
				M = TheEngine->DecayNumberIntegrated(TheBox, Efficiency) - Thr;
			}
		}

		for (double Energy = Start; Energy < End; Energy += EnStep)
		{
			std::vector<double> vlU2(nD, lU2Start);
			TheNu_->SetEnergy(Energy);
			TheNuB->SetEnergy(Energy);

			bool Iter = true;
			unsigned int g = 0;
			while (Iter)
			{
				UeSet = UmSet = UtSet = false;
				for (unsigned int i = 0; i < nD; ++i)
				{
					double Uu = pow(10.0, 0.5 * vlU2.at(i));

					if (UeFlag && !UeSet)
					{
						TheNu_->SetMixings(Uu, TheNu_->Um(), TheNu_->Ut());
						TheNuB->SetMixings(Uu, TheNuB->Um(), TheNuB->Ut());
						UeSet = true;
					}
					else if (UmFlag && !UmSet)
					{
						TheNu_->SetMixings(TheNu_->Ue(), Uu, TheNu_->Ut());
						TheNuB->SetMixings(TheNuB->Ue(), Uu, TheNuB->Ut());
						UmSet = true;
					}
					else if (UtFlag && !UtSet)
					{
						TheNu_->SetMixings(TheNu_->Ue(), TheNu_->Um(), Uu);
						TheNuB->SetMixings(TheNuB->Ue(), TheNuB->Um(), Uu);
						UtSet = true;
					}
				}

				if (vGridlU2.size() < vSignal.size())
					vGridlU2.push_back(vlU2);

				vSignal.at(g++) += EnStep * TheEngine->DecayNumber(TheBox, Efficiency);

				///carry operations
				vlU2.at(0) += lU2Step;
				unsigned int c = 0;
				while (Iter && vlU2.at(c) >= lU2End - lU2Step)
				{
					vlU2.at(c)    = 0;
					vlU2.at(c++) += lU2Step;

					if (c < nD)
						Iter = true;
					else
						Iter = false;
				}
			}
		}
	
		for (unsigned int g = 0; g < vSignal.size(); ++g)
		{
			Out << Mass << "\t";
			for (unsigned int u = 0; u < nD; ++u)
			{
				double Uu2 = pow(10.0, vGridlU2.at(g).at(u));
				Out << Uu2 << "\t";
			}
			Out << vSignal.at(g) << std::endl;
		}
	}

	//OutFile->cd();
	//logCont->Write();
	//Contour->Write();
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
	std::cout <<"\n  -c,  --channel" << std::endl;
	std::cout << "\t\tDecay channel, defaul ALL" << std::endl;
	std::cout <<"\n  -E,  -M,  -T" << std::endl;
	std::cout << "\t\tSelect which mixing element" << std::endl;
	std::cout <<"\n  -h,  --help" << std::endl;
	std::cout << "\t\tPrint this message and exit" << std::endl;
}

void SetMix(Neutrino *Nu_, Neutrino *NuB, const &std::map<char, double> mlU2)
{
	std::map<char, double>::const_iterator iu;
	for (iu = mlU2.begin(); iu != mlU2.end(); ++iu)
	{
		double uu = pow(10, 0.5 * iu->second)
		if (iu->first == 'E')
		{
			TheNu_->SetMixings(uu, TheNu_->Um(), TheNu_->Ut());
			TheNuB->SetMixings(uu, TheNuB->Um(), TheNuB->Ut());
		}
		else if (iu->first == 'M')
		{
			TheNu_->SetMixings(TheNu_->Ue(), uu, TheNu_->Ut());
			TheNuB->SetMixings(TheNuB->Ue(), uu, TheNuB->Ut());
		}
		else if (iu->first == 'T')
		{
			TheNu_->SetMixings(TheNu_->Ue(), TheNu_->Um(), uu);
			TheNuB->SetMixings(TheNuB->Ue(), TheNuB->Um(), uu);
		}
	}
}
