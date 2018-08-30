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
	std::string DetConfig, FluxConfig, FileName;
	std::ofstream OutFile;
	//TFile *OutFile;
	std::string Channel = "ALL";
	bool UeFlag = false, UmFlag = false, UtFlag = false;

	bool Left = false, Right = false;		//default unpolarised
	bool Particle = false, Antipart = false;	//default majorana

	bool Efficiency = false;
	double Thr = 2.44, Qct = 0.0;
	
	while((iarg = getopt_long(argc,argv, "d:f:c:o:WEMTLRAPt:q:h", longopts, &index)) != -1)
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
			case 't':
				Thr = std::strtod(optarg, NULL);
				break;
			case 'q':
				Qct = std::strtod(optarg, NULL);
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

	std::vector<double> vSignal;	//summing over energy, array of Uus
	std::vector<std::vector<double> > vGridlU2;

	double lU2Start = -12.0;
	double lU2End   = 0.0;
	double lU2Step  = (lU2End - lU2Start) / Grid;
	double Ue, Um, Ut;
	bool UeSet, UmSet, UtSet;

	for (double logMass = -2.0; logMass < 0.3; logMass += 2.3/Grid)	//increase mass log
	{
		Mass = pow(10.0, logMass);
		std::cout << "Mass " << Mass << std::endl;

		TheNu->SetMass(Mass);

		if (TheNu->IsDecayAllowed() &&
		    TheNu->IsProductionAllowed())
		{
			TheEngine->MakeFlux();
			TheEngine->ScaleDetector(TheBox);

			for (double lU2 = lU2Start; lU2 < lU2End; lU2 += lU2Step)
			{
				double Uu = pow(10.0, 0.5 * lU2);
				double Ue = 0.0, Um = 0.0, Ut = 0.0;

				for (unsigned int f = 0; f < vFlag.size(); ++f)
				{
					if (vFlag.at(f) == 'E')
						Ue = Uu;
					else if (vFlag.at(f) == 'M')
						Um = Uu;
					else if (vFlag.at(f) == 'T')
						Ut = Uu;
				}

				TheNu->SetMixings(Ue, Um, Ut);

				std::vector<double> vInt;
				double Signal = TheEngine->MakeSampler(TheBox, vInt);
				Out << Mass << "\t" << Uu*Uu << "\t" << Signal << std::endl;
			}
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

/*
			double Start, End;
			double EnStep = TheEngine->RangeWidth(Start, End);
			for (double Energy = Start; Energy < End; Energy += EnStep)
			{
				std::vector<double> vlU2(nD, lU2Start);
				TheNu_->SetEnergy(Energy);
				TheNuB->SetEnergy(Energy);

				bool Iter = true;
				unsigned int g = 0;

			//std::cout << ue << "\t" << um << "\t" << ut << std::endl;
			NN->SetMixings(ue, um, ut);
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

					vSignal.at(g++) += EnStep * TheEngine->DecayNumber(TheBox, Efficiency);

*/
