#include <iostream>
#include <fstream>
#include <cstring>
#include <getopt.h>

#include "EventGenerator.h"
#include "FluxDriver.h"
#include "DecayRates.h"
#include "Detector.h"

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
		{"lambda", 	required_argument,	0, 'L'},
		{"wilson", 	required_argument,	0, 'c'},
		{"scattering", 	no_argument,		0, 'S'},
		{"decay", 	no_argument,		0, 'D'},
		{"eft", 	no_argument,		0, 'E'},
		{"output", 	required_argument,	0, 'o'},
		{"help", 	no_argument,	 	0, 'h'},
		{0,	0, 	0,	0},
	};

	int index; 
	int iarg = 0;
	opterr = 1;
	
	//Initialize variables
	std::string SMConfig, DetConfig;
	std::string FluxConfig, EffFile;
	std::ofstream OutFile;
	//TFile *OutFile;
	double Lambda = 1.0;
	double c_W = 1.0;
	bool DIF = true;
	bool EFT = false;
	
	while((iarg = getopt_long(argc,argv, "s:d:f:c:L:SDEo:h", longopts, &index)) != -1)
	{
		switch(iarg)
		{
			case 's':
				SMConfig.assign(optarg);
				break;
			case 'd':
				DetConfig.assign(optarg);
				break;
			case 'f':
				FluxConfig.assign(optarg);
				break;
			case 'c':
				c_W = strtod(optarg, NULL);
				break;
			case 'L':
				Lambda = strtod(optarg, NULL);
				break;
			case 'o':
				OutFile.open(optarg);
				break;
			case 'S':
				DIF = false;
				break;
			case 'D':
				DIF = true;
				break;
			case 'E':
				EFT = true;
				break;
			case 'h':
				Usage(argv[0]);
				return 1;
			default:
				break;
		}
	}

	std::ostream &Out = (OutFile.is_open()) ? OutFile : std::cout;

	//TH2D * logCont = new TH2D ("logcont", "Above threshold", 100, -2.0, 0.0, 100, -10.0, -4.0);
	//TH2D * Contour = new TH2D ("contour", "Above threshold", 100, 0.01, 1.0, 100, 1.0e-10, 1.0e-4);
	EventGenerator * EvGen = new EventGenerator(SMConfig, DetConfig, FluxConfig);
	
	EvGen->SetChannel("ALL");
	EvGen->SetMass(0);
	EvGen->SetUe(0);
	EvGen->SetUm(0);
	EvGen->SetUt(0);
	
	double Mass, Uu, Gamma;
	std::vector<double> vSignal;	//summing over energy, array of Uus
	int Grid = 100;
	
	for (double logMass = -2.0; logMass < 0.75; logMass += 2.75/Grid)	//increase mass log
	{
		Mass = pow(10.0, logMass);
		std::cout << "Mass " << Mass << std::endl;

		EvGen->SetMass(Mass);
		EvGen->MakeFlux(DIF);

		vSignal.clear();
		vSignal.resize(Grid);	//number of Uus probing

		double Start, End;
		double EnStep = EvGen->GetRange(Start, End)/EvGen->GetBinNumber();
		for (double Energy = Start; Energy < End; Energy += EnStep)
		{
			unsigned int i = 0;
			for (double logUu2 = -10.0; logUu2+1e-6 < 0.0; logUu2 += 10.0/Grid, ++i)	//increase Uu logarithmically
			{
				Uu = pow(10.0, 0.5*logUu2);
				EvGen->SetUe(Uu/sqrt(2.0));
				EvGen->SetUm(Uu/sqrt(2.0));

				if (EFT)
					Gamma = 1e-17 * pow(Mass, 5) * ((1 + c_W/Lambda * (1 + c_W/Lambda)) * Uu*Uu + c_W/Lambda/Lambda * (1 + c_W/Lambda/Lambda));
				else
					Gamma = 1e-17 * pow(Mass, 5) * Uu*Uu;

				EvGen->SetUserData(Gamma);
	
				if (DIF)
					vSignal.at(i) += EnStep * EvGen->DecayNumber(Energy);
				else
					vSignal.at(i) += EnStep * EvGen->ScatterNumber(Energy);
			}
		}

		unsigned int j = 0;
		for (double logUu2 = -10.0; logUu2+1e-6 < 0.0; logUu2 += 10.0/Grid, ++j)	//increase Uu logarithmically
		{
			Uu = pow(10.0, 0.5*logUu2);
			if (EFT)
				Gamma = 1e-17 * pow(Mass, 5) * ((1 + c_W/Lambda * (1 + c_W/Lambda)) * Uu*Uu + c_W/Lambda/Lambda * (1 + c_W/Lambda/Lambda));
			else
				Gamma = 1e-17 * pow(Mass, 5) * Uu*Uu;

			Out << Mass << "\t" << Uu*Uu << "\t" << vSignal.at(j) << "\t" << Const::fhBar / Gamma << "\t" << std::endl;
		}

		Out << std::endl;
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
	std::cout <<"\n  -c,  --channel" << std::endl;
	std::cout << "\t\tDecay channel, defaul ALL" << std::endl;
	std::cout <<"\n  -E,  -M,  -T" << std::endl;
	std::cout << "\t\tSelect which mixing element" << std::endl;
	std::cout <<"\n  -h,  --help" << std::endl;
	std::cout << "\t\tPrint this message and exit" << std::endl;
}
