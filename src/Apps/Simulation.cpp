#include <iostream>
#include <fstream>
#include <vector>
#include <cstring>
#include <getopt.h>

#include "TFile.h"
#include "TH2.h"

#include "Tools.h"
#include "EventGenerator.h"
#include "DecayRates.h"
#include "3Body.h"

void Usage(char* argv0);

int main(int argc, char** argv)
{
	const struct option longopts[] = 
	{
		{"smconfig", 	required_argument, 	0, 's'},
		{"detconfig", 	required_argument,	0, 'd'},
		{"fluxconfig", 	required_argument,	0, 'f'},
		{"channel", 	required_argument,	0, 'c'},
		{"output", 	required_argument,	0, 'o'},
		{"help", 	no_argument,		0, 'h'},
		{0,	0, 	0,	0},
	};

	int index; 
	int iarg = 0;
	opterr = 1;
	
	std::string SMConfig, DetConfig, FluxConfig;
	TFile *OutFile;

	std::string Channel = "ALL";
	
	while((iarg = getopt_long(argc,argv, "s:d:f:c:o:h", longopts, &index)) != -1)	
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
				Channel.assign(optarg);
				break;
			case 'o':
				OutFile = new TFile(optarg, "RECREATE");
				//OutFile = new TFile(optarg, "RECREATE");
				break;
			case 'h':
				Usage(argv[0]);
				return 1;
			default:
				break;
		}
	}

	//std::ostream &Out = (OutFile.is_open()) ? OutFile : std::cout;

	OutFile->cd();
	
	EventGenerator * EvGen = new EventGenerator(SMConfig, DetConfig, FluxConfig);

	EvGen->SetChannel("nEE");
	EvGen->SetEnergy(1);	//1GeV energy

	TH2D *hDalitz = new TH2D("dalitz", "Dalitz", 200,0,0.2, 200,0,0.2);
	TH1D *hEnergy = new TH1D("energy", "Energy", 200, 0, 1);

	int N = 0;
	while (N < 100000)
	{
		double part = EvGen->EventKinematics();
		if (part > 0)
		{
			TLorentzVector *p0 = EvGen->GetDecayProduct(0);
			TLorentzVector *p1 = EvGen->GetDecayProduct(1);
			TLorentzVector *p2 = EvGen->GetDecayProduct(2);

			TLorentzVector m01 = *p0 + *p1;
			//TLorentzVector m12 = *p1 + *p2;

			//hDalitz->Fill(m01.M2(), m12.M2());
			hEnergy->Fill(p1->E());
		
			++N;
		}
	}

	//hDalitz->Write();
	hEnergy->Write();

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
	std::cout <<"\n  -c,  --channel" << std::endl;
	std::cout << "\t\tDecay channel, defaul ALL" << std::endl;
	std::cout <<"\n  -o,  --output" << std::endl;
	std::cout << "\t\tOutput file" << std::endl;
	std::cout <<"\n  -h,  --help" << std::endl;
	std::cout << "\t\tPrint this message and exit" << std::endl;
}
