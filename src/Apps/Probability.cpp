#include <iostream>
#include <fstream>
#include <cstring>
#include <getopt.h>

#include "EventGenerator.h"
#include "FluxDriver.h"
#include "DecayRates.h"
#include "Detector.h"

#include "TH1D.h"

void Usage(char* argv0);

int main(int argc, char** argv)
{
	const struct option longopts[] = 
	{
		{"smconfig", 	required_argument, 	0, 's'},
		{"detconfig", 	required_argument,	0, 'd'},
		{"fluxconfig", 	required_argument,	0, 'f'},
		{"output", 	required_argument,	0, 'o'},
		{"help", 	no_argument,	 	0, 'h'},
		{0,	0, 	0,	0},
	};

	int index; 
	int iarg = 0;
	opterr = 1;
	
	//Initialize variables
	std::string SMConfig, DetConfig;
	std::string FluxConfig;
	TFile *OutFile;
	
	while((iarg = getopt_long(argc,argv, "s:d:f:o:h", longopts, &index)) != -1)
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
			case 'o':
				OutFile = new TFile(optarg, "RECREATE");
				break;
			case 'h':
				Usage(argv[0]);
				return 1;
				return 1;
			default:
				break;
		}
	}

	//To have multiple output, handled by usage
//	std::ostream &Out = (OutFile.is_open()) ? OutFile : std::cout;

	TH1D *Fake = new TH1D("fake", "Fake", 200, 0, 20);
	TH1D *True = new TH1D("true", "True", 200, 0, 20);
	TH1D *Prob = new TH1D("prob", "Prob", 200, 0, 20);
	TH1D *Flux = new TH1D("flux", "Flux", 200, 0, 20);
	EventGenerator * EvGen = new EventGenerator(SMConfig, DetConfig, FluxConfig);
	
	EvGen->SetChannel("ALL");
	EvGen->MakeSterileFlux();

	std::cout << "Mass " << EvGen->GetMass() << std::endl;
	std::cout << "UUU " << EvGen->GetUe() << "\t" << EvGen->GetUm() << "\t" <<  EvGen->GetUt() << std::endl;
	for (double Energy = 0.05; Energy < 20.0; Energy += 0.1)	//increase mass
	{
		EvGen->SetEnergy(Energy);
		//std::cout << Energy << "\t" << EvGen->FluxIntensity() << "\t" << EvGen->EventProbability() << std::endl;
		Flux->Fill(Energy, EvGen->FluxIntensity());
		Prob->Fill(Energy, EvGen->EventProbability());
		Fake->Fill(Energy, EvGen->FluxIntensity()*EvGen->EventProbability());
		for (int i = 0; i < 1000000; ++i)
		{
			if (EvGen->EventInDetector())
				True->Fill(Energy);
		}
	}

	OutFile->cd();
	Prob->Write();
	Flux->Write();
	Fake->Write();
	True->Write();
	OutFile->Close();

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
	std::cout <<"\n  -h,  --help" << std::endl;
	std::cout << "\t\tPrint this message and exit" << std::endl;
}
