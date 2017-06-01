#include <iostream>
#include <fstream>
#include <vector>
#include <cstring>
#include <getopt.h>

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"

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
	
	EventGenerator * EvGen = new EventGenerator(SMConfig, DetConfig, FluxConfig);

	EvGen->SetChannel(Channel);
	EvGen->MakeSterileFlux(1);

	//TH2D *hDalitz = new TH2D("dalitz", "Dalitz", 200,0,0.2, 200,0,0.2);
	TH1D *hEnergy = new TH1D("energy", "Energy", 500, 0, 20);
	TH1D *hAngleS = new TH1D("angles", "Angle Separation", 500, -180, 180);
	TH1D *hAngleN = new TH1D("anglen", "Angle N", 500, -180, 180);
	TH1D *hAngle0 = new TH1D("angle0", "Angle p0", 500, -180, 180);
	TH1D *hAngle1 = new TH1D("angle1", "Angle p1", 500, -180, 180);
	TH1D *hIMassN = new TH1D("imassn", "Invariant Mass of N", 500, 0, 0.5);
	TH1D *hIMass0 = new TH1D("imass0", "Invariant Mass of p0", 500, 0, 0.5);
	TH1D *hIMass1 = new TH1D("imass1", "Invariant Mass of p1", 500, 0, 0.5);

	int N = 0;
	while (N < 100000)
	{
		EvGen->SampleEnergy();

		double part = EvGen->EventKinematics();
		if (part > 0)
		{
			TLorentzVector *p0 = EvGen->GetDecayProduct(0, 1);
			TLorentzVector *p1 = EvGen->GetDecayProduct(1, 1);
			//TLorentzVector *p2 = EvGen->GetDecayProduct(2);

			TLorentzVector *m01 = new TLorentzVector(*p0 + *p1);
			//TLorentzVector m12 = *p1 + *p2;

			//hDalitz->Fill(m01.M2(), m12.M2());
			hEnergy->Fill(m01->E());
			hAngleS->Fill(Const::fDeg * p0->Angle(p1->Vect()) );	//track separation
			hAngleN->Fill(Const::fDeg * m01->Theta() );	//along z-axis
			hAngle0->Fill(Const::fDeg * p0->Theta() );	//along z-axis
			hAngle1->Fill(Const::fDeg * p1->Theta() );	//along z-axis
			hIMassN->Fill(m01->M());	//decaying mass
			hIMass0->Fill(p0->M());		//p0 mass
			hIMass1->Fill(p1->M());		//p1 mass
		
			++N;
		}
	}

	OutFile->cd();

	//EvGen->GetFluxDriverPtr()->GetTotal()->Write();
	hEnergy->Write();
        hAngleS->Write();
        hAngleN->Write();
        hAngle0->Write();
        hAngle1->Write();
        hIMassN->Write();
        hIMass0->Write();
        hIMass1->Write();

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
