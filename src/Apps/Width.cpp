#include <iostream>
#include <fstream>
#include <vector>
#include <cstring>
#include <sstream>
#include <getopt.h>

#include "Tools.h"
#include "ThreeBody.h"

#include "TFile.h"
#include "TH1D.h"

int main(int argc, char** argv)
{
	const struct option longopts[] = 
	{
		{"root", 	required_argument, 	0, 'r'},
		{"output", 	required_argument, 	0, 'o'},
		{"smconfig", 	required_argument, 	0, 's'},
		{"muon", 	no_argument,	 	0, 'm'},
		{"kaon", 	no_argument, 		0, 'k'},
		{"kaon0", 	no_argument,	 	0, 'K'},
		{"help", 	no_argument,		0, 'h'},
		{0,	0, 	0,	0},
	};

	int index; 
	int iarg = 0;
	opterr = 1;
	
	std::ifstream ConfigFile;
	std::ofstream Out;
	TFile *OutFile;
	bool IsElectron = false;
	bool IsMuon = false;
	bool MuonFlag = false;
	bool KaonFlag = false;
	bool Kaon0Flag = false;
	
	while((iarg = getopt_long(argc,argv, "r:o:s:EMmkKh", longopts, &index)) != -1)	
	{
		switch(iarg)
		{
			case 'r':
				OutFile = new TFile(optarg, "RECREATE");
				break;
			case 'o':
				Out.open(optarg);
				break;
			case 'h':
				std::cout << "Compute decay width fro three body decay" << std::endl;
				std::cout << "Usage : " << std::endl;
				std::cout << "Width [OPTIONS]" << std::endl;
				std::cout <<"\n  -r,  --root" << std::endl;
				std::cout << "\t\tOutput file to save plot (ROOT)" << std::endl;
				std::cout <<"\n  -o,  --output" << std::endl;
				std::cout << "\t\tOutput file to save plot (text)" << std::endl;
				std::cout <<"\n  -h,  --help" << std::endl;
				std::cout << "\t\tPrint this message and exit" << std::endl;
				return 1;
			default:
				break;
		}
	}

	//std::ostream &Out = (OutFile.is_open()) ? OutFile : std::cout;

	double M_MSterile; 
	double U_e;
	double U_m;
	double U_t;

	std::string Line, Key, Name;
	std::stringstream ssL;
	double Element;

	/*
	while (std::getline(ConfigFile, Line))
	{
		if (Line[0] == '#') continue;

		ssL.str("");
		ssL.clear();
		ssL << Line;

		ssL >> Key >> Element;
		if (Key == "M_Sterile") M_MSterile = Element;
		if (Key == "U_e") U_e = Element;
		if (Key == "U_m") U_m = Element;
		if (Key == "U_t") U_t = Element;
	}
	ConfigFile.close();
	*/

	int bins = 500;

	TH1D *hMuonMuon = new TH1D("muonmuon", "Muon->Muon", bins,0,0.500);
	TH1D *hMuonElec = new TH1D("muonelec", "Muon->Elec", bins,0,0.500);
	TH1D *hKaonMuon = new TH1D("kaonmuon", "Kaon->Muon", bins,0,0.500);
	TH1D *hKaonElec = new TH1D("kaonelec", "Kaon->Elec", bins,0,0.500);
	TH1D *hKaon0Muon = new TH1D("kaon0muon", "Kaon0->Muon", bins,0,0.500);
	TH1D *hKaon0Elec = new TH1D("kaon0elec", "Kaon0->Elec", bins,0,0.500);

	double a,b,dx;
	ThreeBody *Decay = new ThreeBody("Muon", 0.0, 1.0, 1.0, 1.0);
	Decay->MuonChannel();
	double gMM = Decay->Gamma();
	Decay->ElectronChannel();
	double gME = Decay->Gamma();

	Decay->SetParent("Kaon");
	Decay->MuonChannel();
	double gKM = Decay->Gamma();
	Decay->ElectronChannel();
	double gKE = Decay->Gamma();

	Decay->SetParent("Kaon0");
	Decay->MuonChannel();
	double gK0M = Decay->Gamma();
	Decay->ElectronChannel();
	double gK0E = Decay->Gamma();

	for (double MS = 0.0; MS <= 0.5; MS += 0.5/bins)
	{
		std::cout << "Mass " << MS << std::endl;
		Decay->SetSterileMass(MS);

		Out << MS << "\t";

		Decay->SetParent("Muon");
		Decay->MuonChannel();
		hMuonMuon->Fill(MS, Decay->Gamma()/gMM);
		Out << Decay->Gamma()/gMM << "\t";
		Decay->ElectronChannel();
		hMuonElec->Fill(MS, Decay->Gamma()/gME);
		Out << Decay->Gamma()/gME << "\t";

		Decay->SetParent("Kaon");
		Decay->MuonChannel();
		hKaonMuon->Fill(MS, Decay->Gamma()/gKM);
		Out << Decay->Gamma()/gKM << "\t";
		Decay->ElectronChannel();
		hKaonElec->Fill(MS, Decay->Gamma()/gKE);
		Out << Decay->Gamma()/gKE << "\t";

		Decay->SetParent("Kaon0");
		Decay->MuonChannel();
		hKaon0Muon->Fill(MS, Decay->Gamma()/gK0M);
		Out << Decay->Gamma()/gK0M << "\t";
		Decay->ElectronChannel();
		hKaon0Elec->Fill(MS, Decay->Gamma()/gK0E);
		Out << Decay->Gamma()/gK0E << std::endl;
	}

	hMuonMuon->Write();
        hMuonElec->Write();
        hKaonMuon->Write();
        hKaonElec->Write();
        hKaon0Muon->Write();
        hKaon0Elec->Write();

	OutFile->Close();

	return 0;
}
