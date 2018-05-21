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

	int bins = 2000;

	TH1D *hMuonElec  = new TH1D("muonelec",  "Muon->Elec",  bins,0,2.0);
	TH1D *hMuonMuon  = new TH1D("muonmuon",  "Muon->Muon",  bins,0,2.0);
	TH1D *hKaonElec  = new TH1D("kaonelec",  "Kaon->Elec",  bins,0,2.0);
	TH1D *hKaonMuon  = new TH1D("kaonmuon",  "Kaon->Muon",  bins,0,2.0);
	TH1D *hKaon0Elec = new TH1D("kaon0elec", "Kaon0->Elec", bins,0,2.0);
	TH1D *hKaon0Muon = new TH1D("kaon0muon", "Kaon0->Muon", bins,0,2.0);
	TH1D *hTauEElec  = new TH1D("taueelec",  "TauE->Elec",  bins,0,2.0);
	TH1D *hTauETau   = new TH1D("tauetau",   "TauE->Tau",   bins,0,2.0);
	TH1D *hTauMMuon  = new TH1D("taummuon",  "TauM->Muon",  bins,0,2.0);
	TH1D *hTauMTau   = new TH1D("taumtau",   "TauM->Tau",   bins,0,2.0);

	double a,b,dx;
	ThreeBody *DecayMuon  = new ThreeBody("Muon",  0.0, 1.0, 1.0, 1.0);
	ThreeBody *DecayKaon0 = new ThreeBody("Kaon0", 0.0, 1.0, 1.0, 1.0);
	ThreeBody *DecayKaon  = new ThreeBody("Kaon",  0.0, 1.0, 1.0, 1.0);
	ThreeBody *DecayTauE  = new ThreeBody("TauE",  0.0, 1.0, 1.0, 1.0);
	ThreeBody *DecayTauM  = new ThreeBody("TauM",  0.0, 1.0, 1.0, 1.0);

	DecayMuon->ElectronChannel();
	double gME = DecayMuon->Gamma();
	DecayMuon->MuonChannel();
	double gMM = DecayMuon->Gamma();

	DecayKaon->ElectronChannel();
	double gKE = DecayKaon->Gamma();
	DecayKaon->MuonChannel();
	double gKM = DecayKaon->Gamma();

	DecayKaon0->ElectronChannel();
	double gK0E = DecayKaon0->Gamma();
	DecayKaon0->MuonChannel();
	double gK0M = DecayKaon0->Gamma();

	DecayTauE->ElectronChannel();
	double gTEE = DecayTauE->Gamma();
	DecayTauE->TauChannel();
	double gTET = DecayTauE->Gamma();

	DecayTauM->MuonChannel();
	double gTMM = DecayTauM->Gamma();
	DecayTauM->TauChannel();
	double gTMT = DecayTauM->Gamma();

	for (double MS = 0.0; MS <= 2.0; MS += 0.001)
	{
		DecayMuon->SetSterileMass(MS);
		DecayKaon->SetSterileMass(MS);
		DecayKaon0->SetSterileMass(MS);
		DecayTauE->SetSterileMass(MS);
		DecayTauM->SetSterileMass(MS);

		Out << MS << "\t";

		DecayMuon->ElectronChannel();
		hMuonElec->Fill(MS, DecayMuon->Gamma()/gME);
		Out << DecayMuon->Gamma()/gME << "\t";
		DecayMuon->MuonChannel();
		hMuonMuon->Fill(MS, DecayMuon->Gamma()/gMM);
		Out << DecayMuon->Gamma()/gMM << "\t";

		DecayKaon->ElectronChannel();
		hKaonElec->Fill(MS, DecayKaon->Gamma()/gKE);
		Out << DecayKaon->Gamma()/gKE << "\t";
		DecayKaon->MuonChannel();
		hKaonMuon->Fill(MS, DecayKaon->Gamma()/gKM);
		Out << DecayKaon->Gamma()/gKM << "\t";

		DecayKaon0->ElectronChannel();
		hKaon0Elec->Fill(MS, DecayKaon0->Gamma()/gK0E);
		Out << DecayKaon0->Gamma()/gK0E << "\t";
		DecayKaon0->MuonChannel();
		hKaon0Muon->Fill(MS, DecayKaon0->Gamma()/gK0M);
		Out << DecayKaon0->Gamma()/gK0M << "\t";

		DecayTauE->ElectronChannel();
		hTauEElec->Fill(MS, DecayTauE->Gamma()/gTEE);
		Out << DecayTauE->Gamma()/gTEE << "\t";
		DecayTauE->TauChannel();
		hTauETau->Fill(MS, DecayTauE->Gamma()/gTET);
		Out << DecayTauE->Gamma()/gTET << "\t";

		DecayTauM->MuonChannel();
		hTauMMuon->Fill(MS, DecayTauM->Gamma()/gTMM);
		Out << DecayTauM->Gamma()/gTMM << "\t";
		DecayTauM->TauChannel();
		hTauMTau->Fill(MS, DecayTauM->Gamma()/gTMT);
		Out << DecayTauM->Gamma()/gTMT << "\t";

		Out << std::endl;
	}

        hMuonElec->Write();
	hMuonMuon->Write();
        hKaonElec->Write();
        hKaonMuon->Write();
        hKaon0Elec->Write();
        hKaon0Muon->Write();
        hTauEElec->Write();
	hTauETau->Write();
        hTauMMuon->Write();
	hTauMTau->Write();

	OutFile->Close();

	return 0;
}
