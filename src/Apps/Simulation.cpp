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
		{"nevent", 	required_argument,	0, 'n'},
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
	unsigned int Nevent;

	std::string Channel = "ALL";
	
	while((iarg = getopt_long(argc,argv, "s:d:f:n:c:o:h", longopts, &index)) != -1)	
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
			case 'n':
				Nevent = strtod(optarg, NULL);
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
	EvGen->MakeInDetector();


	//TH2D *hDalitz = new TH2D("dalitz", "Dalitz", 200,0,0.2, 200,0,0.2);
	
	unsigned int ID = 0;
	unsigned int SaveMe = 0;
	double Real;

	double EnergyA, EnergyB, EnergyC, Energy0;
	double energyA, energyB;
	double MomentA, MomentB, MomentC, Moment0;
	double momentA, momentB;
	double TransvA, TransvB, TransvC, Transv0;
	double transvA, transvB;
	double ThetaA, ThetaB, ThetaC, Theta0;
	double thetaA, thetaB;
	double PhiA, PhiB, PhiC, Phi0;
	double phiA, phiB;
	double MassA, MassB, MassC, Mass0;
	double massA, massB;
	double Angle;	//separation between A & B
	double angle;	//separation between A & B

	TTree *tTrack = new TTree("tTrack", "Particle tracks");

	tTrack->Branch("ID", &ID, "iID/i");
	tTrack->Branch("Real", &Real, "fReal/D");

	tTrack->Branch("E_A", &EnergyA, "fEnergyA/D");
	tTrack->Branch("P_A", &MomentA, "fMomentA/D");
	tTrack->Branch("T_A", &TransvA, "fTransvA/D");
	tTrack->Branch("TheA", &ThetaA, "fThetaA/D");
	tTrack->Branch("PhiA", &PhiA, "fPhiA/D");
	tTrack->Branch("M_A", &MassA, "fMassA/D");

	tTrack->Branch("E_B", &EnergyB, "fEnergyB/D");
	tTrack->Branch("P_B", &MomentB, "fMomentB/D");
	tTrack->Branch("T_B", &TransvB, "fTransvB/D");
	tTrack->Branch("TheB", &ThetaB, "fThetaB/D");
	tTrack->Branch("PhiB", &PhiB, "fPhiB/D");
	tTrack->Branch("M_B", &MassB, "fMassB/D");

	//unboosted
	tTrack->Branch("e_A", &energyA, "fenergyA/D");
	tTrack->Branch("p_A", &momentA, "fmomentA/D");
	tTrack->Branch("t_A", &transvA, "ftransvA/D");
	tTrack->Branch("theA", &thetaA, "fthetaA/D");
	tTrack->Branch("phiA", &phiA, "fphiA/D");
	tTrack->Branch("m_A", &massA, "fmassA/D");

	tTrack->Branch("e_B", &energyB, "fenergyB/D");
	tTrack->Branch("p_B", &momentB, "fmomentB/D");
	tTrack->Branch("t_B", &transvB, "ftransvB/D");
	tTrack->Branch("theB", &thetaB, "fthetaB/D");
	tTrack->Branch("phiB", &phiB, "fphiB/D");
	tTrack->Branch("m_B", &massB, "fmassB/D");

	//tTrack->Branch("EnergyC", &EnergyC, "fEnergyC/D");
	//tTrack->Branch("MomentC", &MomentC, "fMomentC/D");
	//tTrack->Branch("TransvC", &TransvC, "fTransvC/D");
	//tTrack->Branch("ThetaC", &ThetaC, "fThetaC/D");
	//tTrack->Branch("PhiC", &PhiC, "fPhiC/D");
	//tTrack->Branch("MassC", &MassC, "fMassC/D");

	tTrack->Branch("Angle", &Angle, "fAngle/D");
	//tTrack->Branch("angle", &angle, "fangle/D");

	/*
	tTrack->Branch("E_0", &Energy0, "fEnergy0/D");
	tTrack->Branch("P_0", &Moment0, "fMoment0/D");
	tTrack->Branch("T_0", &Transv0, "fTransv0/D");
	tTrack->Branch("The0", &Theta0, "fTheta0/D");
	tTrack->Branch("Phi0", &Phi0, "fPhi0/D");
	tTrack->Branch("M_0", &Mass0, "fMass0/D");
	*/

	TH1D *E_All = new TH1D("E_All", "Neutrino energy before cut", 100, 0, 20);
	TH1D *E_Cut = new TH1D("E_Cut", "Neutrino energy after cut", 100, 0, 20);

	TH1D *InvMass = new TH1D("E_Cut", "Neutrino energy after cut", 500, 0, 1);	//Invariant mass for 2body decays

	//TH1D *hEnergy = new TH1D("energy", "Energy", 500, 0, 20);
	//TH1D *hAngleS = new TH1D("angles", "Angle Separation", 500, -180, 180);
	//TH1D *hAngleN = new TH1D("anglen", "Angle N", 500, -180, 180);
	//TH1D *hAngle0 = new TH1D("angle0", "Angle p0", 500, -180, 180);
	//TH1D *hAngle1 = new TH1D("angle1", "Angle p1", 500, -180, 180);
	//TH1D *hIMassN = new TH1D("imassn", "Invariant Mass of N", 500, 0, 0.5);
	//TH1D *hIMass0 = new TH1D("imass0", "Invariant Mass of p0", 500, 0, 0.5);
	//TH1D *hIMass1 = new TH1D("imass1", "Invariant Mass of p1", 500, 0, 0.5);

	double Grad = Const::fPi/180.0;

	//EvGen->SetEnergy(1);
	//fill a tree with simulated neutrinos
	while (ID < Nevent)
	{
		EvGen->SampleInDetector(1);

		Real = EvGen->GetEnergy();

		if (EvGen->EventKinematics())
		{
			//std::cout << "h1" << std::endl;
			TLorentzVector *pA = EvGen->GetDecayProduct(0, 1);
			TLorentzVector *pB = EvGen->GetDecayProduct(1, 1);
		//	TLorentzVector *pC = EvGen->GetDecayProduct(2, 1);
			TLorentzVector *p0 = new TLorentzVector(*pA + *pB);

			EnergyA = pA->E();
			MomentA = pA->P();
			TransvA = pA->Pt();
			ThetaA = pA->Theta();
			PhiA = pA->Phi();
			MassA = pA->M();

			EnergyB = pB->E();
			MomentB = pB->P();
			TransvB = pB->Pt();
			ThetaB = pB->Theta();
			PhiB = pB->Phi();
			MassB = pB->M();

		//	EnergyC = pC->E();
		//	MomentC = pC->P();
		//	TransvC = pC->Pt();
		//	ThetaC = pC->Theta();
		//	PhiC = pC->Phi();
		//	MassC = pC->M();

			Angle = pA->Angle(pB->Vect());

			Energy0 = p0->E();
			Moment0 = p0->P();
			Transv0 = p0->Pt();
			Theta0 = p0->Theta();
			Phi0 = p0->Phi();
			Mass0 = p0->M();

			//unboost
			/*
			pA->Boost(- p0->Vect());
			pB->Boost(- p0->Vect());
			angle = pA->Angle(pB->Vect());
	
			energyA = pA->E();
			momentA = pA->P();
			transvA = pA->Pt();
			thetaA = pA->Theta();
			phiA = pA->Phi();
			massA = pA->M();

			energyB = pB->E();
			momentB = pB->P();
			transvB = pB->Pt();
			thetaB = pB->Theta();
			phiB = pB->Phi();
			massB = pB->M();

			angle = pA->Angle(pB->Vect());
			*/

		//	EnergyC = pC->E();
			
			E_All->Fill(Energy0);
			InvMass->Fill(Mass0);

			tTrack->Fill();

			std::cout << "Particle " << ++ID << std::endl;

			delete p0;

			if (SaveMe++ > 1000)
			{
				OutFile->cd();
				tTrack->Write();
				SaveMe = 0;
			}
		}
	}

	OutFile->cd();
        tTrack->Write();

	for (unsigned int i = 0; i < tTrack->GetEntries(); ++i)
	{
		tTrack->GetEntry(i);

		if (Angle < 40*Grad &&
		    ThetaA < 80*Grad &&
		    ThetaB < 80*Grad && 
		    Theta0 < 4*Grad &&
		    Mass0 < InvMass->GetMean()+3*InvMass->GetRMS() &&
		    Mass0 > InvMass->GetMean()-3*InvMass->GetRMS() )
			E_Cut->Fill(Energy0);
	}

	for (unsigned int Bin = 0; Bin < E_All->GetNbinsX(); ++Bin)
	{
		Eff = E_Cut->GetBinContent(Bin)/E_All->GetBinContent(Bin);
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
	std::cout <<"\n  -c,  --channel" << std::endl;
	std::cout << "\t\tDecay channel, defaul ALL" << std::endl;
	std::cout <<"\n  -o,  --output" << std::endl;
	std::cout << "\t\tOutput file" << std::endl;
	std::cout <<"\n  -h,  --help" << std::endl;
	std::cout << "\t\tPrint this message and exit" << std::endl;
}
