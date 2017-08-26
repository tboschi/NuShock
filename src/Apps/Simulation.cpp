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
#include "Particle.h"

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

	TTree *Data = new TTree("Data", "Particle tracks");

	Data->Branch("ID", &ID, "iID/i");
	Data->Branch("Real", &Real, "fReal/D");

	Data->Branch("E_A", &EnergyA, "fEnergyA/D");
	Data->Branch("P_A", &MomentA, "fMomentA/D");
	Data->Branch("T_A", &TransvA, "fTransvA/D");
	Data->Branch("TheA", &ThetaA, "fThetaA/D");
	Data->Branch("PhiA", &PhiA, "fPhiA/D");
	Data->Branch("M_A", &MassA, "fMassA/D");

	Data->Branch("E_B", &EnergyB, "fEnergyB/D");
	Data->Branch("P_B", &MomentB, "fMomentB/D");
	Data->Branch("T_B", &TransvB, "fTransvB/D");
	Data->Branch("TheB", &ThetaB, "fThetaB/D");
	Data->Branch("PhiB", &PhiB, "fPhiB/D");
	Data->Branch("M_B", &MassB, "fMassB/D");

	Data->Branch("Angle", &Angle, "fAngle/D");

	Data->Branch("E_0", &Energy0, "fEnergy0/D");
	Data->Branch("P_0", &Moment0, "fMoment0/D");
	Data->Branch("T_0", &Transv0, "fTransv0/D");
	Data->Branch("The0", &Theta0, "fTheta0/D");
	Data->Branch("Phi0", &Phi0, "fPhi0/D");
	Data->Branch("M_0", &Mass0, "fMass0/D");

	while (ID < Nevent)
	{
		EvGen->SampleInDetector(1);

		Real = EvGen->GetEnergy();

		if (EvGen->EventKinematics())
		{
			Particle *ParticleA = EvGen->GetDecayProduct(0, 1);
			Particle *ParticleB = EvGen->GetDecayProduct(1, 1);

			EnergyA = ParticleA->E();
			MomentA = ParticleA->P();
			TransvA = ParticleA->Pt();
			ThetaA = ParticleA->Theta();
			PhiA = ParticleA->Phi();
			MassA = ParticleA->M();
		
			EnergyB = ParticleB->E();
			MomentB = ParticleB->P();
			TransvB = ParticleB->Pt();
			ThetaB = ParticleB->Theta();
			PhiB = ParticleB->Phi();
			MassB = ParticleB->M();
	
			Angle = ParticleA->Direction().Angle(ParticleB->Direction());
	
			TLorentzVector Reco = ParticleA->GetP4() + ParticleB->GetP4();
		
			Energy0 = Reco.E();
		        Moment0 = Reco.P();
		        Transv0 = Reco.Pt();
		        Theta0 = Reco.Theta();
		        Phi0 = Reco.Phi();
		        Mass0 = Reco.M();
	
			Data->Fill();

			std::cout << "Particle " << ++ID << std::endl;

			delete ParticleA, ParticleB;

			if (SaveMe++ > 10000)
			{
				OutFile->cd();
				Data->Write();
				SaveMe = 0;
			}
		}
	}

	OutFile->cd();
        Data->Write();

	/*
	for (unsigned int i = 0; i < Data->GetEntries(); ++i)
	{
		Data->GetEntry(i);

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
	*/

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
