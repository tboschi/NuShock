#include <iostream>
#include <fstream>
#include <vector>
#include <cstring>
#include <getopt.h>

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TRandom3.h"

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
		{"number", 	required_argument,	0, 'n'},
		{"channel", 	required_argument,	0, 'c'},
		{"output", 	required_argument,	0, 'o'},
		{"mass", 	required_argument,	0, 'm'},
		{"ueratio", 	required_argument,	0, 'E'},
		{"umratio", 	required_argument,	0, 'M'},
		{"help", 	no_argument,		0, 'h'},
		{0,	0, 	0,	0},
	};

	int index; 
	int iarg = 0;
	opterr = 1;
	
	std::string SMConfig, DetConfig, FluxConfig;
	TFile *OutFile;
	unsigned int Nevent;
	double Mass = 0.0;
	double UeF = 0.0, UmF = 0.0;

	std::string Channel = "ALL";
	
	while((iarg = getopt_long(argc,argv, "s:d:f:m:n:c:o:E:M:h", longopts, &index)) != -1)	
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
			case 'm':
				Mass = strtod(optarg, NULL);
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
			case 'M':
				UmF = strtod(optarg, NULL);
				break;
			case 'E':
				UeF = strtod(optarg, NULL);
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

	if (Mass)
		EvGen->SetMass(Mass);
	UeF *= 1e-5/sqrt(UeF*UeF + UmF*UmF);  //sum^2 is 1e-10
	UmF *= 1e-5/sqrt(UeF*UeF + UmF*UmF);  
	EvGen->SetUe(UeF);
	EvGen->SetUe(UmF);

	EvGen->SetChannel(Channel);
	EvGen->MakeFlux(1);
	EvGen->MakeSampler();

	//TH2D *hDalitz = new TH2D("dalitz", "Dalitz", 200,i0,0.2, 200,0,0.2);
	
	unsigned int ID = 0;
	unsigned int ND = 0;
	unsigned int SaveMe = 0;

	double Real, Delay;

	double EnergyA, EnergyB, EnergyC, Energy0;
	double MomentA, MomentB, MomentC, Moment0;
	double TransvA, TransvB, TransvC, Transv0;
	double ThetaA, ThetaB, ThetaC, Theta0;
	double PhiA, PhiB, PhiC, Phi0;
	double MassA, MassB, MassC, Mass0;
	double InA, OutA, InB, OutB;
	double Angle;	//separation between A & B

	TRandom3 *Rand = new TRandom3(0);
	TTree *Data = new TTree("Data", "Particle tracks");

	Data->Branch("ID", &ID, "iID/i");
	Data->Branch("Real", &Real, "fReal/D");
	Data->Branch("Delay", &Delay, "fDelay/D");

	/*
	Data->Branch("E_A", &EnergyA, "fEnergyA/D");
	Data->Branch("P_A", &MomentA, "fMomentA/D");
	Data->Branch("T_A", &TransvA, "fTransvA/D");
	Data->Branch("TheA", &ThetaA, "fThetaA/D");
	Data->Branch("PhiA", &PhiA, "fPhiA/D");
	Data->Branch("M_A", &MassA, "fMassA/D");
	Data->Branch("In_A", &InA, "fInA/D");
	Data->Branch("Out_A", &OutA, "fOutA/D");

	Data->Branch("E_B", &EnergyB, "fEnergyB/D");
	Data->Branch("P_B", &MomentB, "fMomentB/D");
	Data->Branch("T_B", &TransvB, "fTransvB/D");
	Data->Branch("TheB", &ThetaB, "fThetaB/D");
	Data->Branch("PhiB", &PhiB, "fPhiB/D");
	Data->Branch("M_B", &MassB, "fMassB/D");
	Data->Branch("In_B", &InB, "fInB/D");
	Data->Branch("Out_B", &OutB, "fOutB/D");

	Data->Branch("Angle", &Angle, "fAngle/D");

	Data->Branch("E_0", &Energy0, "fEnergy0/D");
	Data->Branch("P_0", &Moment0, "fMoment0/D");
	Data->Branch("T_0", &Transv0, "fTransv0/D");
	Data->Branch("The0", &Theta0, "fTheta0/D");
	Data->Branch("Phi0", &Phi0, "fPhi0/D");
	Data->Branch("M_0", &Mass0, "fMass0/D");
	*/

	unsigned int nBunch = 80;	//80 bunhces 
	double tBunch = 20e-9;	//20ns between bunches
	unsigned int iBunch = 0;
	double Beta = 1.0;
	double Lc = EvGen->GetDetectorPtr()->GetElement("Baseline")/Const::fC;
	while (ID < Nevent)
	{
		Real = EvGen->SampleEnergy(1);

		iBunch = Rand->Integer(nBunch);
		if (Rand->Rndm() < 2.99)
			Beta = sqrt(1 - Mass*Mass / Real/Real);
		else
			Beta = 1.0;

		//Delay = Lc * (1.0/Beta - 1) + iBunch*tBunch;
		Delay = Rand->Gaus(Lc * (1.0/Beta - 1) + iBunch*tBunch, 1e-9);

		Data->Fill();
		++ID;

		/*
		if (EvGen->EventKinematics())
		{
			EvGen->GeneratePosition();
			Particle *ParticleA = EvGen->GetDecayProduct(0, 1);
			Particle *ParticleB = EvGen->GetDecayProduct(1, 1);

			if (ParticleA != 0 && ParticleA->Pdg() == 111)
				EvGen->Pi0Decay(ParticleA, ParticleA, ParticleB);

			if (ParticleA != 0 && ParticleB != 0)
			{
				EnergyA = ParticleA->E();
				MomentA = ParticleA->P();
				TransvA = ParticleA->Pt();
				ThetaA = ParticleA->Theta();
				PhiA = ParticleA->Phi();
				MassA = ParticleA->M();
				InA = ParticleA->TrackIn();
				OutA = ParticleA->TrackOut();
			
				EnergyB = ParticleB->E();
				MomentB = ParticleB->P();
				TransvB = ParticleB->Pt();
				ThetaB = ParticleB->Theta();
				PhiB = ParticleB->Phi();
				MassB = ParticleB->M();
				InB = ParticleB->TrackIn();
				OutB = ParticleB->TrackOut();
		
				Angle = ParticleA->Direction().Angle(ParticleB->Direction());

				TLorentzVector Reco = ParticleA->GetP4() + ParticleB->GetP4();
			
				Energy0 = Reco.E();
			        Moment0 = Reco.P();
			        Transv0 = Reco.Pt();
			        Theta0 = Reco.Theta();
			        Phi0 = Reco.Phi();
			        Mass0 = Reco.M();

				//std::cout << "Mass0 " << Mass0 << std::endl;
				//std::cout << "E02 " << Energy0 << std::endl;
				//std::cout << "P02 " << Moment0 << std::endl;
		
				Data->Fill();
				++ID;
			}
			else ++ND;

			delete ParticleA, ParticleB;
			ParticleA = 0;
			ParticleB = 0;

			Data->Fill();
			if (SaveMe++ > 10000)
			{
				OutFile->cd();
				Data->Write();
				SaveMe = 0;
			}
		}
			*/
	}

	OutFile->cd();
        Data->Write();

	std::cout << "There are " << ND << " (" << double(ND)/(ND+ID) * 100;
	std::cout << "\%) particles under detection threshold." << std::endl;

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
	std::cout <<"\n  -n,  --number" << std::endl;
	std::cout << "\t\tNumber of events" << std::endl;
	std::cout <<"\n  -o,  --output" << std::endl;
	std::cout << "\t\tOutput file" << std::endl;
	std::cout <<"\n  -h,  --help" << std::endl;
	std::cout << "\t\tPrint this message and exit" << std::endl;
}
