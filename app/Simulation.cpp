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
#include "Physics.h"
#include "Detector.h"
#include "Flux.h"

void Usage(char* argv0);

int main(int argc, char** argv)
{
	const struct option longopts[] = 
	{
		{"detconfig", 	required_argument,	0, 'd'},
		{"fluxconfig", 	required_argument,	0, 'f'},
		{"channel", 	required_argument,	0, 'c'},
		{"mass", 	required_argument,	0, 'm'},
		{"number", 	required_argument,	0, 'n'},
		{"output", 	required_argument,	0, 'o'},
		{"left", 	no_argument,		0, 'L'},
		{"right", 	no_argument,		0, 'R'},
		{"help", 	no_argument,		0, 'h'},
		{0,	0, 	0,	0},
	};

	int index; 
	int iarg = 0;
	opterr = 1;
	
	std::string SMConfig, DetConfig, FluxConfig;
	TFile *OutFile;
	unsigned int Nevent = 1000;
	double Mass = 0.0;
	bool Left = false, Right = false;	//default unpolarised

	std::string Channel = "ALL";
	
	while((iarg = getopt_long(argc,argv, "d:f:m:n:c:o:LRh", longopts, &index)) != -1)	
	{
		switch(iarg)
		{
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
				break;
			case 'L':
				Left = true;
				Right = false;
				break;
			case 'R':
				Left = false;
				Right = true;
				break;
			case 'h':
				Usage(argv[0]);
				return 1;
			default:
				break;
		}
	}

	//std::ostream &Out = (OutFile.is_open()) ? OutFile : std::cout;
	
	Neutrino *TheNu;
	if (Left)
		TheNu = new Neutrino(0, Neutrino::Dirac | Neutrino::Left );
	else if (Right)
		TheNu = new Neutrino(0, Neutrino::Dirac | Neutrino::Right );
	else
		TheNu = new Neutrino(0, Neutrino::Dirac | Neutrino::Unpolarised );

	Tracker *TheBox = new Tracker(DetConfig);
	Engine *TheEngine = new Engine(FluxConfig, 1, 1);	//creating 1FHC and 1RHC fluxedrivers

	//Binding neutrino to driver
	TheEngine->BindNeutrino(TheNu, Engine::FHC, 0);
	TheEngine->BindNeutrino(TheNu, Engine::RHC, 0);

	if (Mass)
		TheNu->SetMass(Mass);

	TheNu->SetMixings(1e-5/sqrt(3.0), 1e-5/sqrt(3.0), 1e-5/sqrt(3.0));
	TheNu->SetDecayChannel(Channel);

	TheEngine->MakeFlux();
	TheEngine->MakeSampler(TheBox);

	Particle *ParticleA, *ParticleB;
	int iA, iB;
	unsigned int SaveMe = 0;

	double Real, Delay;

	double EnergyA, EnergyB, EnergyC, Energy0;
	double MomentA, MomentB, MomentC, Moment0;
	double TransvA, TransvB, TransvC, Transv0;
	double ThetaA, ThetaB, ThetaC, Theta0;
	double PhiA, PhiB, PhiC, Phi0;
	double MassA, MassB, MassC, Mass0;
	double InnA, BakA, OutA;
	double InnB, BakB, OutB;
	double Angle;				//separation between A & B

	//TRandom3 *Rand = new TRandom3(0);
	TTree *Data = new TTree("Data", "Particle tracks");

	Data->Branch("Real",  &Real, "  fReal/D"    );
	//Data->Branch("Delay", &Delay, "fDelay/D");

	Data->Branch("E_A",   &EnergyA, "fEnergyA/D");
	Data->Branch("P_A",   &MomentA, "fMomentA/D");
	Data->Branch("T_A",   &TransvA, "fTransvA/D");
	Data->Branch("TheA",  &ThetaA,  "fThetaA/D" );
	Data->Branch("PhiA",  &PhiA,    "fPhiA/D"   );
	Data->Branch("M_A",   &MassA,   "fMassA/D"  );
	Data->Branch("Inn_A", &InnA,    "fInnA/D"   );
	Data->Branch("Bak_A", &BakA,    "fBakA/D"   );
	Data->Branch("Out_A", &OutA,    "fOutA/D"   );

	Data->Branch("E_B",   &EnergyB, "fEnergyB/D");
	Data->Branch("P_B",   &MomentB, "fMomentB/D");
	Data->Branch("T_B",   &TransvB, "fTransvB/D");
	Data->Branch("TheB",  &ThetaB,  "fThetaB/D" );
	Data->Branch("PhiB",  &PhiB,    "fPhiB/D"   );
	Data->Branch("M_B",   &MassB,   "fMassB/D"  );
	Data->Branch("Inn_B", &InnB,    "fInnB/D"   );
	Data->Branch("Bak_B", &BakB,    "fBakB/D"   );
	Data->Branch("Out_B", &OutB,    "fOutB/D"   );

	Data->Branch("Angle", &Angle,   "fAngle/D"  );

	Data->Branch("E_0",   &Energy0, "fEnergy0/D");
	Data->Branch("P_0",   &Moment0, "fMoment0/D");
	Data->Branch("T_0",   &Transv0, "fTransv0/D");
	Data->Branch("The0",  &Theta0,  "fTheta0/D" );
	Data->Branch("Phi0",  &Phi0,    "fPhi0/D"   );
	Data->Branch("M_0",   &Mass0,   "fMass0/D"  );

	//unsigned int nBunch = 80;	//80 bunhces 
	//double tBunch = 20e-9;	//20ns between bunches
	//unsigned int iBunch = 0;
	//double Beta = 1.0;
	//double Lc = EvGen->GetDetectorPtr()->GetElement("Baseline")/Const::fC;
	
	unsigned int ND = 0, ID;
	for (ID = 0; ID < Nevent; ++ND)
	{
		std::cout << "Sampling " << ND << std::endl;
		//RealFHC = TheEngine->SampleEnergy(Engine::FHC, 0, 1);
		Real = TheEngine->SampleEnergy(Engine::FHC, 0, 1);

		/*
		iBunch = Rand->Integer(nBunch);
		if (Rand->Rndm() < 2.99)
			Beta = sqrt(1 - Mass*Mass / Real/Real);
		else
			Beta = 1.0;

		//Delay = Lc * (1.0/Beta - 1) + iBunch*tBunch;
		Delay = Rand->Gaus(Lc * (1.0/Beta - 1) + iBunch*tBunch, 1e-9);


		*/


		std::vector<Particle*> vParticle = TheNu->DecayPS();		//it should not be empty

		ParticleA = 0, ParticleB = 0;
		iA = -1, iB = -1;
		std::vector<Particle*>::iterator iP;
		//for (unsigned int i = 0; i < vParticle.size(); ++i)
		for (iP = vParticle.begin(); iP != vParticle.end(); ++iP)
		{
			if ((*iP)->Pdg() == 12)	//neutrino is invibisle
				continue;
			else if ((*iP)->Pdg() == 111)		//pi0, must decay rn
				TheBox->Pi0Decay(*iP, ParticleA, ParticleB);
			else if ((*iP)->Pdg() == 22)		//nu gamma decay
				ParticleA = ParticleB = *iP;
			else
			{
				//if (iA < 0)
				if (!ParticleA)
					ParticleA = *iP;
				else if (!ParticleB)
					ParticleB = *iP;
				else
				{
					delete *iP;
					*iP = 0;
				}
			}
		}

		TheBox->TrackReconstruct(ParticleA);
		TheBox->TrackReconstruct(ParticleB);
		if (ParticleA != 0 && ParticleB != 0)
		{
			std::cout << "\tParticle " << ID << std::endl;
			EnergyA = ParticleA->Energy();
			MomentA = ParticleA->Momentum();
			TransvA = ParticleA->Transverse();
			ThetaA  = ParticleA->Theta();
			PhiA    = ParticleA->Phi();
			MassA   = ParticleA->Mass();
			InnA    = ParticleA->TrackIn();
			BakA    = ParticleA->TrackBack();
			OutA    = ParticleA->TrackOut();
		
			EnergyB = ParticleB->Energy();
			MomentB = ParticleB->Momentum();
			TransvB = ParticleB->Transverse();
			ThetaB  = ParticleB->Theta();
			PhiB    = ParticleB->Phi();
			MassB   = ParticleB->Mass();
			InnB    = ParticleB->TrackIn();
			BakB    = ParticleB->TrackBack();
			OutB    = ParticleB->TrackOut();
	
			Angle = ParticleA->Direction().Angle(ParticleB->Direction());

			TLorentzVector Reco = ParticleA->FourVector() + ParticleB->FourVector();
		
			Energy0 = Reco.E();
			Moment0 = Reco.P();
			Transv0 = Reco.Pt();
			Theta0  = Reco.Theta();
			Phi0    = Reco.Phi();
			Mass0   = Reco.M();

			Data->Fill();
			++ID;
		}

		if (SaveMe++ > 10000)
		{
			OutFile->cd();
			Data->Write();
			SaveMe = 0;
		}

		delete ParticleA, ParticleB;
		ParticleA = 0, ParticleB = 0;
	}

	std::cout << "Above detection thresholds there are " << (100.0*ID)/ND << "\% of simulated particles" << std::endl;
	OutFile->cd();
        Data->Write();
	OutFile->Close();

	return 0;
}

void Usage(char* argv0)
{
	std::cout << "Description" << std::endl;
	std::cout << "Usage : " << std::endl;
	std::cout << argv0 << " [OPTIONS]" << std::endl;
	std::cout <<"\n  -d,  --detconfig" << std::endl;
	std::cout << "\t\tDetector configuration file" << std::endl;
	std::cout <<"\n  -f,  --fluxconfig" << std::endl;
	std::cout << "\t\tFlux configuration file" << std::endl;
	std::cout <<"\n  -m,  --mass" << std::endl;
	std::cout << "\t\tNeutrino mass" << std::endl;
	std::cout <<"\n  -c,  --channel" << std::endl;
	std::cout << "\t\tDecay channel, defaul ALL" << std::endl;
	std::cout <<"\n  -n,  --number" << std::endl;
	std::cout << "\t\tNumber of events" << std::endl;
	std::cout <<"\n  -o,  --output" << std::endl;
	std::cout << "\t\tOutput root file" << std::endl;
	std::cout <<"\n  -L,  --left" << std::endl;
	std::cout << "\t\tNeutrino polarisation left, default unpolarised" << std::endl;
	std::cout <<"\n  -R,  --right" << std::endl;
	std::cout << "\t\tNeutrino polarisation right, default unpolarised" << std::endl;
	std::cout <<"\n  -h,  --help" << std::endl;
	std::cout << "\t\tPrint this message and exit" << std::endl;
}
