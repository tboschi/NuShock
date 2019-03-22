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
		{"dirac", 	no_argument,		0, 'D'},
		{"majorana", 	no_argument,		0, 'M'},
		{"help", 	no_argument,		0, 'h'},
		{0,	0, 	0,	0},
	};

	int index; 
	int iarg = 0;
	opterr = 1;
	
	std::string OutName, DetConfig, FluxConfig;
	unsigned int Nevent = 10000;
	double Mass = 0.0;
	bool Left = false, Right = false;	//default unpolarised
	bool Dirac = true;	//default unpolarised

	std::string Channel = "ALL";
	
	while((iarg = getopt_long(argc,argv, "d:f:m:n:c:o:LRDMh", longopts, &index)) != -1)	
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
				OutName.assign(optarg);
				break;
			case 'L':
				Left = true;
				Right = false;
				break;
			case 'R':
				Left = false;
				Right = true;
				break;
			case 'D':
				Dirac = true;
				break;
			case 'M':
				Dirac = false;
				break;
			case 'h':
				return 1;
			default:
				break;
		}
	}

	//std::ostream &Out = (OutFile.is_open()) ? OutFile : std::cout;
	
	TRandom3 *Ran = new TRandom3(0);
	TFile *File = new TFile(OutName.c_str(), "RECREATE");
	TTree *Data = new TTree("Data", "Particle tracks");

	double EnergyA, EnergyB, Energy0;
	double MomentA, MomentB, Moment0;
	double TransvA, TransvB, Transv0;
	double ThetaA, ThetaB, Theta0;
	double PhiA, PhiB, Phi0;
	double MassA, MassB, Mass0;
	double Angle;				//separation between A & B

	Data->Branch("E_A",   &EnergyA, "fEnergyA/D");
	Data->Branch("P_A",   &MomentA, "fMomentA/D");
	Data->Branch("T_A",   &TransvA, "fTransvA/D");
	Data->Branch("TheA",  &ThetaA,  "fThetaA/D" );
	Data->Branch("PhiA",  &PhiA,    "fPhiA/D"   );
	Data->Branch("M_A",   &MassA,   "fMassA/D"  );

	Data->Branch("E_B",   &EnergyB, "fEnergyB/D");
	Data->Branch("P_B",   &MomentB, "fMomentB/D");
	Data->Branch("T_B",   &TransvB, "fTransvB/D");
	Data->Branch("TheB",  &ThetaB,  "fThetaB/D" );
	Data->Branch("PhiB",  &PhiB,    "fPhiB/D"   );
	Data->Branch("M_B",   &MassB,   "fMassB/D"  );

	Data->Branch("Angle", &Angle,   "fAngle/D"  );

	PhaseSpace *ThePS = new PhaseSpace();

	int Hel;
	if (Left)
		Hel = -1;
	else if (Right)
		Hel = 1;
	else
		Hel = 0;

	bool Ferm;
	if (Dirac)
		Ferm = true;
	else
		Ferm = false;

	double *mix = new double[3];
	mix[0] = 1e-3;
	mix[1] = 1e-3;
	mix[2] = 1e-3;
	ThePS->SetNeutrino(Mass, mix, Ferm, 1, Hel);
	ThePS->SetRest(Mass);

	for (unsigned int loop = 0; loop < Nevent;++loop)
	{
		std::vector<Particle*> vDaughter;
		if (ThePS->Generate(Amplitude::_EPI))
			for (unsigned int i = 0; i < ThePS->Daughters(); ++i)
				vDaughter.push_back(ThePS->Daughter(i, PhaseSpace::RestFrame));


		Particle *ParticleA = 0, *ParticleB = 0;

		std::vector<Particle*>::iterator iP;
		for (iP = vDaughter.begin(); iP != vDaughter.end(); ++iP)
		{
			if ((*iP)->Pdg() == 12)	//neutrino is invibisle
				continue;
			//else if ((*iP)->Pdg() == 111)		//pi0, must decay rn
			//	TheBox->Pi0Decay(*iP, ParticleA, ParticleB);
			//else if ((*iP)->Pdg() == 22)		//nu gamma decay
			//	ParticleA = ParticleB = *iP;
			else
			{
				if (!ParticleA)
					ParticleA = *iP;
				else if (!ParticleB)
					ParticleB = *iP;
				else
				{
					//delete (*iP);
					//(*iP) = 0;
				}
			}
		}

		if (ParticleA && ParticleB)
		{
			EnergyA = ParticleA->Energy();
			MomentA = ParticleA->Momentum();
			TransvA = ParticleA->Transverse();
			ThetaA  = ParticleA->Theta();
			PhiA    = ParticleA->Phi();
			MassA   = ParticleA->Mass();
		
			EnergyB = ParticleB->Energy();
			MomentB = ParticleB->Momentum();
			TransvB = ParticleB->Transverse();
			ThetaB  = ParticleB->Theta();
			PhiB    = ParticleB->Phi();
			MassB   = ParticleB->Mass();

			Angle = ParticleA->Direction().Angle(ParticleB->Direction());
			Data->Fill();
			//++loop;
		}

		for (unsigned int i = 0; i < vDaughter.size(); ++i)
			delete vDaughter.at(i);
	}
	Data->Write();
	File->Close();

	return 0;
}
