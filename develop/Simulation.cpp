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
	double mix = 1e-8;
	bool ueFlag = false, umFlag = false, utFlag = false;

	std::string Channel = "ALL";
	
	while((iarg = getopt_long(argc,argv, "d:f:m:n:c:o:x:EUTLRDMh", longopts, &index)) != -1)	
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
			case 'x':
				mix = strtod(optarg, NULL);
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
			case 'E':
				ueFlag = true;
				break;
			case 'U':
				umFlag = true;
				break;
			case 'T':
				utFlag = true;
				break;
			case 'h':
				Usage(argv[0]);
				return 1;
			default:
				break;
		}
	}

	//std::ostream &Out = (OutFile.is_open()) ? OutFile : std::cout;
	
	Neutrino *TheNu0_L, *TheNu0_R, *TheNuB_L, *TheNuB_R;
	unsigned int OptHel, OptFerm;

	if (Left)
		OptHel = Neutrino::Left;
	else if (Right)
		OptHel = Neutrino::Right;
	else
		OptHel = Neutrino::Unpolarised;

	if (Dirac)
	{
		TheNu0_L = new Neutrino(Mass, Neutrino::Left  | Neutrino::Dirac);
		TheNuB_L = new Neutrino(Mass, Neutrino::Left  | Neutrino::Dirac | Neutrino::Antiparticle);
		TheNu0_R = new Neutrino(Mass, Neutrino::Right | Neutrino::Dirac);
		TheNuB_R = new Neutrino(Mass, Neutrino::Right | Neutrino::Dirac | Neutrino::Antiparticle);
		OutName += "_dirac_";
	}
	else
	{
		TheNu0_L = TheNuB_L = new Neutrino(Mass, Neutrino::Left  | Neutrino::Majorana);
		TheNu0_R = TheNuB_R = new Neutrino(Mass, Neutrino::Right | Neutrino::Majorana);
		OutName += "_major_";
	}

	mix = sqrt(mix);
	if (ueFlag)
	{
		TheNu0_L->SetMixings(mix, mix/10., mix/10.);
		TheNuB_L->SetMixings(mix, mix/10., mix/10.);
		TheNu0_R->SetMixings(mix, mix/10., mix/10.);
		TheNuB_R->SetMixings(mix, mix/10., mix/10.);
	}
	else if (umFlag)
	{
		TheNu0_L->SetMixings(mix/10., mix, mix/10.);
		TheNuB_L->SetMixings(mix/10., mix, mix/10.);
		TheNu0_R->SetMixings(mix/10., mix, mix/10.);
		TheNuB_R->SetMixings(mix/10., mix, mix/10.);
	}
	else if (utFlag)
	{
		TheNu0_L->SetMixings(mix/10., mix/10., mix);
		TheNuB_L->SetMixings(mix/10., mix/10., mix);
		TheNu0_R->SetMixings(mix/10., mix/10., mix);
		TheNuB_R->SetMixings(mix/10., mix/10., mix);
	}
	else
	{
		TheNu0_L->SetMixings(mix, mix, mix);
		TheNuB_L->SetMixings(mix, mix, mix);
		TheNu0_R->SetMixings(mix, mix, mix);
		TheNuB_R->SetMixings(mix, mix, mix);
	}

	TheNu0_L->SetDecayChannel(Channel);
	TheNuB_L->SetDecayChannel(Channel);
	TheNu0_R->SetDecayChannel(Channel);
	TheNuB_R->SetDecayChannel(Channel);

	if (!TheNu0_L->IsDecayAllowed() || !TheNu0_L->IsProductionAllowed() ||
	    !TheNuB_L->IsDecayAllowed() || !TheNuB_L->IsProductionAllowed() ||
	    !TheNu0_R->IsDecayAllowed() || !TheNu0_R->IsProductionAllowed() ||
	    !TheNuB_R->IsDecayAllowed() || !TheNuB_R->IsProductionAllowed() )
	{
		std::cerr << "Decay " << Channel << " is not allowed for a mass of " << Mass << " GeV" << std::endl;
		return 1;
	}

	
	std::string FileName;
	std::stringstream ssL;
	ssL << OutName << std::setfill('0') << std::setw(4) << Mass*1000 << ".root";
	std::cout << "Saving in " << ssL.str() << std::endl;
	TFile *OutFile = new TFile(ssL.str().c_str(), "RECREATE");

	Tracker *TheBox = new Tracker(DetConfig);
	Engine *TheEngine = new Engine(FluxConfig, 2, 2);	//creating 1FHC and 1RHC fluxedrivers

	//Binding neutrino to driver
	TheEngine->BindNeutrino(TheNu0_L, Engine::FHC, 0);
	TheEngine->BindNeutrino(TheNu0_R, Engine::FHC, 1);
	TheEngine->BindNeutrino(TheNuB_L, Engine::RHC, 0);
	TheEngine->BindNeutrino(TheNuB_R, Engine::RHC, 1);

	TheEngine->MakeFlux();
	TheEngine->ScaleDetector(TheBox);

	double S, E;
	std::vector<double> vWeight;
	double Total = TheEngine->MakeSampler(TheBox, vWeight);
	std::cout << "total number of events for " << OutFile->GetName() << "is " << Total << std::endl;
	std::vector<double> vRatio = vWeight;
	for (unsigned int i = 0; i < vWeight.size(); ++i)
	{
		vRatio.at(i) = vWeight.at(i) / Total;
		vWeight.at(i) /= Nevent * TheEngine->RangeWidth(S, E);
	}
	//if (!Dirac)
	//	for (unsigned int i = vWeight.size()/2; i > 0; --i)
	//	{
	//		vWeight.at(i-1) += vWeight.at(i-1 + vWeight.size()/2);
	//		vWeight.pop_back();
	//	}

	double True, W, R;
	bool H, P;
	double EnergyA, EnergyB, Energy0;
	double MomentA, MomentB, Moment0;
	double TransvA, TransvB, Transv0;
	double ThetaA, ThetaB, Theta0;
	double PhiA, PhiB, Phi0;
	double MassA, MassB, Mass0;
	double InnA, OutA, TotA;
	double InnB, OutB, TotB;
	double Angle;				//separation between A & B

	TRandom3 *Ran = new TRandom3(0);
	TTree *Data = new TTree("Data", "Particle tracks");

	Data->Branch("True",   &True,   "fTrue/D");
	Data->Branch("W", &W, "fW/D");
	Data->Branch("R", &R, "fR/D");
	Data->Branch("H", &H, "bH/O");
	Data->Branch("P", &P, "bP/O");

	Data->Branch("E_A",   &EnergyA, "fEnergyA/D");
	Data->Branch("P_A",   &MomentA, "fMomentA/D");
	Data->Branch("T_A",   &TransvA, "fTransvA/D");
	Data->Branch("TheA",  &ThetaA,  "fThetaA/D" );
	Data->Branch("PhiA",  &PhiA,    "fPhiA/D"   );
	Data->Branch("M_A",   &MassA,   "fMassA/D"  );
	Data->Branch("Inn_A", &InnA,    "fInnA/D"   );
	Data->Branch("Out_A", &OutA,    "fOutA/D"   );
	Data->Branch("Tot_A", &TotA,    "fTotA/D"   );

	Data->Branch("E_B",   &EnergyB, "fEnergyB/D");
	Data->Branch("P_B",   &MomentB, "fMomentB/D");
	Data->Branch("T_B",   &TransvB, "fTransvB/D");
	Data->Branch("TheB",  &ThetaB,  "fThetaB/D" );
	Data->Branch("PhiB",  &PhiB,    "fPhiB/D"   );
	Data->Branch("M_B",   &MassB,   "fMassB/D"  );
	Data->Branch("Inn_B", &InnB,    "fInnB/D"   );
	Data->Branch("Out_B", &OutB,    "fOutB/D"   );
	Data->Branch("Tot_B", &TotB,    "fTotB/D"   );

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
		std::vector<double> vEnergy, vIntens;
		TheEngine->SampleEnergy(vEnergy, vIntens);

		Neutrino *TheNu;
		for (unsigned int i = 0; i < vEnergy.size(); ++i)
		{
			if (vEnergy.at(i) < 0)
				continue;

			if (Dirac)
			{
				True = vEnergy.at(i);
				W = vWeight.at(i);
				R = vRatio.at(i);
				switch (i)
				{
					case 0:
						TheNu = TheNu0_L;
						H = false;
						P = true;
						break;
					case 1:
						TheNu = TheNu0_R;
						H = true;
						P = true;
						break;
					case 2:
						TheNu = TheNuB_L;
						H = false;
						P = false;
						break;
					case 3:
						TheNu = TheNuB_R;
						H = true;
						P = false;
						break;
				}
			}
			else if (i < 2)
			{
				W = (vWeight.at(i) + vWeight.at(i+2))/2.0;
				R = vRatio.at(i) + vRatio.at(i+2);
				double i0 = vIntens.at(i);
				double iB = vIntens.at(i+2); 
				switch (i)
				{
					case 0:
						if (Ran->Rndm() < i0/(iB+i0))
						{
							TheNu = TheNu0_L;
							True = vEnergy.at(i);
						}
						else
						{
							TheNu = TheNuB_L;
							True = vEnergy.at(i+2);
						}
						H = false;
						break;
					case 1:
						if (Ran->Rndm() < i0/(iB+i0))
						{
							TheNu = TheNu0_R;
							True = vEnergy.at(i);
						}
						else
						{
							TheNu = TheNuB_R;
							True = vEnergy.at(i);
						}
						H = true;
						break;
				}
			}
			else
				continue;

			TheNu->SetEnergy(vEnergy.at(i));
			std::vector<Particle*> vParticle = TheNu->DecayPS();

			Particle *ParticleA = 0, *ParticleB = 0;

			std::vector<Particle*>::iterator iP;
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

			TheBox->TrackReconstruct(ParticleA);
			TheBox->TrackReconstruct(ParticleB);

			if (ParticleA && ParticleB)
			{
				EnergyA = ParticleA->Energy();
				MomentA = ParticleA->Momentum();
				TransvA = ParticleA->Transverse();
				ThetaA  = ParticleA->Theta();
				PhiA    = ParticleA->Phi();
				MassA   = ParticleA->Mass();
				InnA    = ParticleA->TrackIn();
				OutA    = ParticleA->TrackOut();
				TotA    = ParticleA->TrackTot();
			
				EnergyB = ParticleB->Energy();
				MomentB = ParticleB->Momentum();
				TransvB = ParticleB->Transverse();
				ThetaB  = ParticleB->Theta();
				PhiB    = ParticleB->Phi();
				MassB   = ParticleB->Mass();
				InnB    = ParticleB->TrackIn();
				OutB    = ParticleB->TrackOut();
				TotB    = ParticleB->TrackTot();
		
				Angle   = ParticleA->Direction().Angle(ParticleB->Direction());

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

			if (ID % 10000 == 0)
			{
				OutFile->cd();
				Data->Write();
			}

			//delete ParticleA, ParticleB;
			ParticleA = 0, ParticleB = 0;
			for (unsigned k = 0; k < vParticle.size(); ++k)
				delete vParticle.at(k);
			vParticle.clear();
		}
	}

	unsigned int Size = vWeight.size() / (Dirac ? 1 : 2);
	std::cout << "Above detection thresholds there are " << (ID * 100.0 / Size)/ND << "\% of simulated particles" << std::endl;
	OutFile->cd();
        Data->Write();
	OutFile->Close();

	//delete Ran;
	//Ran = 0;
	//delete TheNu0_L;
	//TheNu0_L = 0;
	//delete TheNuB_L;
	//TheNuB_L = 0;
	//delete TheNu0_R;
	//TheNu0_R = 0;
	//delete TheNuB_R;
	//TheNuB_R = 0;

	//delete TheBox;
	//TheBox = 0;
	//delete TheEngine;
	//TheEngine = 0;

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
