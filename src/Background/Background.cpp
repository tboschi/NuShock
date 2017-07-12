#include "Background.h"

Background::Background(std::string EventDB, std::string DetectorConfig)	: //Decay rates calculator
	M_Neutrino(0.0),
	M_Photon(0.0),
	M_Electron(Const::fMElectron),
	M_Muon(Const::fMMuon),
	M_Pion(Const::fMPion),
	M_Pion0(Const::fMPion0),
	M_Kaon(Const::fMKaon),
	M_Kaon0(Const::fMKaon0)
{
	std::cout << "H0" << std::endl;
	//GetCommandLineArgs (argc, argv);	//Useful?
	std::cout << "H1" << std::endl;
	//genie::NtpMCTreeHeader *Gthdr = 0;
 	InFile = new TFile(EventDB.c_str(), "READ");
	InFile->cd();
	std::cout << "H2" << std::endl;
	genie::NtpMCTreeHeader *thdr = 0;
 	Genie = dynamic_cast<TTree*> (InFile->Get("gtree"));
	thdr = dynamic_cast<genie::NtpMCTreeHeader *> (InFile->Get("header"));
	std::cout << "H3" << std::endl;
 	if(!Genie)
		std::cout << "No tree found in genie file" << std::endl;
	std::cout << Genie->GetEntries() << std::endl;
	std::cout << "H4" << std::endl;

	std::cout << "H5" << std::endl;
	//gEvRec = 0;
 	//Genie->SetBranchAddress("gmcrec", &gEvRec);
	//std::cout << "mcrec " << gEvRec << std::endl;
 	//Genie->SetBranchStatus("*", 0);
 	//Genie->SetBranchStatus("gmcrec", 1);
	std::cout << "H6" << std::endl;
 	// Get the nbr of evts to analyse (-n argument)
 	//NEV = (gOptNEvt > 0) ? TMath::Min(gOptNEvt, Gtree->GetEntries()) : (int)Gtree->GetEntries();
	std::cout << "H7" << std::endl;
	std::cout << Genie << std::endl;
 	NEvt = Genie->GetEntries();
	std::cout << "H8" << NEvt<< std::endl;

	TheBox = new Detector(DetectorConfig);
	std::cout << "H9" << std::endl;

	GenMT = new TRandom3;
	std::cout << "H10" << std::endl;

	/*
	Data = new TTree("Data", "All data from GENIE");
	std::cout << "H11" << std::endl;

	Data->Branch("ID", &ID, "iID/i");	//global event ID
	std::cout << "H12" << std::endl;
	Data->Branch("NH", &Hadron, "iHadron/i");	//global event ID

	//Particle A
	Data->Branch("E_A", &EnergyA, "fEnergyA/D");
	Data->Branch("P_A", &MomentA, "fMomentA/D");
	Data->Branch("T_A", &TransvA, "fTransvA/D");
	Data->Branch("TheA", &ThetaA, "fThetaA/D");
	Data->Branch("PhiA", &PhiA, "fPhiA/D");
	Data->Branch("M_A", &MassA, "fMassA/D");

	//Particle B
	Data->Branch("E_B", &EnergyB, "fEnergyB/D");
	Data->Branch("P_B", &MomentB, "fMomentB/D");
	Data->Branch("T_B", &TransvB, "fTransvB/D");
	Data->Branch("TheB", &ThetaB, "fThetaB/D");
	Data->Branch("PhiB", &PhiB, "fPhiB/D");
	Data->Branch("M_B", &MassB, "fMassB/D");

	//Particle 0 = A + B
	Data->Branch("E_0", &Energy0, "fEnergy0/D");
	Data->Branch("P_0", &Moment0, "fMoment0/D");
	Data->Branch("T_0", &Transv0, "fTransv0/D");
	Data->Branch("The0", &Theta0, "fTh0ta0/D");
	Data->Branch("Phi0", &Phi0, "fPhi0/D");
	Data->Branch("M_0", &Mass0, "fMass0/D");
	std::cout << "H13" << std::endl;
	*/
}

//Initialisation of map
void Background::MapInit()
{
	mapChan["ALL"] = _ALL;
	mapChan["nnn"] = _nnn;
	mapChan["nGAMMA"] = _nGAMMA;
	mapChan["nEE"] = _nEE;
	mapChan["nEMU"] = _nEMU;
	mapChan["nMUE"] = _nMUE;
	mapChan["nPI0"] = _nPI0;
	mapChan["EPI"] = _EPI;
	mapChan["nMUMU"] = _nMUMU;
	mapChan["MUPI"] = _MUPI;
	mapChan["EKA"] = _EKA;
	mapChan["nKA0"] = _nKA0;
}

TTree *Background::GetTree()
{
	//return Data;
}

//Load tree
void Background::LoadTree()
{
	if (ParticleA == 0 || ParticleB == 0)
		std::cout << "You lost a particle!" << std::endl;
	else 
	{
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

		TLorentzVector Reco = ParticleA->GetP4() + ParticleB->GetP4();
	
		Energy0 = Reco.E();
	        Moment0 = Reco.P();
	        Transv0 = Reco.Pt();
	        Theta0 = Reco.Theta();
	        Phi0 = Reco.Phi();
	        Mass0 = Reco.M();

		//Data->Fill();
	}
}


void Background::Loop(std::string Channel)
{
	std::cout << Genie << std::endl;
	std::cout << Genie->GetEntries() << std::endl;
	std::cout << "o0" << std::endl;
	genie::NtpMCEventRecord * gEvRec = 0;
	std::cout << Genie->GetEntries() << std::endl;
	std::cout << "o1" << std::endl;
 	Genie->SetBranchAddress("gmcrec", &gEvRec);
	std::cout << "o2" << std::endl;
	std::cout << "NEvt " << NEvt << std::endl;
	std::cout << "o3" << std::endl;
	for (ID = 1; ID < NEvt; ++ID)
	{
		std::cout << "L0" << std::endl;
		std::cout << "GetEntry " << Genie->GetEntry(ID) << std::endl;;	//get event from ID
		std::cout << "L1" << std::endl;

		/*
		//random position for event
		double PosX = GenMT->Uniform(TheBox->GetXsize());
		double PosY = GenMT->Uniform(TheBox->GetYsize());
		double PosZ = GenMT->Uniform(TheBox->GetZsize());
		std::cout << "L2" << std::endl;

		vParticle.clear();	//reset particle array at each loop
		std::cout << "L3" << std::endl;

		genie::EventRecord & gEvent = *(gEvRec->event);
		std::cout << "L4" << std::endl;
		genie::GHepParticle * neu = gEvent.Probe();	//get probe
		std::cout << "L5" << std::endl;
		genie::GHepParticle * Hep = 0;	//particle pointer for loop in event
		std::cout << "L6" << std::endl;
		TIter EvIter(&gEvent);		//iterate inside the same event, skip the probe
		std::cout << "L7" << std::endl;

		while((Hep = dynamic_cast<genie::GHepParticle *>(EvIter.Next())))	//loop on all particles
		{									//inside the event
			std::cout << "W0" << std::endl;
			if (Hep->Status() == 1)
			{
			std::cout << "W1" << std::endl;
				if (abs(Hep->Pdg()) == 111)		//special treatments for pi0
				{					//almost 100% into 2y
			std::cout << "W2" << std::endl;
					Particle *PhotonA, *PhotonB;
					Pi0Decay(CreateParticle(Hep, PosX, PosY, PosZ), PhotonA, PhotonB);	
					vParticle.push_back(PhotonA);
					vParticle.push_back(PhotonB);
				}
				else if (abs(Hep->Pdg()) < 1e9)	//no nucleus
					vParticle.push_back(CreateParticle(Hep, PosX, PosY, PosZ));		//Everything detectable
			std::cout << "W3" << std::endl;
			}			//sould contains muons, electron, pions, protons, kaons and other strange and charmed kaons
		}
		*/
		//Run through the collected particles
		std::cout << "L8" << std::endl;
		//CountParticles(Channel);
		std::cout << "L9" << std::endl;
		//Identify(Channel);
		std::cout << "L10" << std::endl;
	}
}

Particle* Background::CreateParticle(genie::GHepParticle *Hep, double PosX, double PosY, double PosZ)
{
	TVector3 Pos(PosX, PosY, PosZ);
	TLorentzVector P4(*Hep->P4());
	Particle *P = new Particle(Hep->Pdg(), Hep->Charge(), P4, Pos);

	TheBox->SignalSmearing(GenMT, P);

	return P;
}


int Background::Count(std::string PartName, int N)
{
	if (N < 0)
		pCount[PartName] = 0;
	else pCount[PartName] += N;

	return pCount[PartName];
}

void Background::ListCount()
{
	std::map<std::string, int>::iterator im = pCount.begin();
	for ( ; im != pCount.end(); ++im)
	{
		std::cout << im->first << ": " << im->second << "\t";
	}
	std::cout << std::endl;
}


/***** Identification of single particles and counting *****/
bool Background::CountParticles(std::string Channel)	//need to add smearing
{
	ParticleA = 0, ParticleB = 0;

	for (iP = vParticle.begin(); iP != vParticle.end(); ++iP)
	{
		if ((*iP)->Pdg() == 13)		//Muon
		{
			if ((*iP)->Ekin() > TheBox->GetElement("ThrMuon"))	//can see it, else gnorri
			{
				if (TheBox->TrackLength(*iP) < 0.50 && GenMT->Rndm() > 0.5)	//likely a pion, 50:50 chance can be misidentified
				{
					(*iP)->SetPdg(211);
					(*iP)->SetMass(M_Pion);
					--iP;			//must recheck
				}
				else		//quite a muon!
				{
					Count("Muon");
					switch (mapChan[Channel])
					{
						case _nEMU:
							ParticleB = *iP;
							break;
						case _nMUE:
							ParticleB = *iP;
							break;
						case _nMUMU:
							!ParticleA ? ParticleA = *iP : ParticleB = *iP;	//Fill first PA, then PB
							break;
						case _MUPI:
							ParticleA = *iP;
							break;
						default:
							break;
					}
				}
			}
		}
		else if ((*iP)->Pdg() == 211)	//Pion+
		{
			if ((*iP)->Ekin() > TheBox->GetElement("ThrPion"));	//can't detect
			{
				if (TheBox->TrackLength(*iP) > 0.50 && GenMT->Rndm() > 0.5)	//likely a muon, 50:50 chance can be misidentified
				{
					(*iP)->SetPdg(13);
					(*iP)->SetMass(M_Muon);
					--iP;			//must recheck
				}
				else 		//quite a pion!
				{
					Count("Pion");
					switch (mapChan[Channel])
					{
						case _EPI:
							ParticleB = *iP;
							break;
						case _MUPI:
							ParticleB = *iP;
							break;
						default:
							break;
					}
				}
			}
		}
		else if ((*iP)->Pdg() == 11)	//Electron
		{
			if ((*iP)->Ekin() > TheBox->GetElement("ThrGamma"));	//quite an electron!
			{
				Count("Electron");
				switch (mapChan[Channel])
				{
					case _nEE:
						!ParticleA ? ParticleA = *iP : ParticleB = *iP;	//Fill first PA, then PB
						break;
					case _nEMU:
						ParticleA = *iP;
						break;
					case _nMUE:
						ParticleA = *iP;
						break;
					case _EPI:
						ParticleB = *iP;
						break;
					default:
						break;
				}
			}
		}
		else if ((*iP)->Pdg() == 22)	//Photons
		{
			if ((*iP)->Ekin() > 2*M_Electron)	//pair production can occur
			{
				if (GammaDecay() < 0.03)	//this photon could be an electron!
				{
					(*iP)->SetPdg(11);
					(*iP)->SetMass(M_Electron);
					--iP;			//must recheck
				}
				else
				{
					Count("Photon");
					switch (mapChan[Channel])
					{
						case _nPI0:
							!ParticleA ? ParticleA = *iP : ParticleB = *iP;	//Fill first PA, then PB
							break;
						default:
							break;
					}
				}
			}
		}
		else IdentifyHadron(*iP);
	}
}


/***** Purification of events, will also add background contribution ********/
//if event matach cut, then tree is filled
bool Background::Identify(std::string Channel)
{
	switch(mapChan[Channel])
	{
		/*
		case _ALL:
			Result = Total();
			break;
		case _nnn:
			Result = nnn();
			break;
		case _nGAMMA:
			Result = nGAMMA();
			break;
		*/
		case _nEE:
			if (IdentifynEE())
				LoadTree();
			break;
		case _nEMU:
			if (IdentifynEMU())
				LoadTree();
			break;
		case _nMUE:
			if (IdentifynEMU())
				LoadTree();
			break;
		case _nPI0:
			if (IdentifynPI0())
				LoadTree();
			break;
		case _EPI:
			if (IdentifyEPI())
				LoadTree();
			break;
		case _nMUMU:
			if (IdentifynMUMU())
				LoadTree();
			break;
		case _MUPI:
			if (IdentifyMUPI())
				LoadTree();
			break;
		/*
		case _EKA:
			Result = EKA();
			break;
		case _nKA0:
			Result = nKA0();
		*/
		default:
			++iP;
			break;
	}
}

/***** Not interesting hadrons, mesons are counted as kaon and charm, baryons as a whole *****/
void Background::IdentifyHadron(Particle *iP)	//here everything not selected by cut
{
	if (iP->Ekin() > TheBox->GetElement("ThrHadron"))	//hadron (proton) activity
	{
		if (iP->Pdg() == 130 || (iP->Pdg() > 300 && iP->Pdg() < 400))
			Count("Kaon");
		else if (iP->Pdg() > 400 && iP->Pdg() < 500)
			Count("Charm");
		else if (iP->Charge() != 0 && iP->Pdg() > 1000)		//I can see charged hadrons
			Count("Hadron");					//light, strange and charmed baryons
	}
}

/***** single channel check for pcount map, exact number needed *****/
bool Background::IdentifynEE()	//need to add smearing
{
	if (pCount["Electron"]	== 2	&&
	    pCount["Photon"]	== 0	&& 
	    pCount["Muon"]	== 0	&& 
	    pCount["Pion"]	== 0	&& 
	    pCount["Kaon"]	== 0	&& 
	    pCount["Charm"]	== 0	  )
		return true;
	else return false;
}

bool Background::IdentifynEMU()	//need to add smearing
{
	if (pCount["Electron"]	== 1	&&
	    pCount["Photon"]	== 0	&& 
	    pCount["Muon"]	== 1	&& 
	    pCount["Pion"]	== 0	&& 
	    pCount["Kaon"]	== 0	&& 
	    pCount["Charm"]	== 0	  )
		return true;
	else return false;
}

bool Background::IdentifynPI0()	//need to add smearing
{
	if (pCount["Electron"]	== 0	&&
	    pCount["Photon"]	== 2	&& 
	    pCount["Muon"]	== 0	&& 
	    pCount["Pion"]	== 0	&& 
	    pCount["Kaon"]	== 0	&& 
	    pCount["Charm"]	== 0	  )
		return true;
	else return false;
}

bool Background::IdentifyEPI()	//need to add smearing
{
	if (pCount["Electron"]	== 1	&&
	    pCount["Photon"]	== 0	&& 
	    pCount["Muon"]	== 0	&& 
	    pCount["Pion"]	== 1	&& 
	    pCount["Kaon"]	== 0	&& 
	    pCount["Charm"]	== 0	  )
		return true;
	else return false;
}

bool Background::IdentifynMUMU()	//need to add smearing
{
	if (pCount["Electron"]	== 0	&&
	    pCount["Photon"]	== 0	&& 
	    pCount["Muon"]	== 2	&& 
	    pCount["Pion"]	== 0	&& 
	    pCount["Kaon"]	== 0	&& 
	    pCount["Charm"]	== 0	  )
		return true;
	else return false;
}

bool Background::IdentifyMUPI()	//need to add smearing
{
	if (pCount["Electron"]	== 0	&&
	    pCount["Photon"]	== 0	&& 
	    pCount["Muon"]	== 1	&& 
	    pCount["Pion"]	== 1	&& 
	    pCount["Kaon"]	== 0	&& 
	    pCount["Charm"]	== 0	  )
		return true;
	else return false;
}
	


//Distance travelled by photons before conversion
double Background::GammaDecay()
{
	double Path = 9.0/7.0 * TheBox->InteractionLength();
	return GenMT->Exp(Path);	//in cm
}
	
//Treat pi0 decay into 2 photons
void Background::Pi0Decay(Particle *Pi0, Particle *PA, Particle *PB)
{
	//in rest frame
	TLorentzVector GammaA(0, 0, M_Pion0/2.0, M_Pion0/2.0); 
	TLorentzVector GammaB(0, 0, -M_Pion0/2.0, M_Pion0/2.0); 

	TVector3 vBoost(Pi0->GetP4().BoostVector());
	double Theta = GenMT->Uniform(-Const::fPi, Const::fPi);
	double Phi = GenMT->Uniform(-Const::fPi, Const::fPi);

	double Tau = 1.0 / (8.52e-17 * Pi0->GetP4().Beta() * Pi0->GetP4().Gamma());	//should be the relativistic tau
	double Travel = GenMT->Exp(Tau);
	TVector3 Move(Pi0->GetP4().Vect().Unit());	//direction of motion
	Move *= Travel;
	Move += Pi0->Position();

	GammaA.SetTheta(Theta);
	GammaB.SetTheta(Const::fPi + Theta);
	GammaA.SetPhi(Phi);
	GammaB.SetPhi(Phi);

	GammaA.Boost(vBoost);
	GammaB.Boost(vBoost);

	PA = new Particle(22, 0, GammaA, Move);	//here are the photons
	PB = new Particle(22, 0, GammaB, Move);	//position should be different

	delete Pi0;
}
