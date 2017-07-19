#include "Background.h"

//Background::Background(std::string EventDB, std::string DetectorConfig, std::string RootFile, std::string Channel)	: //Decay rates calculator
Background::Background(std::string EventDB, std::string DetectorConfig, std::string Channel)	: //Decay rates calculator
	M_Neutrino(0.0),
	M_Photon(0.0),
	M_Electron(Const::fMElectron),
	M_Muon(Const::fMMuon),
	M_Pion(Const::fMPion),
	M_Pion0(Const::fMPion0),
	M_Kaon(Const::fMKaon),
	M_Kaon0(Const::fMKaon0),
	Global(0)
{
 	InFile = new TFile(EventDB.c_str(), "READ");	//Don't close it until the end
	InFile->cd();
 	Genie = dynamic_cast<TTree*> (InFile->Get("gtree"));
	genie::NtpMCTreeHeader *thdr = dynamic_cast<genie::NtpMCTreeHeader*> (InFile->Get("header"));
 	if(!Genie)
		std::cout << "No tree found in genie file" << std::endl;
 	NEvt = Genie->GetEntries();
	gEvRec = 0;
 	Genie->SetBranchAddress("gmcrec", &gEvRec);	//fetch event

	TheBox = new Detector(DetectorConfig);

	GenMT = new TRandom3;

	TheChan.assign(Channel);

	//OutFile = new TFile(RootFile.c_str(), "RECREATE");	//Output file, keep open
	//OutFile->cd();

//	InitTree();
	InitMap();
}

Background::~Background()
{
	/*
	OutFile->Close();
	InFile->Close();

	for (iP = vParticle.begin() ; iP != vParticle.end(); ++iP)
		delete (*iP);
	vParticle.clear();	//reset particle array at each loop
	pCount.clear();

	delete TheBox;
	delete GenMT;
	*/
}

void Background::InitTree()
{
	Data = new TTree("Data", "All data from GENIE");	//Tree with candidate events

	Data->Branch("ID", &ID, "iID/i");	//global event ID
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
}

//Initialisation of map
void Background::InitMap()
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
	return Data;
}

std::string Background::GetChannel()
{
	return TheChan;
}

//Load tree
void Background::LoadTree()
{
	if (ParticleA == 0 || ParticleB == 0)
		std::cout << "You lost a particle!" << std::endl;
	else 
	{
		Hadron = pCount["Hadron"];

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

		Data->Fill();
	}
}


void Background::Loop(unsigned int Save)
{
	unsigned int Span = NEvt/Save;
	for (ID = Global; ID < Global+Span && ID < NEvt; ++ID)
	{
		std::cout << "Entry " << ID << std::endl;
		Genie->GetEntry(ID);	//get event from ID
		
		//random position for event
		double PosX = GenMT->Uniform(TheBox->GetXsize());
		double PosY = GenMT->Uniform(TheBox->GetYsize());
		double PosZ = GenMT->Uniform(TheBox->GetZsize());

		for (iP = vParticle.begin() ; iP != vParticle.end(); ++iP)
			delete (*iP);
		vParticle.clear();	//reset particle array at each loop
		pCount.clear();

		genie::EventRecord & gEvent = *(gEvRec->event);		//Get event
		genie::GHepParticle * neu = gEvent.Probe();	//get probe
		genie::GHepParticle * Hep = 0;	//particle pointer for loop in event
		TIter EvIter(&gEvent);		//iterate inside the same event, skip the probe

		std::cout << "Analysing " << gEvent.GetEntries() << " particles" << std::endl;
		while((Hep = dynamic_cast<genie::GHepParticle *>(EvIter.Next())))	//loop on all particles
		{									//inside the event
			if (Hep->Status() == 1)
			{
				if (abs(Hep->Pdg()) == 111)		//special treatments for pi0
				{					//almost 100% into 2y
					Particle *PhotonA, *PhotonB;
					Pi0Decay(CreateParticle(Hep, PosX, PosY, PosZ), PhotonA, PhotonB);	//passing reference to pointers
					vParticle.push_back(PhotonA);
					vParticle.push_back(PhotonB);
				}
				else if (abs(Hep->Pdg()) < 1e9)	//no nucleus
					vParticle.push_back(CreateParticle(Hep, PosX, PosY, PosZ));		//Everything detectable
			}			//sould contains muons, electron, pions, protons, kaons and other strange and charmed kaons
		}
		
		//Run through the collected particles
		CountParticles();
		ListCount();
		Identify();

		std::cout << std::endl;
	}

	//Data->Write();
	Global = ID;	//ID should be +1 
}

Particle* Background::CreateParticle(genie::GHepParticle *Hep, double PosX, double PosY, double PosZ)
{
	TVector3 Pos(PosX, PosY, PosZ);
	TLorentzVector P4(*Hep->P4());
	Particle *P = new Particle(Hep->Pdg(), Hep->Charge(), P4, Pos);

	if (P->Charge() != 0)
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
	if (pCount.size() == 0)
		std::cout << "No interesting particles " << std::endl;
	else
	{
		std::map<std::string, int>::iterator im = pCount.begin();
		for ( ; im != pCount.end(); ++im)
		{
			std::cout << im->first << ": " << im->second << "\t";
		}
		std::cout << std::endl;
	}
}


/***** Identification of single particles and counting *****/
bool Background::CountParticles()	//need to add smearing
{
	ParticleA = 0, ParticleB = 0;

	int i = 0;
	for (iP = vParticle.begin(); iP != vParticle.end(); ++iP, ++i)
	{
		if (TheBox->IsDetectable(*iP))
		{
			if ((*iP)->Pdg() == 13)		//Muon
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
					switch (mapChan[GetChannel()])
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
			else if ((*iP)->Pdg() == 211)	//Pion+
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
					switch (mapChan[GetChannel()])
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
			else if ((*iP)->Pdg() == 11)	//Electron
			{
				Count("Electron");
				switch (mapChan[GetChannel()])
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
						ParticleA = *iP;
						break;
					default:
						break;
				}
			}
			else if ((*iP)->Pdg() == 22)	//Photons
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
					switch (mapChan[GetChannel()])
					{
						case _nPI0:
							!ParticleA ? ParticleA = *iP : ParticleB = *iP;	//Fill first PA, then PB
							break;
						default:
							break;
					}
				}
			}
			else IdentifyHadron(*iP);
		}
	}
}


/***** Purification of events, will also add background contribution ********/
//if event matach cut, then tree is filled
bool Background::Identify()
{
	switch(mapChan[GetChannel()])
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
			if (IdentifynEE())	//A = e, B = e
				LoadTree();
			break;
		case _nEMU:
			if (IdentifynEMU())	//A = e, B = mu
				LoadTree();
			break;
		case _nMUE:
			if (IdentifynEMU())	//A = e, B = mu
				LoadTree();
			break;
		case _nPI0:
			if (IdentifynPI0())	//A = y, B = y
				LoadTree();
			break;
		case _EPI:
			if (IdentifyEPI())	//A = e, B = pi
				LoadTree();
			break;
		case _nMUMU:
			if (IdentifynMUMU())	//A = mu, B = mu
				LoadTree();
			break;
		case _MUPI:
			if (IdentifyMUPI())	//A = mu, B = pi
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
	if (iP->Pdg() == 130 || (iP->Pdg() > 300 && iP->Pdg() < 400))
		Count("Kaon");
	else if (iP->Pdg() > 400 && iP->Pdg() < 500)
		Count("Charm");
	else if (iP->Charge() != 0 && iP->Pdg() > 1000)		//I can see charged hadrons
		Count("Hadron");					//light, strange and charmed baryons
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
void Background::Pi0Decay(Particle *Pi0, Particle *&PA, Particle *&PB)
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
