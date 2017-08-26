#include "Background.h"

//Background::Background(std::string EventDB, std::string DetectorConfig, std::string RootFile, std::string Channel)	: //Decay rates calculator
Background::Background(std::ostream &OutStream, std::string EventDB, std::string DetectorConfig, std::string Channel)	: //Decay rates calculator
	M_Neutrino(0.0),
	M_Photon(0.0),
	M_Electron(Const::fMElectron),
	M_Muon(Const::fMMuon),
	M_Pion(Const::fMPion),
	M_Pion0(Const::fMPion0),
	M_Kaon(Const::fMKaon),
	M_Kaon0(Const::fMKaon0),
	Global(0),
	Out(OutStream)
{
 	InFile = new TFile(EventDB.c_str(), "READ");	//Don't close it until the end
	InFile->cd();
 	Genie = dynamic_cast<TTree*> (InFile->Get("gtree"));
	genie::NtpMCTreeHeader *thdr = dynamic_cast<genie::NtpMCTreeHeader*> (InFile->Get("header"));
 	if(!Genie)
		std::cerr << "No tree found in genie file" << std::endl;
 	NEvt = Genie->GetEntries();
 	//NEvt = 1e4;
	gEvRec = 0;
 	Genie->SetBranchAddress("gmcrec", &gEvRec);	//fetch event

	TheBox = new Detector(DetectorConfig);

	GenMT = new TRandom3;
	Vertex = new TVector3();

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
	Data->Branch("L_A", &LengthA, "fLengthA/D");
	Data->Branch("Lo_A", &LengthoA, "fLengthoA/D");

	//Particle B
	Data->Branch("E_B", &EnergyB, "fEnergyB/D");
	Data->Branch("P_B", &MomentB, "fMomentB/D");
	Data->Branch("T_B", &TransvB, "fTransvB/D");
	Data->Branch("TheB", &ThetaB, "fThetaB/D");
	Data->Branch("PhiB", &PhiB, "fPhiB/D");
	Data->Branch("M_B", &MassB, "fMassB/D");
	Data->Branch("L_B", &LengthB, "fLengthB/D");
	Data->Branch("Lo_B", &LengthoB, "fLengthoB/D");

	Data->Branch("Angle", &Angle, "fAngle/D");

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
		std::cerr << "You lost a particle there!" << std::endl;
	else 
	{
		Hadron = pCount["Hadron"];

		EnergyA = ParticleA->E();
		MomentA = ParticleA->P();
		TransvA = ParticleA->Pt();
		ThetaA = ParticleA->Theta();
		PhiA = ParticleA->Phi();
		MassA = ParticleA->M();
		LengthA = ParticleA->TrackIn();
		LengthoA = ParticleA->TrackOut();
	
		EnergyB = ParticleB->E();
		MomentB = ParticleB->P();
		TransvB = ParticleB->Pt();
		ThetaB = ParticleB->Theta();
		PhiB = ParticleB->Phi();
		MassB = ParticleB->M();
		LengthB = ParticleB->TrackIn();
		LengthoB = ParticleB->TrackOut();

		Angle = ParticleA->Direction().Angle(ParticleB->Direction());

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
		Out << "Entry " << ID << std::endl;
		Genie->GetEntry(ID);	//get event from ID

		for (iP = vParticle.begin() ; iP != vParticle.end(); ++iP)
			delete (*iP);
		vParticle.clear();	//reset particle array at each loop
		pCount.clear();

		genie::EventRecord & gEvent = *(gEvRec->event);		//Get event
		genie::GHepParticle * neu = gEvent.Probe();	//get probe
		genie::GHepParticle * Hep = 0;	//particle pointer for loop in event
		TIter EvIter(&gEvent);		//iterate inside the same event, skip the probe

		//random position for event
		double PosX = GenMT->Uniform(TheBox->GetXstart(), TheBox->GetXend());
		double PosY = GenMT->Uniform(TheBox->GetYstart(), TheBox->GetYend());
		double PosZ = GenMT->Uniform(TheBox->GetZstart(), TheBox->GetZend());
		Vertex->SetXYZ(PosX, PosY, PosZ);

		Out << "Event contains " << gEvent.GetEntries() << " particles" << std::endl;
		while((Hep = dynamic_cast<genie::GHepParticle *>(EvIter.Next())))	//loop on all particles
		{									//inside the event
			if (Hep->Status() == 1)
			{
				if (abs(Hep->Pdg()) == 111)		//special treatments for pi0
				{					//almost 100% into 2y
					Particle *PhotonA, *PhotonB;
					Pi0Decay(CreateParticle(Hep, Vertex), PhotonA, PhotonB);	//passing reference to pointers
					vParticle.push_back(PhotonA);
					vParticle.push_back(PhotonB);
				}
				else if (abs(Hep->Pdg()) < 1e9)	//no nucleus
					vParticle.push_back(CreateParticle(Hep, Vertex));		//Everything detectable
			}			//sould contains muons, electron, pions, protons, kaons and other strange and charmed kaons

		}
		
		//Run through the collected particles
		CountParticles();
		ListCount();
		Identify();

		Out << std::endl;
	}

	//Data->Write();
	Global = ID;	//ID should be +1 
}

Particle* Background::CreateParticle(genie::GHepParticle *Hep, TVector3* Pos)
{
	TVector3 RefPos(Pos->X(),Pos->Y(),Pos->Z());
	TLorentzVector P4(*Hep->P4());
	Particle *P = new Particle(Hep->Pdg(), P4, RefPos);

	if (P->Charge() != 0)
		TheBox->SignalSmearing(GenMT, P);

	return P;
}


int Background::Count(std::string PartName, int N)
{
	if (N < -111)
		pCount[PartName] = 0;
	else pCount[PartName] += N;

	return pCount[PartName];
}

void Background::ListCount()
{
	if (pCount.size() == 0)
		Out << "No interesting particles " << std::endl;
	else
	{
		std::map<std::string, int>::iterator im = pCount.begin();
		for ( ; im != pCount.end(); ++im)
		{
			Out << im->first << ": " << im->second << "\t";
		}
		Out << std::endl;
	}
}

//Need to recalculate this factors..
bool Background::MuonOrPion(Particle *P)	//T if muon, F if pion
{
	if (P->TrackIn() < 1.0)
	{
	       if (GenMT->Rndm() < 0.3906445328)
		       return true;	
	       else return false;
	}
	else if (P->TrackIn() > 1.0 && P->TrackIn() < 2.0)
	{
	       if (GenMT->Rndm() < 0.549552089)
		       return true;	
	       else return false;
	}
	else if (P->TrackIn() > 2.0 && P->TrackIn() < 3.0)
	{
	       if (GenMT->Rndm() < 0.6748498302)
		       return true;	
	       else return false;
	}
	else if (P->TrackIn() > 3.0 && P->TrackIn() < 4.0)
	{
	       if (GenMT->Rndm() < 0.7482980534)
		       return true;	
	       else return false;
	}
	else if (P->TrackIn() > 4.0)
	{
	       if (GenMT->Rndm() < 0.7192982456)
		       return true;	
	       else return false;
	}
}

bool Background::IsPhoton(Particle *P)	//T if photon, F can be electron
{
	TVector3 Travel = (P->Position() - *Vertex);
	if (P->TrackIn() > 0.0 || Travel.Mag() > 0.03)	//displacement of 3cm
		return true;
	else return false;
}

bool Background::IseePair(Particle *P, Particle *Main)	//T if pair, F single electron
{
    	if (P->Direction().Angle(Main->Direction()) < 3.0/Const::fDeg)	//can't distinguish electrons 3deg close
		return true;
	else
		return false;
}

/***** Identification of single particles and counting *****/
bool Background::CountParticles()	//need to add smearing
{
	ParticleA = 0, ParticleB = 0;
	std::vector<Particle*> vElectron, vPhoton;
	std::vector<Particle*>::iterator ie;

	Out << "Particles selected: " << vParticle.size() << std::endl;

	int p = 0;
	for (iP = vParticle.begin(); iP != vParticle.end(); ++iP, ++p)
	{
		if (TheBox->IsDetectable(*iP) && TheBox->IsInside(*iP))
		{
			if ((*iP)->Pdg() == 13)		//Muon
			{
				if (!MuonOrPion(*iP))
				{
					(*iP)->SetPdg(211);
					(*iP)->SetMass(M_Pion);
					--iP;			//must recheck (for thr for instance)
					--p;
				}
				else 
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
				if (MuonOrPion(*iP))
				{
					(*iP)->SetPdg(13);
					(*iP)->SetMass(M_Muon);
					--iP;			//must recheck (for thr for instance)
					--p;
				}
				else 
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
				bool TrueElectron = true;
				for (ie = vElectron.begin(); ie != vElectron.end(); )
				{
					if (IseePair(*iP, *ie))
					{
						TVector3 Dir((*iP)->Px() + (*ie)->Px(),
							     (*iP)->Py() + (*ie)->Py(),
							     (*iP)->Pz() + (*ie)->Pz());
						TLorentzVector GammaReco(Dir, (*iP)->E() + (*ie)->E());
						(*iP)->SetPdg(22);	//could be a photon close to vertex
						(*iP)->SetP4(GammaReco);
						(*iP)->SetMass(0);
						(*iP)->SetTrackIn(abs((*iP)->TrackIn())+abs((*ie)->TrackIn()));

						--iP;			//must recheck (for thr for instance)
						--p;
						ie = vElectron.erase(ie);
						TrueElectron = false;
						Count("Electron", -1);
						break;
					}
					else ++ie;
				}

				if (TrueElectron && GenMT->Rndm() < TheBox->GetElement("ElecEfficiency"))
				{
					Count("Electron");
					vElectron.push_back(*iP);
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
			}
			else if ((*iP)->Pdg() == 22)	//Photons
			{
				if(!IsPhoton(*iP) || GenMT->Rndm() < TheBox->GetElement("GammaRejection"))
				{		//conversion before 3cm 
					(*iP)->SetPdg(11);
					(*iP)->SetMass(M_Electron);
					TheBox->TrackLength(GenMT, *iP);
					--iP;			//must recheck (for thr for instance)
					--p;
				}
				else
				{
					Count("Photon");
					switch (mapChan[GetChannel()])
					{
						case _nGAMMA:
							ParticleA = *iP;
							ParticleB = *iP;
							break;
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
		*/
		case _nGAMMA:
			if (IdentifynGAMMA())	//A = e, B = e
				LoadTree();
			break;
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
bool Background::IdentifynGAMMA()	//need to add smearing
{
	if (pCount["Electron"]	== 0	&&
	    pCount["Photon"]	== 1	&& 
	    pCount["Muon"]	== 0	&& 
	    pCount["Pion"]	== 0	&& 
	    pCount["Kaon"]	== 0	&& 
	    pCount["Charm"]	== 0	  )
		return true;
	else return false;
}

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
	
//Treat pi0 decay into 2 photons
void Background::Pi0Decay(Particle *Pi0, Particle *&PA, Particle *&PB)	//will travle length O(10nm)
{
	//in rest frame
	TLorentzVector GammaA(0, 0, M_Pion0/2.0, M_Pion0/2.0); 
	TLorentzVector GammaB(0, 0, -M_Pion0/2.0, M_Pion0/2.0); 

	TVector3 vBoost(Pi0->GetP4().BoostVector());
	TVector3 Start(Pi0->Position());		//starting point is vertex
	double Theta = GenMT->Uniform(-Const::fPi, Const::fPi);
	double Phi = GenMT->Uniform(-Const::fPi, Const::fPi);

	/*
	//double dx = TheBox->GetElement("Vertex");	//1 mm
	double dx = 1e-9;	//1 mm
	double Move = 0.0;
	Out << "before while " << dx/(Const::fC * Pi0->LabSpace()) <<  std::endl;
	while (GenMT->Rndm() < exp(-dx/(Const::fC * Pi0->LabSpace())))	//check every 1mm if particle has decayed
		Move += dx;
	Start *= (Start.Mag() + Move) / Start.Mag();
	Out << "after while " << Start.Mag() << std::endl;
	Out << "after while " << Move << std::endl;
	*/

	GammaA.SetTheta(Theta);
	GammaB.SetTheta(Const::fPi + Theta);
	GammaA.SetPhi(Phi);
	GammaB.SetPhi(Phi);

	GammaA.Boost(vBoost);
	GammaB.Boost(vBoost);

	TVector3 MoveA(GammaA.Vect().Unit());
	TVector3 MoveB(GammaB.Vect().Unit());
	MoveA *= GammaDecay();
	MoveB *= GammaDecay();
	MoveA += Start;
	MoveB += Start;

	PA = new Particle(22, GammaA, MoveA);	//here are the photons
	PB = new Particle(22, GammaB, MoveB);	//position should be different

	delete Pi0;
}

//Distance travelled by photons before conversion
double Background::GammaDecay()
{
	double Path = 9.0/7.0 * TheBox->InteractionLength(1);
	double dx = 0.001;	//1 mm
	double Length = 0.0;
	while (GenMT->Rndm() < exp(-dx/Path))	//check every 1mm if particle has decayed
		Length += dx;
	return Length;	//in cm
}
