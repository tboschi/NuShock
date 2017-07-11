#include "Background.h"


Background::Background(std::string EventDB, std::string DetectorConfig, int Nevt)	: //Decay rates calculator
	M_Neutrino(0.0),
	M_Photon(0.0),
	M_Electron(Const::fMElectron),
	M_Muon(Const::fMMuon),
	M_Pion(Const::fMPion),
	M_Pion0(Const::fMPion0),
	M_Kaon(Const::fMKaon),
	M_Kaon0(Const::fMKaon0)
{
	//GetCommandLineArgs (argc, argv);	//Useful?
	TTree *Gtree = 0;
	genie::NtpMCTreeHeader *Gthdr = 0;
 	TFile InFile(Channel.c_str(), "READ");
 	Gtree = dynamic_cast<TTree*>           (InFile.Get("gtree"));
 	Gthdr = dynamic_cast<NtpMCTreeHeader*> (InFile.Get("header"));
 	if(!Gtree)
	       return 1;

	genie::NtpMCEventRecord *Gmcrec = 0;
 	Gtree->SetBranchAddress("gmcrec", &Gmcrec);
 	// Get the nbr of evts to analyse (-n argument)
 	//NEV = (gOptNEvt > 0) ? TMath::Min(gOptNEvt, Gtree->GetEntries()) : (int)Gtree->GetEntries();
 	NEvt = Gtree->GetEntries();

	TheBox = new Detector(DetectorConfig);

	GenMT = new TRandom3;

	Data = new TTree("Data", "All data from GENIE");

	Data->Branch("ID", &ID, "iID/i");	//global event ID

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

void Background::Loop(std::string Channel)
{
	for (unsigned int i = 0; i < NEvt; ++i)
	{
		Gtree->GetEntry(i);
		EventRecord & Gevent = *(Gmcrec->event);
		GHepParticle * neu = Gevent.Probe();

		//TLorentzVector & neu4vec = *(neu->P4());
		//double total_xs = (XS_graph == 0) ? 1 : XS_graph->Eval(neu->E());
		//hinProbe->Fill(neu4vec.E());

		GHepParticle * p = 0;
		TIter event_iter(&Gevent);

		double PosX = GenMT->Uniform(TheBox->GetXsize());
		double PosY = GenMT->Uniform(TheBox->GetYsize());
		double PosZ = GenMT->Uniform(TheBox->GetZsize());
	
		vParticle.clear();	//reset particle array at each loop
		mapCount.clear();	//reset particle counting at every event
		while((CdHep = dynamic_cast<genie::GHepParticle *>(event_iter.Next())))	//loop on all particles
		{									//inside the event
			if (CdHep->Status() == 1)
			{
				if (CdHep->Pdg() == 211 ||		//special treatments for pi0
				    CdHep->Pdg() == -211)
				{
					Particle *PhotonA, *PhotonB;
					Pi0Decay(CreateParticle(CdHep, PosX, PosY, PosZ), PhotonA, PhotonB);
					vParticle.push_back(PhotonA);
					vParticle.push_back(PhotonB);
				}
				else if (CdHep->Pdg() < 1000000000 && 	//no nucleus
				         CdHep->Pdg() > -1000000000) 	//no nucleus
					vParticle.push_back(CreateParticle(CdHep, PosX, PosY, PosZ));		//Everything detectable
			}			//sould contains muons, electron, pions, protons
		}

		Purify(Channel);
		Identify(Channel, i);

		//Ready for output
	}
}

Particle* Background::CreateParticle(GHepParticle *Hep, double PosX, double PosY, double PosZ)
{
	Particle *P = new Particle(Hep->Pdg(), Hep->Charge(), Hep->P4(), PosX, PosY, PosZ);

	TheBox->SignalSmearing(GenMT, P);

	return P;
}

int Background::Count(std::string PartName, int N)
{
	if (N < 0)
		mapCount[PartName] = 0;
	else mapCount[PartName] += N;

	return mapCount[PartName];
}

/***** Purification of events, will also add background contribution ********/
void Background::Purify(std::string Channel)
{
	for (iP = vParticle.begin(); iP != vParticle.end(); )
	{
		//if(PurifyHadron(*iP))	//check if hadron first, if so I can remove it
		//	iP = vParticle.erase(iP);

		//else switch(mapChan[Channel])
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
				if (!PurifynEE(*iP))
					iP = vParticle.erase(iP);
				else ++iP;
				break;
			case _nEMU:
				if (!PurifynEMU(*iP))
					iP = vParticle.erase(iP);
				else ++iP;
				break;
			case _nMUE:
				if (!PurifynEMU(*iP))
					iP = vParticle.erase(iP);
				else ++iP;
				break;
			case _nPI0:
				if (!PurifynPI0(*iP))
					iP = vParticle.erase(iP);
				else ++iP;
				break;
			case _EPI:
				Result = EPI();
				if (!PurifyEPI(*iP))
					iP = vParticle.erase(iP);
				else ++iP;
				break;
			case _nMUMU:
				Result = nMUMU();
				if (!PurifynMUMU(*iP))
					iP = vParticle.erase(iP);
				else ++iP;
				break;
			case _MUPI:
				if (!PurifyMUPI(*iP))
					iP = vParticle.erase(iP);
				else ++iP;
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
}

bool Background::PurifyHadron(Particle *iP)	//here everything not selected by cut
{
	if (iP->E() > TheBox->GetElement("ThrHadron"))	//hadron (proton) activity
	{
		if (iP->PdgCode() == 311)	//Count kaons
			Count("Kaon", 1);
		else if (iP->Charge() != 0 && iP->PdgCode() > 100)		//I can see charged hadrons
			Count("Hadron", 1);	//count the protons // will need reliable efficiencies
		else Count("Uknown", 1);
	}
	else Count("UnderThr", 1);
}

/*
bool Background::PurifynEE(Particle *iP)
{
	bool Ret = true;

	if (iP->PdgCode == 11)		//Muon
	{
		if (GenMT->Rndm() < Alpha)
			Ret = false;
		else Count("Electron", 1);
	}
	else Ret = false;

	return Ret;
}

bool Background::PurifynEMU(Particle *iP)
{
	bool Ret = true;

	if (iP->PdgCode == 11)		//Electron
	{
		if (GenMT->Rndm() < Alpha)
			Ret = false;
		else Count("Electron", 1);
	}
	else if (iP->PdgCode == 13)	//Muon
	{
		if (GenMT->Rndm() < Alpha)
			Ret = false;
		else Count("Muon", 1);
	}
	else Ret = false;

	return Ret;
}

bool Background::PurifynPI0(Particle *iP)
{
	bool Ret = true;

	if (iP->PdgCode == 11)		//Electron
	{
		if (GenMT->Rndm() < Alpha)
			Ret = false;
		else Count("Electron", 1);
	}
	else if (iP->PdgCode == 111)	//Pion0
	{
		if (GenMT->Rndm() < Alpha)
			Ret = false;
		else Count("Pion0", 1);
	}
	else if (iP->PdgCode == 22)	//Photon
	{
		if (GenMT->Rndm() < Alpha)
			Ret = false;
		else Count("Photon", 1);
	}
	else Ret = false;

	return Ret;
}

bool Background::PurifyEPI(Particle *iP)	//Possible background is photon?
{
	bool Ret = true;

	if (iP->PdgCode == 11)		//Electron
	{
		if (GenMT->Rndm() < Alpha)
			Ret = false;
		else Count("Electron", 1);
	}
	else if (iP->PdgCode == 211)	//Pion+
	{
		if (GenMT->Rndm() < Alpha)
			Ret = false;
		else Count("Pion", 1);
	}
	else if (iP->PdgCode == 111) 	//Pion0 to be misidentified
	{
		Particle *PhotonA = new Particle();
		Particle *PhotonB = new Particle();
		if (GenMT->Rndm() < Alpha)
			Ret = false;
		else Count("Photon", 1);
	}
	else Ret = false;

	return Ret;
}

bool Background::PurifynMUMU(Particle *iP)
{
	bool Ret = true;

	if (iP->PdgCode == 13)		//Muon
	{
		if (GenMT->Rndm() < TheBox->GetElement("ThrMuon"));	//can't detect
			Ret = false;
		else Count("Muon", 1);
	}
	else if (iP->PdgCode == 211)	//Pion+
	{
		if (iP->E() < TheBox->GetElement("ThrPion"))
		       Ret = false;
		else if (TheBox->TrackLength(iP) > 0.50)	//track length is muon-like
		{
			iP->SetPdg(13);
			Count("Muon", 1);	//thus I count a muon
		}
	}
	else Ret = false;

	return Ret;
}
*/
bool Background::PurifyMUPI(Particle *iP)	//need to add smearing
{
	bool Ret = true;

	if (iP->PdgCode() == 13)		//Muon
	{
		if (iP->E() < TheBox->GetElement("ThrMuon"));	//can't detect
			Ret = false;
		else
		{
			if (TheBox->TrackLength(iP) < 0.50)
			{
				iP->SetPdg(211);
				iP->SetMass(M_Pion);
				Count("Pion", 1);
			}
			else Count("Muon", 1);	//is really a muon
		}
	}
	else if (iP->PdgCode() == 211)	//Pion+
	{
		if (iP->E() < TheBox->GetElement("ThrPion"));	//can't detect
			Ret = false;
		else
		{
			if (TheBox->TrackLength(iP) > 0.50)	//track length is muon-like
			{
				iP->SetPdg(13);
				iP->SetMass(M_Muon);
				Count("Muon", 1);
			}
			else Count("Pion", 1);	//is really a pion
	}
	else 	//check if is charged hadron
	{
		PurifyHadron(iP);
		Ret = false;
	}

	return Ret;
}

/***** Idenfitication of events, by counting the tracks ********/
bool Background::Identify(std::string Channel, int NEvt)
{
	switch(mapChan[Channel])
	{
		/*
		case _ALL:
			break;
		case _nnn:
			break;
		case _nGAMMA:
			break;
		*/
		case _nEE:
			return IdentifynEE();
		case _nEMU:
			return IdentifynEMU();
		case _nMUE:
			return IdentifynMUE();
		case _nPI0:
			return IdentifynPI0();
		case _EPI:
			return IdentifyEPI();
		case _nMUMU:
			return IdentifynMUMU();
		case _MUPI:
			if (IdentifyMUPI())
				LoadTree(NEvt, 13, 211);	//Load the mu and the pi
			break;
		/*
		case _EKA:
			break;
		case _nKA0:
			break;
		*/
	}
}

/*
bool Background::IdentifynEE()
{
	bool Ret = true;

	std::map<std::string, int>::iterator im = mapCount.begin();
	for ( ; im != mapCount.end(); ++im)
	{
		if (im->first == "Hadron")	//keep hadrons
			continue;

		if (im->second == 2)		//Looking for this: 2 e
		{
			if (im->first != "Electron")
				Ret = false;
			//else is true
		}
		else if (im->second == 0)
		{
			if (im->first == "Electron")
				Ret = false;
			//else is true
		}
		else Ret = false;
	}

	return Ret;
}

bool Background::IdentifynEMU()
{
	bool Ret = true;

	std::map<std::string, int>::iterator im = mapCount.begin();
	for ( ; im != mapCount.end(); ++im)
	{
		if (im->first == "Hadron")	//keep hadrons
			continue;

		if (im->second == 1)
		{
			if (im->first == "Muon")
				Ret *= true;
			else if (im->first == "Pion")
				Ret *= true;
			else false;
		}
		else Ret = false;
	}

	return Ret;
}

bool Background::IdentifynPI0()
{
	bool Ret = true;

	std::map<std::string, int>::iterator im = mapCount.begin();
	for ( ; im != mapCount.end(); ++im)
	{
		if (im->second == 1)
		{
			if (im->first == "Muon")
				Ret *= true;
			else if (im->first == "Pion")
				Ret *= true;
			else false;
		}
		else Ret = false;
	}

	return Ret;
}

bool Background::IdentifyEPI()	//Possible background is photon?
{
	bool Ret = true;

	std::map<std::string, int>::iterator im = mapCount.begin();
	for ( ; im != mapCount.end(); ++im)
	{
		if (im->second == 1)
		{
			if (im->first == "Muon")
				Ret *= true;
			else if (im->first == "Pion")
				Ret *= true;
			else false;
		}
		else Ret = false;
	}

	return Ret;
}

bool Background::IdentifynMUMU()
{
	bool IsOther = false, IsMuon = false;

	std::map<std::string, int>::iterator im = mapCount.begin();
	for ( ; im != mapCount.end(); ++im)
	{
		if (im->first == "Muon")	//Look for 2 muon
		{
			if (im->second == 2)
				IsMuon = true;
			else IsMuon = false;
		}
		else if (im->first == "Nuclear")	//what should I do?
		{
		}
		else 
		{
			if (im->second > 0)
				IsOther += true;
		}
	}

	return !IsOther * IsMuon; 	//0 others, 1 pi, 2 mu
}
*/
bool Background::IdentifyMUPI()	//need to add smearing
{
	bool IsOther = false, IsMuon = false, IsPion = false;

	std::map<std::string, int>::iterator im = mapCount.begin();
	for ( ; im != mapCount.end(); ++im)
	{
		if (im->first == "Muon")	//Look for 1 muon
		{
			if (im->second == 1)
				IsMuon = true;
			else IsMuon = false;
		}
		else if (im->first == "Pion")	//Look for 1 pion
		{
			if (im->second == 1)
				IsPion = true;
			else IsPion = false;
		}
		else if (im->first == "Nuclear")	//what should I do?
		{
		}
		else 
		{
			if (im->second > 0)
				IsOther += true;
		}
	}

	return !IsOther * IsMuon * IsPion; 	//0 others, 1 pi, 1 mu
}

void Background::LoadTree(int idEvent, int pdgA, int pdgB)
{
	ID = idEvent; 

	TLorentzVector InA, InB;
	for (iP = vParticle.begin(); iP != vParticle.end(); ++iP)
	{
		if (iP->PdgCode() == pdgA)
		{
			EnergyA = iP->E();
			MomentA = iP->P();
			TransvA = iP->Pt();
			ThetaA = iP->Theta();
			PhiA = iP->Phi();
			MassA = iP->M();
			InA = iP->GetP4();
		}
		else if (iP->PdgCode() == pdgB)
		{
			EnergyB = iP->E();
			MomentB = iP->P();
			TransvB = iP->Pt();
			ThetaB = iP->Theta();
			PhiB = iP->Phi();
			MassB = iP->M();
			InB = iP->GetP4();
		}
	}

	TLorentzVector Reco = InA+InB;

	Energy0 = Reco.E();
        Moment0 = Reco.P();
        Transv0 = Reco.Pt();
        Theta0 = Reco.Theta(
        Phi0 = Reco.Phi();
        Mass0 = Reco.M();

	Data->Fill();
}


//Treat pi0 decay into 2 photons
void Background::Pi0Decay(Particle *Pi0, Particle *PA, Paticle *PB)
{
	TLorentzVector GammaA(0, 0, M_Pion0/2.0, M_Pion0/2.0); 
	TLorentzVector GammaA(0, 0, -M_Pion0/2.0, M_Pion0/2.0); 

	TVector3 vBoost(Pi0->GetP4()->BoostVector());
	double Theta = GenMT->Uniform(-Const::fPi, Const::fPi);
	double Phi = GenMT->Uniform(-Const::fPi, Const::fPi);

	GammaA->SetTheta(Theta);
	GammaB->SetTheta(Const::fPi + Theta);
	GammaA->SetPhi(Phi);
	GammaB->SetPhi(Phi);

	GammaA.Boost(vBoost);
	GammaB.Boost(vBoost);

	PA = new Particle(22, GammaA, PosX, PosY, PosZ);	//here are the photons
	PB = new Particle(22, GammaB, PosX, PosY, PosZ);	//position should be different

	delete Pi0;
}
