#include "Cut.h"


Cut::Cut(std::string EventDB, std::string DetectorConfig, int Nevt)	: //Decay rates calculator
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
 	NEV = Gtree->GetEntries();

	TheBox = new Detector(DetectorConfig);

	GenMT = new TRandom3;

	Alpha = 0.05;
}

//Initialisation of map
void Cut::MapInit()
{
	mapCut["ALL"] = _ALL;
	mapCut["nnn"] = _nnn;
	mapCut["nGAMMA"] = _nGAMMA;
	mapCut["nEE"] = _nEE;
	mapCut["nEMU"] = _nEMU;
	mapCut["nMUE"] = _nMUE;
	mapCut["nPI0"] = _nPI0;
	mapCut["EPI"] = _EPI;
	mapCut["nMUMU"] = _nMUMU;
	mapCut["MUPI"] = _MUPI;
	mapCut["EKA"] = _EKA;
	mapCut["nKA0"] = _nKA0;
}

void Cut::Loop(std::string Channel)
{
	for (unsigned int i = 0; i < NEV; ++i)
	{
		Gtree->GetEntry(i);
		EventRecord & Gevent = *(Gmcrec->event);
		GHepParticle * neu = Gevent.Probe();

		TLorentzVector & neu4vec = *(neu->P4());
		double total_xs = (XS_graph == 0) ? 1 : XS_graph->Eval(neu->E());
		hinProbe->Fill(neu4vec.E());

		GHepParticle * p = 0;
		TIter event_iter(&Gevent);

		vParticle.clear();	//reset particle array at each loop
		mapCount.clear();	//reset particle counting at every event
		while((Cd = dynamic_cast<genie::GHepParticle *>(event_iter.Next())))	//loop on all particles
		{									//inside the event
			double PosX = GenMT->Uniform(TheBox->GetXsize());
			double PosY = GenMT->Uniform(TheBox->GetYsize());
			double PosZ = GenMT->Uniform(TheBox->GetZsize());
	
			if (Cd->Status() == 1)
			{
				if (Cd->Pdg() == 211 ||		//special treatments for pi0
				    Cd->Pdg() == -211)
				{
					Particle *PhotonA, *PhotonB;
					Pi0Decay(CreateParticle(Cd->P4(), Cd->Pdg(), PosX, PosY, PosZ), PhotonA, PhotonB);
					vParticle.push_back(PhotonA);
					vParticle.push_back(PhotonB);
				}
				else if (Cd->Pdg() < 1000000000 && 	//not nucleus
				         Cd->Pdg() > -1000000000 && 	//not nucleus
				         Cd->Pdg() != 2112)		//not neutron
					vParticle.push_back(CreateParticle(Cd->P4(), Cd->Pdg(), PosX, PosY, PosZ));		//Everything detectable
			}			//sould contains muons, electron, pions, protons
		}

		Purify(Channel);
		Identify(Channel);

		//Ready for output
	}
}

Particle* Cut::CreateParticle(TLorentzVector *N, int Pdg, double PosX, double PosY, double PosZ)
{
	Particle *P = new Particle(Pdg, N, PosX, PosY, PosZ);

	TheBox->SignalSmearing(GenMT, P);

	return P;
}

int Cut::Count(std::string PartName, int N)
{
	if (N < 0)
		mapCount[PartName] = 0;
	else mapCount[PartName] += N;

	return mapCount[PartName];
}

/***** Purification of events, will also add background contribution ********/
void Cut::Purify(std::string Channel)
{
	for (iP = vParticle.begin(); iP != vParticle.end(); )
	{
		if ((*iP)->E() < TheBox->GetElement("ThrHadron"))	//hadron (proton) activity
		{
			if ((*iP)->PdgCode == 2212 ||
			    (*iP)->PdgCode == 22222222 )
			Count("Hadron", 1);	//count the protons // will need reliable efficiencies
		}
		else switch(GetChannel())
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
				break;
		}
	}
}

/*
bool Cut::PurifynEE(Particle *iP)
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

bool Cut::PurifynEMU(Particle *iP)
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

bool Cut::PurifynPI0(Particle *iP)
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

bool Cut::PurifyEPI(Particle *iP)	//Possible background is photon?
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

bool Cut::PurifynMUMU(Particle *iP)
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
bool Cut::PurifyMUPI(Particle *iP)	//need to add smearing
{
	bool Ret = true;

	if (iP->PdgCode == 13)		//Muon
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
	else if (iP->PdgCode == 211)	//Pion+
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
	else Ret = false;

	return Ret;
}

/***** Idenfitication of events, by counting the tracks ********/
bool Cut::Identify(std::string Channel)
{
	bool Ret = true;

	std::map<std::string, int>::iterator im = mapCount.begin();

	while (Ret || im != mapCount.end())
	{
		if (im->first != "Hadron")	//keep hadrons	
			switch(GetChannel())
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
					if (!IdentifynEE(im->first, im->second))
						Ret = false;
					break;
				case _nEMU:
					if (!IdentifynEMU(im->first, im->second))
						Ret = false;
					break;
				case _nMUE:
					if (!IdentifynEMU(im->first, im->second))
						Ret = false;
					else ++iP;
					break;
				case _nPI0:
					if (!IdentifynPI0(im->first, im->second))
						Ret = false;
					break;
				case _EPI:
					if (!IdentifyEPI(im->first, im->second))
						Ret = false;
					else ++iP;
					break;
				case _nMUMU:
					if (!IdentifynMUMU(im->first, im->second))
						Ret = false;
					break;
				case _MUPI:
					if (!IdentifyMUPI(im->first, im->second))
						Ret = false;
					break;
				/*
				case _EKA:
					break;
				case _nKA0:
					break;
				*/
			}

			++im;	//increment pointer
		}
	}	
}
/*
bool Cut::IdentifynEE(Particle *iP)
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

bool Cut::IdentifynEMU(Particle *iP)
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

bool Cut::IdentifynPI0(Particle *iP)
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

bool Cut::IdentifyEPI(Particle *iP)	//Possible background is photon?
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

bool Cut::IdentifynMUMU(Particle *iP)
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
*/
bool Cut::IdentifyMUPI(std::string PartName, int Num)	//need to add smearing
{
	bool Ret = true;
	bool IsMuon, IsPion;

	std::map<std::string, int>::iterator im = mapCount.begin();
	for ( ; im != mapCount.end(); ++im)
	{
		if (im->second == 1)
		{
			if (im->first == "Muon") 	//1 Muon : OK
				IsMuon = true;
			else if (im->first == "Pion")	//1 Pion : OK
				IsPion = true;
			else Ret = false;		//1 other : X
		}
		else if (im->second == 0)
		{
			if (im->first != "Muon" && im->first != "Pion") //0 other : OK
				Ret = true;
			else false;			//0 Muon or 0 Pion : X
		}
		else Ret = false;
	}

	return Ret*IsMuon*IsPion;
}


//Treat pi0 decay into 2 photons
void Cut::Pi0Decay(Particle *Pi0, Particle *PA, Paticle *PB)
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
