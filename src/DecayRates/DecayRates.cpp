#include "DecayRates.h"

Decay::Decay(double MSterile, double Ue, double Um, double Ut)	: //Decay rates calculator
	M_Neutrino(0.0),
	M_Photon(0.0),
	M_Electron(Const::fMElectron),
	M_Muon(Const::fMMuon),
	M_Pion(Const::fMPion),
	M_Pion0(Const::fMPion0),
	M_Kaon(Const::fMKaon),
	M_Kaon0(Const::fMKaon0),
	M_Eta(Const::fMEta),
        M_Rho(Const::fMRho),
        M_Rho0(Const::fMRho0),
        M_Omega(Const::fMOmega),
        M_Kaonx(Const::fMKaonx),
        M_Kaon0x(Const::fMKaon0x),
        M_Etai(Const::fMEta),
        M_Phi(Const::fMPhi),
        M_Tau(Const::fMTau),
        M_Charm(Const::fMCharm)
{
	M_Sterile_prev = -1.0;

	TheSpace = new ThreeBody("");

	SetMass(MSterile);
	SetUe(Ue);
	SetUm(Um);
	SetUt(Ut);

	MapInit();
	SetEnhancement();
	Event = new TGenPhaseSpace;
	N_vec = new TLorentzVector;
	N_rest = new TLorentzVector;
}

//Initialisation of map
void Decay::MapInit()
{
	mapChannel["ALL"]    = _ALL;

	//unclassified
	mapChannel["nnn"]    = _nnn;
	mapChannel["nGAMMA"] = _nGAMMA;

	//pure leptonic
	mapChannel["nEE"]    = _nEE;
	mapChannel["nEMU"]   = _nEMU;
	mapChannel["nMUE"]   = _nMUE;
	mapChannel["nMUMU"]  = _nMUMU;
	mapChannel["nET"]    = _nET;
	mapChannel["nTE"]    = _nTE;
	mapChannel["nMUT"]   = _nMUT;
	mapChannel["nTMU"]   = _nTMU;

	//pion
	mapChannel["nPI0"]   = _nPI0;
	mapChannel["EPI"]    = _EPI;
	mapChannel["MUPI"]   = _MUPI;
	mapChannel["TPI"]    = _TPI;

	//kaon
	mapChannel["EKA"]    = _EKA;
	mapChannel["MUKA"]   = _MUKA;

	//rho
	mapChannel["nRHO0"]  = _nRHO0;
	mapChannel["ERHO"]   = _ERHO;
	mapChannel["MURHO"]  = _MURHO;

	//kaon*
	mapChannel["EKAx"]   = _EKAx;
	mapChannel["nKA0x"]  = _nKA0x;
	mapChannel["MUKAx"]  = _MUKAx;

	//other (eta, phi, omega.. )
	mapChannel["nETA"]   = _nETA;
	mapChannel["nETAi"]  = _nETAi;
	mapChannel["nOMEGA"] = _nOMEGA;
	mapChannel["nPHI"]   = _nPHI;

	//charm
	mapChannel["ECHARM"] = _ECHARM;

}

double Decay::Gamma(std::string Channel, double B)
{
	//SetEnhancement(Channel, B);

	double Result;
	switch(mapChannel[Channel])
	{
		case _ALL:
			Result = Total();
			break;
		case _nnn:
			Result = nnn();
			break;
		case _nGAMMA:
			Result = nGAMMA();
			break;
		case _nEE:
			Result = nEE();
			break;
		case _nEMU:
			Result = nEMU();
			break;
		case _nMUE:
			Result = nMUE();
			break;
		case _nMUMU:
			Result = nMUMU();
			break;
		case _nET:
			Result = nET();
			break;
		case _nTE:
			Result = nTE();
			break;
		case _nMUT:
			Result = nMUT();
			break;
		case _nTMU:
			Result = nTMU();
			break;
		case _nPI0:
			Result = nPI0();
			break;
		case _EPI:
			Result = EPI();
			break;
		case _MUPI:
			Result = MUPI();
			break;
		case _TPI:
			Result = TPI();
			break;
		case _EKA:
			Result = EKA();
			break;
		case _MUKA:
			Result = MUKA();
			break;
		case _EKAx:
			Result = EKAx();
			break;
		case _nKA0x:
			Result = nKA0x();
			break;
		case _MUKAx:
			Result = MUKAx();
			break;
		case _nRHO0:
			Result = nRHO0();
			break;
		case _ERHO:
			Result = ERHO();
			break;
		case _MURHO:
			Result = MURHO();
			break;
		case _nETA:
			Result = nETA();
			break;
		case _nETAi:
			Result = nETAi();
			break;
		case _nOMEGA:
			Result = nOMEGA();
			break;
		case _nPHI:
			Result = nPHI();
			break;
		case _ECHARM:
			Result = ECHARM();
			break;
		default:
			Result = -1.0;
			break;
	}

	//SetEnhancement();
	return Result;
}

double Decay::Other(std::string Channel, double A)
{
	//SetEnhancement("ALL", A);

	double Result;
	switch(mapChannel[Channel])
	{
		case _ALL:
			Result = 0.0;
			break;
		case _nnn:
			Result = Total() - nnn();
			break;
		case _nGAMMA:
			Result = Total() - nGAMMA();
			break;
		case _nEE:
			Result = Total() - nEE();
			break;
		case _nEMU:
			Result = Total() - nEMU();
			break;
		case _nMUE:
			Result = Total() - nMUE();
			break;
		case _nMUMU:
			Result = Total() - nMUMU();
			break;
		case _nET:
			Result = Total() - nET();
			break;
		case _nTE:
			Result = Total() - nTE();
			break;
		case _nMUT:
			Result = Total() - nMUT();
			break;
		case _nTMU:
			Result = Total() - nTMU();
			break;
		case _nPI0:
			Result = Total() - nPI0();
			break;
		case _EPI:
			Result = Total() - EPI();
			break;
		case _MUPI:
			Result = Total() - MUPI();
			break;
		case _TPI:
			Result = Total() - TPI();
			break;
		case _EKA:
			Result = Total() - EKA();
			break;
		case _MUKA:
			Result = Total() - MUKA();
			break;
		case _EKAx:
			Result = Total() - EKAx();
			break;
		case _nKA0x:
			Result = Total() - nKA0x();
			break;
		case _MUKAx:
			Result = Total() - MUKAx();
			break;
		case _nRHO0:
			Result = Total() - nRHO0();
			break;
		case _ERHO:
			Result = Total() - ERHO();
			break;
		case _MURHO:
			Result = Total() - MURHO();
			break;
		case _nETA:
			Result = Total() - nETA();
			break;
		case _nETAi:
			Result = Total() - nETAi();
			break;
		case _nOMEGA:
			Result = Total() - nOMEGA();
			break;
		case _nPHI:
			Result = Total() - nPHI();
			break;
		case _ECHARM:
			Result = Total() - ECHARM();
			break;
		default:
			Result = -1.0;
			break;
	}

	//SetEnhancement();
	return Result;
}

double Decay::Branch(std::string Channel, double A, double B)
{
	//SetEnhancement("ALL", A);
	//SetEnhancement(Channel, B);

	double Result;
	switch(mapChannel[Channel])
	{
		case _ALL:
			Result = 1.0;
			break;
		case _nnn:
			Result = nnn()/Total();
			break;
		case _nGAMMA:
			Result = nGAMMA()/Total();
			break;
		case _nEE:
			Result = nEE()/Total();
			break;
		case _nEMU:
			Result = nEMU()/Total();
			break;
		case _nMUE:
			Result = nMUE()/Total();
			break;
		case _nMUMU:
			Result = nMUMU()/Total();
			break;
		case _nET:
			Result = nET()/Total();
			break;
		case _nTE:
			Result = nTE()/Total();
			break;
		case _nMUT:
			Result = nMUT()/Total();
			break;
		case _nTMU:
			Result = nTMU()/Total();
			break;
		case _nPI0:
			Result = nPI0()/Total();
			break;
		case _EPI:
			Result = EPI()/Total();
			break;
		case _MUPI:
			Result = MUPI()/Total();
			break;
		case _TPI:
			Result = TPI()/Total();
			break;
		case _EKA:
			Result = EKA()/Total();
			break;
		case _MUKA:
			Result = MUKA()/Total();
			break;
		case _EKAx:
			Result = EKAx()/Total();
			break;
		case _nKA0x:
			Result = nKA0x()/Total();
			break;
		case _MUKAx:
			Result = MUKAx()/Total();
			break;
		case _nRHO0:
			Result = nRHO0()/Total();
			break;
		case _ERHO:
			Result = ERHO()/Total();
			break;
		case _MURHO:
			Result = MURHO()/Total();
			break;
		case _nETA:
			Result = nETA()/Total();
			break;
		case _nETAi:
			Result = nETAi()/Total();
			break;
		case _nOMEGA:
			Result = nOMEGA()/Total();
			break;
		case _nPHI:
			Result = nPHI()/Total();
			break;
		case _ECHARM:
			Result = ECHARM()/Total();
			break;
		default:
			Result = -1.0;
			break;
	}

	//SetEnhancement();
	return Result;
}

int Decay::PhaseSpace(std::string Channel, double &Weight)	//Return number of products 
{								//0 if decay not valid
	//SetEnhancement();
	TheSpace->SetParent(Channel);

	double Mass[3];
	int Products;

	switch(mapChannel[Channel])
	{
		/* Invisible channels
		case _nnn:
			Mass[0] = M_Neutrino;
			Mass[1] = M_Neutrino;
			Mass[2] = M_Neutrino;
			Products = 3 * Event->SetDecay(*N_rest, 3, Mass);
			Weight = Event->Generate();
			break;
		*/		
		case _nGAMMA:
			Mass[0] = M_Photon;
			Mass[1] = M_Neutrino;
			PdgCode[0] = 22;
			PdgCode[1] = 12;
			Products = 2 * Event->SetDecay(*N_rest, 2, Mass);
			Weight = Event->Generate();
			break;
		
		case _nEE:
			Mass[0] = M_Electron;
			Mass[1] = M_Electron;
			Mass[2] = M_Neutrino;
			PdgCode[0] = 11;
			PdgCode[1] = 11;
			PdgCode[2] = 12;
			Products = 3 * Event->SetDecay(*N_rest, 3, Mass);
			Weight = Event->Generate();

			TheSpace->SetEnergyX(Event->GetDecay(0)->E());
			TheSpace->SetEnergyY(Event->GetDecay(1)->E());
			Weight = TheSpace->ddGamma()/TheSpace->MaxGamma();
			break;

		case _nEMU:		//whata about n mu e?
		case _nMUE:		//whata about n mu e?
			Mass[0] = M_Muon;
			Mass[1] = M_Electron;
			Mass[2] = M_Neutrino;
			PdgCode[0] = 13;
			PdgCode[1] = 11;
			PdgCode[2] = 12;
			Products = 3 * Event->SetDecay(*N_rest, 3, Mass);
			Weight = Event->Generate();

			TheSpace->SetEnergyX(Event->GetDecay(0)->E());
			TheSpace->SetEnergyY(Event->GetDecay(1)->E());
			Weight = TheSpace->ddGamma()/TheSpace->MaxGamma();
			break;
	
		case _nPI0:
			Mass[0] = M_Pion0;
			Mass[1] = M_Neutrino;		//Invisible particle
			PdgCode[0] = 111;
			PdgCode[1] = 12;
			Products = 2 * Event->SetDecay(*N_rest, 2, Mass);
			Weight = Event->Generate();
			break;

		case _EPI:
			Mass[0] = M_Electron;
			Mass[1] = M_Pion;
			PdgCode[0] = 11;
			PdgCode[1] = 211;
			Products = 2 * Event->SetDecay(*N_rest, 2, Mass);
			Weight = Event->Generate();
			break;

		case _nMUMU:
			Mass[0] = M_Muon;
			Mass[1] = M_Muon;
			Mass[2] = M_Neutrino;
			PdgCode[0] = 13;
			PdgCode[1] = 13;
			PdgCode[2] = 12;
			Products = 3 * Event->SetDecay(*N_rest, 3, Mass);
			Weight = Event->Generate();

			TheSpace->SetEnergyX(Event->GetDecay(0)->E());
			TheSpace->SetEnergyY(Event->GetDecay(1)->E());
			Weight = TheSpace->ddGamma()/TheSpace->MaxGamma();
			break;

		case _MUPI:
			Mass[0] = M_Muon;
			Mass[1] = M_Pion;
			PdgCode[0] = 13;
			PdgCode[1] = 211;
			Products = 2 * Event->SetDecay(*N_rest, 2, Mass);
			Weight = Event->Generate();
			break;

		case _EKA:
			Mass[0] = M_Electron;
			Mass[1] = M_Kaon;
			PdgCode[0] = 11;
			PdgCode[1] = 321;
			Products = 2 * Event->SetDecay(*N_rest, 2, Mass);
			Weight = Event->Generate();
			break;

		default:
			Products = 0;
			break;
	}

	return Products;
}

//Controller of decay enhancement
void Decay::SetEnhancement(std::string Channel, double K)
{
	std::map<std::string, ChannelName>::iterator it;
	if (Channel == "ALL")
		for (it = mapChannel.begin(); it != mapChannel.end(); ++it)
			mapEnhance[it->first] = K;
	else
	{
		it = mapChannel.find(Channel);
		if (it != mapChannel.end())
			mapEnhance[it->first] = K;
	}
}


//total decay width
double Decay::Total(double A)
{
	return A * ( nnn() + nGAMMA() +
		     nEE() + nEMU() + nMUE() + nMUMU() + nET() + nTE() + nMUT() + nTMU() +
		     nPI0() + EPI() + + MUPI() + TPI() +
		     EKA() + MUKA() + 
		     nRHO0() + ERHO() + MURHO() +
		     EKAx() + nKA0x() + MUKAx() + 
		     nETA() + nETAi() + nOMEGA() + nPHI() +
		     ECHARM()	);
}


//individual decay channels
//all mixing factors are factorised out
double Decay::nnn()
{
	if (fnnn < 0 || IsChanged())
	{
		if (M_Sterile >= 3 * M_Neutrino)
		{
			fnnn = mapEnhance["nnn"] * Const::fGF2 * pow(M_Sterile, 5) /
			       (96.0 * Const::fPi3);
		}
		else fnnn = 0.0;
	}

	return fnnn * (Ue*Ue + Um*Um + Ut*Ut);
}

double Decay::nGAMMA()
{
	if (fnGAMMA < 0 || IsChanged())
	{
		if (M_Sterile >= M_Neutrino + M_Photon)
		{
			double AemPi = Const::fAem / Const::fPi;
			fnGAMMA = mapEnhance["nGAMMA"] * Const::fGF2 * pow(M_Sterile, 5) *
			       (27.0/32.0 * AemPi) / (192.0 * Const::fPi3);
		}
		else fnGAMMA = 0.0;
	}

	return fnGAMMA * (Ue*Ue + Um*Um + Ut*Ut);
}

//M_Sterile > 2 M_Electron (always)
double Decay::nEE()
{
	if (fnEE_e < 0 || fnEE_mt < 0 || IsChanged())
		NeutrinoLeptonAA(fnEE_e, fnEE_mt, M_Electron);

	return mapEnhance["nEE"] * fnEE_e * Ue*Ue + fnEE_mt * (Um*Um + Ut*Ut);
}

//M_Sterile > M_Muon + M_Electron
double Decay::nEMU()	//Antiparticle is Elec
{
	if (fnEMU < 0 || IsChanged())
		fnEMU = NeutrinoLeptonAB(M_Muon, M_Electron);

	return mapEnhance["nEMU"] * fnEMU * Ue*Ue;
}

double Decay::nMUE()	//Anti is Muon
{
	if (fnMUE < 0 || IsChanged())
		fnMUE = NeutrinoLeptonAB(M_Electron, M_Muon);

	return mapEnhance["nMUE"] * fnMUE * Um*Um;
}

//M_Sterile > 2 M_Muon
double Decay::nMUMU()
{
	if (fnMUMU_m < 0 || fnMUMU_et < 0 || IsChanged())
		NeutrinoLeptonAA(fnMUMU_m, fnMUMU_et, M_Muon);

	return mapEnhance["nMUMU"] * (fnMUMU_m * Um*Um + fnMUMU_et * (Ue*Ue + Ut*Ut));
}

//M_Sterile > M_Tau + M_Electron
double Decay::nET()	//Antiparticle is Elec
{
	if (fnET < 0 || IsChanged())
		fnET = NeutrinoLeptonAB(M_Tau, M_Electron);

	return mapEnhance["nET"] * fnET * Ut*Ut;
}

double Decay::nTE()	//Anti is Tau
{
	if (fnTE < 0 || IsChanged())
		fnTE = NeutrinoLeptonAB(M_Electron, M_Tau);

	return mapEnhance["nTE"] * fnTE * Ue*Ue;
}

//M_Sterile > M_Tau + M_Muon
double Decay::nMUT()	//Antiparticle is Muon
{
	if (fnMUT < 0 || IsChanged())
		fnMUT = NeutrinoLeptonAB(M_Tau, M_Muon);

	return mapEnhance["nMUT"] * fnMUT * Ut*Ut;
}

double Decay::nTMU()	//Anti is Tau
{
	if (fnTMU < 0 || IsChanged())
		fnTMU = NeutrinoLeptonAB(M_Muon, M_Tau);

	return mapEnhance["nTMU"] * fnTE * Um*Um;
}

//M_Sterile > M_Pion0
double Decay::nPI0()
{
	if (fnPI0 < 0 || IsChanged())
		fnPI0 = NeutrinoPseudoMeson(M_Pion, Const::fDPion02);	//check

	return mapEnhance["nPI0"] * fnPI0 * (Ue*Ue + Um*Um + Ut*Ut);
}

//M_Sterile > M_Pion
double Decay::EPI()
{
	if (fEPI < 0 || IsChanged())
		fEPI = LeptonPseudoMeson(M_Electron, M_Pion, Const::fU_ud, Const::fDPion2);
	
	return mapEnhance["EPI"] * fEPI * Ue*Ue;
}

//M_Sterile > M_Pion + M_Muon
double Decay::MUPI()
{
	if (fMUPI < 0 || IsChanged())
		fMUPI = LeptonPseudoMeson(M_Muon, M_Pion, Const::fU_ud, Const::fDPion2);
	
	return mapEnhance["MUPI"] * fMUPI * Um*Um;
}

//M_Sterile > M_Tau + M_Pion
double Decay::TPI()
{
	if (fTPI < 0 || IsChanged())
		fTPI = LeptonPseudoMeson(M_Tau, M_Pion, Const::fU_ud, Const::fDPion2);
	
	return mapEnhance["TPI"] * fTPI * Ut*Ut;
}

//M_Sterile > M_Kaon + M_Electron
double Decay::EKA()
{
	if (fEKA < 0 || IsChanged())
		fEKA = LeptonPseudoMeson(M_Electron, M_Kaon, Const::fU_us, Const::fDKaon2);

	return mapEnhance["EKA"] * fEKA * Ue*Ue;
}

//M_Sterile > M_Kaon + M_Muon
double Decay::MUKA()
{
	if (fMUKA < 0 || IsChanged())
		fMUKA = LeptonPseudoMeson(M_Muon, M_Kaon, Const::fU_us, Const::fDKaon2);

	return mapEnhance["MUKA"] * fMUKA * Um*Um;
}

//M_Sterile > M_Rho
double Decay::nRHO0()
{
	if (fnRHO0 < 0 || IsChanged())
		fnRHO0 = NeutrinoVectorMeson(M_Rho0, Const::fDRho02, Const::fVLight);	//check

	return mapEnhance["nRHO0"] * fnRHO0 * (Ue*Ue + Um*Um + Ut*Ut);
}

//M_Sterile > M_Rho + M_Electron 
double Decay::ERHO()
{
	if (fERHO < 0 || IsChanged())
		fERHO = LeptonVectorMeson(M_Electron, M_Rho, Const::fU_ud, Const::fDRho2, Const::fVLight);	//check

	return mapEnhance["ERHO"] * fERHO * Ue*Ue;
}

//M_Sterile > M_Rho + M_Muon 
double Decay::MURHO()
{
	if (fMURHO < 0 || IsChanged())
		fMURHO = LeptonVectorMeson(M_Muon, M_Rho, Const::fU_ud, Const::fDRho2, Const::fVLight);	//check

	return mapEnhance["MURHO"] * fMURHO * Um*Um;
}

//M_Sterile > M_Kaon* + M_Electron 
double Decay::EKAx()
{
	if (fEKAx < 0 || IsChanged())
		fEKAx = LeptonVectorMeson(M_Electron, M_Kaonx, Const::fU_us, Const::fDKaonx2, Const::fVStrange);	//check

	return mapEnhance["EKAx"] * fEKAx * Ue*Ue;
}

//M_Sterile > M_Kaon0* 
double Decay::nKA0x()
{
	if (fnKA0x < 0 || IsChanged())
		fnKA0x = NeutrinoVectorMeson(M_Kaon0x, Const::fDKaon0x2, Const::fVStrange);	//check

	return mapEnhance["nKA0x"] * fnKA0x * (Ue*Ue + Um*Um + Ut*Ut);
}

//M_Sterile > M_Kaon* + M_Muon 
double Decay::MUKAx()
{
	if (fMUKAx < 0 || IsChanged())
		fMUKAx = LeptonVectorMeson(M_Muon, M_Kaonx, Const::fU_us, Const::fDKaonx2, Const::fVStrange);	//check

	return mapEnhance["MUKAx"] * fMUKAx * Um*Um;
}

//M_Sterile > M_Eta
double Decay::nETA()
{
	if (fnETA < 0 || IsChanged())
		fnETA = NeutrinoPseudoMeson(M_Eta, Const::fDEta2);	//check

	return mapEnhance["nETA"] * fnETA * (Ue*Ue + Um*Um + Ut*Ut);
}

//M_Sterile > M_Eta'
double Decay::nETAi()
{
	if (fnETAi < 0 || IsChanged())
		fnETAi = NeutrinoPseudoMeson(M_Etai, Const::fDEtai2);	//check

	return mapEnhance["nETA"] * fnETA * (Ue*Ue + Um*Um + Ut*Ut);
}

//M_Sterile > M_Omega 
double Decay::nOMEGA()
{
	if (fnOMEGA < 0 || IsChanged())
		fnOMEGA = NeutrinoVectorMeson(M_Omega, Const::fDOmega2, Const::fVLight);	//check

	return mapEnhance["nOMEGA"] * fnOMEGA * (Ue*Ue + Um*Um + Ut*Ut);
}

//M_Sterile > M_Phi 
double Decay::nPHI()
{
	if (fnPHI < 0 || IsChanged())
		fnPHI = NeutrinoVectorMeson(M_Phi, Const::fDPhi2, Const::fVStrange);	//check

	return mapEnhance["nPHI"] * fnPHI * (Ue*Ue + Um*Um + Ut*Ut);
}

//M_Sterile > M_Charm + M_Electron
double Decay::ECHARM()
{
	if (fECHARM < 0 || IsChanged())
		fECHARM = LeptonPseudoMeson(M_Electron, M_Charm, Const::fU_cd, Const::fDCharm2);

	return mapEnhance["ECHARM"] * fECHARM * Ue*Ue;
}

/////////////////
//Generic decay//
/////////////////
double Decay::LeptonPseudoMeson(double M_Lepton, double M_Meson, double vCKM, double fDecay2)
{
	if (M_Sterile >= M_Lepton + M_Meson)
	{
		double dML2 = M_Lepton*M_Lepton/M_Sterile/M_Sterile;
		double dMM2 = M_Meson*M_Meson/M_Sterile/M_Sterile;

		return 2.0 * Const::fGF2 * pow(M_Sterile, 3) *
		       pow(vCKM, 2.0) * fDecay2 * Kine::I1_xy(dML2, dMM2) /
		       (16.0 * Const::fPi);
	}
	else return 0.0;
}

double Decay::NeutrinoPseudoMeson(double M_Meson, double fDecay2)
{
	if (M_Sterile >= M_Neutrino + M_Meson)
	{
		double dMM2 = M_Meson*M_Meson/M_Sterile/M_Sterile;

		return  Const::fGF2 * pow(M_Sterile, 3) *
			fDecay2 * pow((1.0-dMM2), 2.0) / 
			(64.0 * Const::fPi);
	}
	else return 0.0;
}

double Decay::LeptonVectorMeson(double M_Lepton, double M_Meson, double vCKM, double fDecay2, double fVector)
{
	if (M_Sterile >= M_Lepton + M_Meson)
	{
		double dML2 = M_Lepton*M_Lepton/M_Sterile/M_Sterile;
		double dMM2 = M_Meson*M_Meson/M_Sterile/M_Sterile;
		fVector *= fVector;

		return 2.0 * Const::fGF2 * pow(M_Sterile, 3) *
		       pow(vCKM, 2.0) * fDecay2 * fVector * Kine::I2_xy(dML2, dMM2) /
		       (16.0 * Const::fPi);
	}
	else return 0.0;
}

double Decay::NeutrinoVectorMeson(double M_Meson, double fDecay2, double fVector)
{
	if (M_Sterile >= M_Neutrino + M_Meson)
	{
		double dMN2 = M_Neutrino*M_Neutrino/M_Sterile/M_Sterile;
		double dMM2 = M_Meson*M_Meson/M_Sterile/M_Sterile;
		fVector *= fVector;

		return  Const::fGF2 * pow(M_Sterile, 3) *
			fDecay2 * fVector * Kine::I3_xy(dMN2, dMM2) / 
			(2.0 * Const::fPi);
	}
	else return 0.0;
}

double Decay::NeutrinoLeptonAA(double &fCC, double &fNC, double M_Lepton)
{
	if (M_Sterile >= M_Neutrino + 2 * M_Lepton)
	{
		double dMm = M_Lepton / M_Sterile;
		double dMn = M_Neutrino / M_Sterile;
		double gL = -0.5 + Const::fSin2W;
		double gR = Const::fSin2W;
		double Int1 = Kine::I1_xyz(dMn, dMm, dMm);
		double Int2 = Kine::I2_xyz(dMn, dMm, dMm);
		double KF_CC = (gL*gR + gR) * Int2 + (gL*gL + gR*gR + (1+2*gL))*Int1;
		double KF_NC = (gL*gR) * Int2 + (gL*gL + gR*gR)*Int1;

		fCC = Const::fGF2 * pow(M_Sterile, 5) * KF_CC / (96.0 * Const::fPi3);
		fNC = Const::fGF2 * pow(M_Sterile, 5) * KF_NC / (96.0 * Const::fPi3);
	}
	else
	{
		fCC = 0.0;
		fNC = 0.0;
	}
	return 0.0;
}

double Decay::NeutrinoLeptonAB(double M_LeptonA, double M_LeptonB)	//it is doubled for the cc conjugate
{
	if (M_Sterile >= M_Neutrino + M_LeptonA + M_LeptonB)
	{
		double dM1 = M_LeptonA / M_Sterile;
		double dM2 = M_LeptonB / M_Sterile;
		double dMn = M_Neutrino / M_Sterile;

		return 2.0 * Const::fGF2 * pow(M_Sterile, 5) * Kine::I1_xyz(dM1, dMn, dM2) /
			(192.0 * Const::fPi3);
	}
	else return 0.0;
}       
/////////////////////////
//end of generic decays//
/////////////////////////


std::vector<std::string> Decay::ListChannels()
{
	std::map<std::string, ChannelName>::iterator it;
	std::vector<std::string> vList;
	for (it = mapChannel.begin(); it != mapChannel.end(); ++it)
		vList.push_back(it->first);
	return vList;
}

std::vector<ChannelName> Decay::ListNameKeys()
{
	std::map<std::string, ChannelName>::iterator it;
	std::vector<ChannelName> vList;
	for (it = mapChannel.begin(); it != mapChannel.end(); ++it)
		vList.push_back(it->second);
	return vList;
}

bool Decay::IsChanged()
{
	bool Ret = (fabs(M_Sterile - M_Sterile_prev) > 1e-9);

	M_Sterile_prev = M_Sterile;

	//Reset decay widths if changed
	if (Ret)
	{
		fTotal    = -1.0;
		fnnn      = -1.0;
                fnGAMMA   = -1.0;
                fnEE_e    = -1.0;
                fnEE_mt   = -1.0;
                fnEMU     = -1.0;
                fnMUE     = -1.0;
                fnMUMU_m  = -1.0;
                fnMUMU_et = -1.0;
                fnET      = -1.0;
                fnTE      = -1.0;
                fnMUT     = -1.0;
                fnTMU     = -1.0;
                fnPI0     = -1.0;
                fEPI      = -1.0;
                fMUPI     = -1.0;
                fTPI      = -1.0;
                fEKA      = -1.0;
                fMUKA     = -1.0;
                fnRHO0    = -1.0;
                fERHO     = -1.0;
                fMURHO    = -1.0;
                fEKAx     = -1.0;
                fnKA0x    = -1.0;
                fMUKAx    = -1.0;
                fnETA     = -1.0;
                fnETAi    = -1.0;
                fnOMEGA   = -1.0;
                fnPHI     = -1.0;
                fECHARM   = -1.0;
	}

	return Ret;
}

//Get functions
TLorentzVector *Decay::GetNvec()
{
	return N_vec;
}

TLorentzVector Decay::GetDecayProduct(int i, int &ID)
//Particle *Decay::GetDecayProduct(int i, int &ID)
{
	ID = PdgCode[i];
	TLorentzVector Daughter = *(Event->GetDecay(i));
	Daughter.Boost(N_vec->BoostVector());

	return Daughter;
}

double Decay::GetMass()
{
	return M_Sterile;
}

double Decay::GetUe()
{
	return Ue;
}

double Decay::GetUm()
{
	return Um;
}

double Decay::GetUt()
{
	return Ut;
}

//Set functions
void Decay::SetNvec(TLorentzVector &X)
{
	*N_vec = X;
	N_rest->SetPxPyPzE(0, 0, 0, N_vec->M());
}

void Decay::SetMass(double X)
{
	M_Sterile = X;
	TheSpace->SetSterileMass(X);
}

void Decay::SetUe(double X)
{
	Ue = X;
	TheSpace->SetUe(X);
}

void Decay::SetUm(double X)
{
	Um = X;
	TheSpace->SetUm(X);
}

void Decay::SetUt(double X)
{
	Ut = X;
	TheSpace->SetUt(X);
}
