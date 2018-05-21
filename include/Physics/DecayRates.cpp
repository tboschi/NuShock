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
        M_Etai(Const::fMEtai),
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

	Event = new TGenPhaseSpace;
	N_vec = new TLorentzVector;
	N_rest = new TLorentzVector;
}

//return the decay width (Gamma)
//
double Decay::Gamma(Channel Name, double B)
{
	double Result;
	switch(ChanName)
	{
		case Channel::_ALL:
			Result = Total();
			break;
		case Channel::_nnn:
			Result = nnn();
			break;
		case Channel::_nGAMMA:
			Result = nGAMMA();
			break;
		case Channel::_nEE:
			Result = nEE();
			break;
		case Channel::_nEMU:
			Result = nEMU();
			break;
		case Channel::_nMUE:
			Result = nMUE();
			break;
		case Channel::_nMUMU:
			Result = nMUMU();
			break;
		case Channel::_nET:
			Result = nET();
			break;
		case Channel::_nTE:
			Result = nTE();
			break;
		case Channel::_nMUT:
			Result = nMUT();
			break;
		case Channel::_nTMU:
			Result = nTMU();
			break;
		case Channel::_nPI0:
			Result = nPI0();
			break;
		case Channel::_EPI:
			Result = EPI();
			break;
		case Channel::_MUPI:
			Result = MUPI();
			break;
		case Channel::_TPI:
			Result = TPI();
			break;
		case Channel::_EKA:
			Result = EKA();
			break;
		case Channel::_MUKA:
			Result = MUKA();
			break;
		case Channel::_EKAx:
			Result = EKAx();
			break;
		case Channel::_MUKAx:
			Result = MUKAx();
			break;
		case Channel::_nRHO0:
			Result = nRHO0();
			break;
		case Channel::_ERHO:
			Result = ERHO();
			break;
		case Channel::_MURHO:
			Result = MURHO();
			break;
		case Channel::_nETA:
			Result = nETA();
			break;
		case Channel::_nETAi:
			Result = nETAi();
			break;
		case Channel::_nOMEGA:
			Result = nOMEGA();
			break;
		case Channel::_nPHI:
			Result = nPHI();
			break;
		case Channel::_ECHARM:
			Result = ECHARM();
			break;
		case Channel::_ExpALL:
			Result = ExpALL();
			break;
		default:
			Result = -1.0;
			break;
	}

	return Result;
}

//Return Gamma_tot - Gamma of interest
//
double Decay::Other(Channel Name)
{
	if (Gamma(Name) < 0.0)
		return -1.0;
	else return Gamma(_ALL) - Gamma(Name);
}

//Return the branching ration
//
double Decay::Branch(Channel Name)
{
	if (Gamma(Name) < 0.0)
		return -1.0;
	else return Gamma(Name)/Gamma(Channel::_ALL);
}

//Compute phase space for 3 body decay channels
int Decay::PhaseSpace(Channel Name, double &Weight)	//Return number of products 
{								//0 if decay not valid
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
		case Channel::_nGAMMA:
			Mass[0] = M_Photon;
			Mass[1] = M_Neutrino;
			PdgCode[0] = 22;
			PdgCode[1] = 12;
			Products = 2 * Event->SetDecay(*N_rest, 2, Mass);
			Weight = Event->Generate();
			break;
		
		case Channel::_nEE:
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

		case Channel::_nEMU:		//whata about n mu e?
		case Channel::_nMUE:		//whata about n mu e?
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
	
		case Channel::_nPI0:
			Mass[0] = M_Pion0;
			Mass[1] = M_Neutrino;		//Invisible particle
			PdgCode[0] = 111;
			PdgCode[1] = 12;
			Products = 2 * Event->SetDecay(*N_rest, 2, Mass);
			Weight = Event->Generate();
			break;

		case Channel::_EPI:
			Mass[0] = M_Electron;
			Mass[1] = M_Pion;
			PdgCode[0] = 11;
			PdgCode[1] = 211;
			Products = 2 * Event->SetDecay(*N_rest, 2, Mass);
			Weight = Event->Generate();
			break;

		case Channel::_nMUMU:
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

		case Channel::_MUPI:
			Mass[0] = M_Muon;
			Mass[1] = M_Pion;
			PdgCode[0] = 13;
			PdgCode[1] = 211;
			Products = 2 * Event->SetDecay(*N_rest, 2, Mass);
			Weight = Event->Generate();
			break;

		case Channel::_EKA:
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

//Check if some decay is allowed (is mass threshold)
//
bool Decay::IsAllowed(Channel Name)
{
	double Limit = 0.0;

	switch(Name)
	{
		case Channel::_ALL:
		case Channel::_nnn:
			Limit = 3*M_Neutrino;
			break;
		case Channel::_nGAMMA:
			Limit = 3*M_Neutrino;
			break;
		case Channel::_nEE:
		case Channel::_ExpALL:
			Limit = M_Neutrino + 2*M_Electron;
			break;
		case Channel::_nEMU:
		case Channel::_nMUE:
			Limit = M_Neutrino + M_Electron + M_Muon;
			break;
		case Channel::_nMUMU:
			Limit = M_Neutrino + 2*M_Muon;
			break;
		case Channel::_nET:
		case Channel::_nTE:
			Limit = M_Neutrino + M_Electron + M_Tau;
			break;
		case Channel::_nMUT:
		case Channel::_nTMU:
			Limit = M_Neutrino + M_Muon + M_Tau;
			break;
		case Channel::_nPI0:
			Limit = M_Neutrino + M_Pion0;
			break;
		case Channel::_EPI:
			Limit = M_Electron + M_Pion;
			break;
		case Channel::_MUPI:
			Limit = M_Muon + M_Pion;
			break;
		case Channel::_TPI:
			Limit = M_Tau + M_Pion;
			break;
		case Channel::_EKA:
			Limit = M_Electron + M_Kaon;
			break;
		case Channel::_MUKA:
			Limit = M_Muon + M_Kaon;
			break;
		case Channel::_EKAx:
			Limit = M_Electron + M_Kaonx;
			break;
		case Channel::_MUKAx:
			Limit = M_Muon + M_Kaonx;
			break;
		case Channel::_nRHO0:
			Limit = M_Neutrino + M_Rho0;
			break;
		case Channel::_ERHO:
			Limit = M_Electron + M_Rho;
			break;
		case Channel::_MURHO:
			Limit = M_Muon + M_Rho;
			break;
		case Channel::_nETA:
			Limit = M_Neutrino + M_Eta;
			break;
		case Channel::_nETAi:
			Limit = M_Neutrino + M_Etai;
			break;
		case Channel::_nOMEGA:
			Limit = M_Neutrino + M_Omega;
			break;
		case Channel::_nPHI:
			Limit = M_Neutrino + M_Phi;
			break;
		case Channel::_ECHARM:
			Limit = M_Electron + M_Charm;
			break;
		default:
			Limit = 2*M_Sterile;
			break;
	}

	return (M_Sterile >= Limit);
}

//total decay width
double Decay::Total(double A)
{
	return (nnn() + nGAMMA() +
		nEE() + nEMU() + nMUE() + nMUMU() + nET() + nTE() + nMUT() + nTMU() +
		nPI0() + EPI() + + MUPI() + TPI() +
		EKA() + MUKA() + 
		nRHO0() + ERHO() + MURHO() +
		EKAx() + MUKAx() + 
		nETA() + nETAi() + nOMEGA() + nPHI() +
		ECHARM() );
}

//special here
double Decay::ExpALL()
{
	return (nEE() + nMUE() + nMUMU() +
		EPI() + MUPI() +
		EKA() + MUKA() +
		ERHO() + MURHO() );
}

//individual decay channels
//all mixing factors are factorised out
double Decay::nnn()
{
	if (fnnn < 0 || IsChanged())
	{
		if (IsAllowed(Channel::_nnn))
		{
			fnnn = Const::fGF2 * pow(M_Sterile, 5) /
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
		if (IsAllowed(Channel::_nGAMMA))
		{
			double AemPi = Const::fAem / Const::fPi;
			fnGAMMA = Const::fGF2 * pow(M_Sterile, 5) *
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

	return fnEE_e * Ue*Ue + fnEE_mt * (Um*Um + Ut*Ut);
}

//M_Sterile > M_Muon + M_Electron
double Decay::nEMU()	//Antiparticle is Elec
{
	if (fnEMU < 0 || IsChanged())
		fnEMU = NeutrinoLeptonAB(M_Electron, M_Muon);

	return fnEMU * Ue*Ue;
}

double Decay::nMUE()	//Anti is Muon
{
	if (fnMUE < 0 || IsChanged())
		fnMUE = NeutrinoLeptonAB(M_Muon, M_Electron);

	return fnMUE * Um*Um;
}

//M_Sterile > 2 M_Muon
double Decay::nMUMU()
{
	if (fnMUMU_m < 0 || fnMUMU_et < 0 || IsChanged())
		NeutrinoLeptonAA(fnMUMU_m, fnMUMU_et, M_Muon);

	return (fnMUMU_m * Um*Um + fnMUMU_et * (Ue*Ue + Ut*Ut));
}

//M_Sterile > M_Tau + M_Electron
double Decay::nET()	//Antiparticle is Elec
{
	if (fnET < 0 || IsChanged())
		fnET = NeutrinoLeptonAB(M_Electron, M_Tau);

	return fnET * Ue*Ue;
}

double Decay::nTE()	//Anti is Tau
{
	if (fnTE < 0 || IsChanged())
		fnTE = NeutrinoLeptonAB(M_Tau, M_Electron);

	return fnTE * Ut*Ut;
}

//M_Sterile > M_Tau + M_Muon
double Decay::nMUT()	//Antiparticle is Muon
{
	if (fnMUT < 0 || IsChanged())
		fnMUT = NeutrinoLeptonAB(M_Muon, M_Tau);

	return fnMUT * Um*Um;
}

double Decay::nTMU()	//Anti is Tau
{
	if (fnTMU < 0 || IsChanged())
		fnTMU = NeutrinoLeptonAB(M_Tau, M_Muon);

	return fnTMU * Ut*Ut;
}

//M_Sterile > M_Pion0
double Decay::nPI0()
{
	if (fnPI0 < 0 || IsChanged())
		fnPI0 = NeutrinoPseudoMeson(M_Pion, Const::fDPion02);	//check

	return fnPI0 * (Ue*Ue + Um*Um + Ut*Ut);
}

//M_Sterile > M_Pion
double Decay::EPI()
{
	if (fEPI < 0 || IsChanged())
		fEPI = LeptonPseudoMeson(M_Electron, M_Pion, Const::fU_ud, Const::fDPion2);
	
	return fEPI * Ue*Ue;
}

//M_Sterile > M_Pion + M_Muon
double Decay::MUPI()
{
	if (fMUPI < 0 || IsChanged())
		fMUPI = LeptonPseudoMeson(M_Muon, M_Pion, Const::fU_ud, Const::fDPion2);
	
	return fMUPI * Um*Um;
}

//M_Sterile > M_Tau + M_Pion
double Decay::TPI()
{
	if (fTPI < 0 || IsChanged())
		fTPI = LeptonPseudoMeson(M_Tau, M_Pion, Const::fU_ud, Const::fDPion2);
	
	return fTPI * Ut*Ut;
}

//M_Sterile > M_Kaon + M_Electron
double Decay::EKA()
{
	if (fEKA < 0 || IsChanged())
		fEKA = LeptonPseudoMeson(M_Electron, M_Kaon, Const::fU_us, Const::fDKaon2);

	return fEKA * Ue*Ue;
}

//M_Sterile > M_Kaon + M_Muon
double Decay::MUKA()
{
	if (fMUKA < 0 || IsChanged())
		fMUKA = LeptonPseudoMeson(M_Muon, M_Kaon, Const::fU_us, Const::fDKaon2);

	return fMUKA * Um*Um;
}

//M_Sterile > M_Rho
double Decay::nRHO0()
{
	if (fnRHO0 < 0 || IsChanged())
		fnRHO0 = NeutrinoVectorMeson(M_Rho0, Const::fDRho02, Const::fVLight);	//check

	return fnRHO0 * (Ue*Ue + Um*Um + Ut*Ut);
}

//M_Sterile > M_Rho + M_Electron 
double Decay::ERHO()
{
	if (fERHO < 0 || IsChanged())
		fERHO = LeptonVectorMeson(M_Electron, M_Rho, Const::fU_ud, Const::fDRho2, Const::fVLight);	//check

	return fERHO * Ue*Ue;
}

//M_Sterile > M_Rho + M_Muon 
double Decay::MURHO()
{
	if (fMURHO < 0 || IsChanged())
		fMURHO = LeptonVectorMeson(M_Muon, M_Rho, Const::fU_ud, Const::fDRho2, Const::fVLight);	//check

	return fMURHO * Um*Um;
}

//M_Sterile > M_Kaon* + M_Electron 
double Decay::EKAx()
{
	if (fEKAx < 0 || IsChanged())
		fEKAx = LeptonVectorMeson(M_Electron, M_Kaonx, Const::fU_us, Const::fDKaonx2, Const::fVStrange);	//check

	return fEKAx * Ue*Ue;
}

//M_Sterile > M_Kaon* + M_Muon 
double Decay::MUKAx()
{
	if (fMUKAx < 0 || IsChanged())
		fMUKAx = LeptonVectorMeson(M_Muon, M_Kaonx, Const::fU_us, Const::fDKaonx2, Const::fVStrange);	//check

	return fMUKAx * Um*Um;
}

//M_Sterile > M_Eta
double Decay::nETA()
{
	if (fnETA < 0 || IsChanged())
		fnETA = NeutrinoPseudoMeson(M_Eta, Const::fDEta2);	//check

	return fnETA * (Ue*Ue + Um*Um + Ut*Ut);
}

//M_Sterile > M_Eta'
double Decay::nETAi()
{
	if (fnETAi < 0 || IsChanged())
		fnETAi = NeutrinoPseudoMeson(M_Etai, Const::fDEtai2);	//check

	return fnETAi * (Ue*Ue + Um*Um + Ut*Ut);
}

//M_Sterile > M_Omega 
double Decay::nOMEGA()
{
	if (fnOMEGA < 0 || IsChanged())
		fnOMEGA = NeutrinoVectorMeson(M_Omega, Const::fDOmega2, Const::fVLight);	//check

	return fnOMEGA * (Ue*Ue + Um*Um + Ut*Ut);
}

//M_Sterile > M_Phi 
double Decay::nPHI()
{
	if (fnPHI < 0 || IsChanged())
		fnPHI = NeutrinoVectorMeson(M_Phi, Const::fDPhi2, Const::fVStrange);	//check

	return fnPHI * (Ue*Ue + Um*Um + Ut*Ut);
}

//M_Sterile > M_Charm + M_Electron
double Decay::ECHARM()
{
	if (fECHARM < 0 || IsChanged())
		fECHARM = LeptonPseudoMeson(M_Electron, M_Charm, Const::fU_cd, Const::fDCharm2);

	return fECHARM * Ue*Ue;
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
		       pow(vCKM, 2.0) * fDecay2 * I1_xy(dML2, dMM2) /
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
		       pow(vCKM, 2.0) * fDecay2 * fVector * I2_xy(dML2, dMM2) /
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
			fDecay2 * fVector * I3_xy(dMN2, dMM2) / 
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
		double Int1 = I1_xyz(dMn, dMm, dMm);
		double Int2 = I2_xyz(dMn, dMm, dMm);
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
		double dMA = M_LeptonA / M_Sterile;
		double dMB = M_LeptonB / M_Sterile;
		double dMn = M_Neutrino / M_Sterile;

		return 2.0 * Const::fGF2 * pow(M_Sterile, 5) * I1_xyz(dMA, dMn, dMB) /
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

bool Decay::IsChanged()
{
	bool Ret = (fabs(M_Sterile - M_Sterile_prev) > 1e-9);

	M_Sterile_prev = M_Sterile;

	//Reset decay widths if changed
	if (Ret)
	{
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

////////////////////////////
//function for decay rates//
////////////////////////////
//
double Decay::I1_f(double t)	//To be integrated, using global variables fX, fY, fZ
{
	return (t-fX*fX-fY*fY)*(1+fZ*fZ-t)*sqrt(Kine::Lambda(t, fX*fX, fY*fY)*Kine::Lambda(1, t, fZ*fZ)) / t;
}

double Decay::I1_xyz(double X, double Y, double Z)
{
	fFunc = &I2_f(double);
	fA = (1-Z)*(1-Z);
	fB = (X+Y)*(X+Y);
	if (fA > fB)
	{
		double tmp = fB;
		fB = fA;
		fA = tmp;
	}

	fX = X;
	fY = Y;
	fZ = Z;

	return 12.0*(fB-fA)*Inte::BooleIntegration(this);

	double a, b;
	double h = (B-A)/1000;
	double Integral = 0;	//Boole's method for integration
	for (a = A; b < B; a = b)
	{
		b = a + h;
		Integral += h/90.0 * (7*I1_f(a) + 
				      32*I1_f((3*a+b)/4.0) + 
				      12*I1_f((a+b)/2.0) + 
				      32*I1_f((a+3*b)/4.0) +
				      7*I1_f(b));
	}	
	return 12.0 * Integral;
}

double Decay::I1_xy(double X, double Y)
{
	return ((1+X-Y)*(1+X) - 4*X) * sqrt(Lambda(1.0, X, Y));
}

double Decay::I2_f(double t)	//To be integrated
{
	return (1+fX*fX-t)*sqrt(Lambda(t, fY*fY, fZ*fZ))*sqrt(Lambda(1.0, t, fX*fX)) / t;
}

double Decay::I2_xyz(double X, double Y, double Z)
{
	fFunc = &I2_f(double);
	fA = (1-X)*(1-X);
	fB = (Y+Z)*(Y+Z);
	if (fA > fB)
	{
		double tmp = fB;
		fB = fA;
		fA = tmp;
	}

	fX = X;
	fY = Y;
	fZ = Z;

	return 24.0*fY*fZ* (fB-fA)*Inte::BooleIntegration(this);

	///this remove
	double a, b;
	double h = (B-A)/1000;
	double Integral = 0.0;	//Boole's method for integration
	for (a = A; b < B; a = b)
	{
		b = a + h;
		Integral += h/90.0 * (7*I2_f(a, X, Y, Z) + 
				      32*I2_f((3*a+b)/4.0, X, Y, Z) + 
				      12*I2_f((a+b)/2.0, X, Y, Z) + 
				      32*I2_f((a+3*b)/4.0, X, Y, Z) +
				      7*I2_f(b, X, Y, Z));
	}	
	return 24.0*Y*Z * Integral;
}

double Decay::I2_xy(double X, double Y)
{
	return ((1+X-Y)*(1+X+2*Y) - 4*X) * sqrt(Lambda(1.0, X, Y));
}

double Decay::I3_xy(double X, double Y)
{
	return (1+2*Y)*(1-Y) * sqrt(Lambda(1.0, X, Y));
}

double Decay::Variable(double t)
{
	return (*fFunc)(fA + (fB-fA)*t);
}

