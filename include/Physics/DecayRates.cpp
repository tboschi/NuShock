#include "DecayRates.h"

DecayRates::DecayRates(double Mass, double Ue, double Um, double Ut)	: //Decay rates calculator
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
        M_Etai(Const::fMEtai),
        M_Phi(Const::fMPhi),
        M_Tau(Const::fMTau),
        M_Charm(Const::fMCharm)
{
	M_Sterile_prev = -1.0;

	SetMass(Mass);
	SetUe(Ue);
	SetUm(Um);
	SetUt(Ut);

	Event = new TGenPhaseSpace;
	N_vec = new TLorentzVector;
	N_rest = new TLorentzVector;
}

//return the decay width (Gamma)
//
double DecayRates::Gamma(Channel Name)
{
	double Result = 0.0;

	if (GetHelicity() == 0)
	{
		SetHelicity(1);
		Result += Gamma(Name);
		
		SetHelicity(1);
		Result += Gamma(Name);

		SetHelicity(1);
		return Result / 2.0;
	}
	else
	{
		switch(Name)
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
			case _ExpALL:
				Result = ExpALL();
				break;
			default:
				Result = 0.0;
				break;
		}
		return Result;
	}
}

//Return Gamma_tot - Gamma of interest
//
double DecayRates::Other(Channel Name)
{
	if (Gamma(Name) < 0.0)
		return -1.0;
	else return Gamma(_ALL) - Gamma(Name);
}

//Return the branching ration
//
double DecayRates::Branch(Channel Name)
{
	if (Gamma(Name) < 0.0)
		return -1.0;
	else return Gamma(Name)/Gamma(_ALL);
}

//Compute phase space for 3 body decay channels
int DecayRates::PhaseSpace(Channel Name, double &Weight)	//Return number of products 
{								//0 if decay not valid
	TheSpace->SetParent(Channel);

	double Mass[3];
	int Products;

	switch(Name)
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

//Check if some decay is allowed (is mass threshold)
//
bool DecayRates::IsAllowed(Channel Name)
{
	double Limit = 0.0;

	switch(Name)
	{
		case _ALL:
		case _nnn:
			Limit = 3*M_Neutrino;
			break;
		case _nGAMMA:
			Limit = 3*M_Neutrino;
			break;
		case _nEE:
		case _ExpALL:
			Limit = M_Neutrino + 2*M_Electron;
			break;
		case _nEMU:
		case _nMUE:
			Limit = M_Neutrino + M_Electron + M_Muon;
			break;
		case _nMUMU:
			Limit = M_Neutrino + 2*M_Muon;
			break;
		case _nET:
		case _nTE:
			Limit = M_Neutrino + M_Electron + M_Tau;
			break;
		case _nMUT:
		case _nTMU:
			Limit = M_Neutrino + M_Muon + M_Tau;
			break;
		case _nPI0:
			Limit = M_Neutrino + M_Pion0;
			break;
		case _EPI:
			Limit = M_Electron + M_Pion;
			break;
		case _MUPI:
			Limit = M_Muon + M_Pion;
			break;
		case _TPI:
			Limit = M_Tau + M_Pion;
			break;
		case _EKA:
			Limit = M_Electron + M_Kaon;
			break;
		case _MUKA:
			Limit = M_Muon + M_Kaon;
			break;
		case _EKAx:
			Limit = M_Electron + M_Kaonx;
			break;
		case _MUKAx:
			Limit = M_Muon + M_Kaonx;
			break;
		case _nRHO0:
			Limit = M_Neutrino + M_Rho0;
			break;
		case _ERHO:
			Limit = M_Electron + M_Rho;
			break;
		case _MURHO:
			Limit = M_Muon + M_Rho;
			break;
		case _nETA:
			Limit = M_Neutrino + M_Eta;
			break;
		case _nETAi:
			Limit = M_Neutrino + M_Etai;
			break;
		case _nOMEGA:
			Limit = M_Neutrino + M_Omega;
			break;
		case _nPHI:
			Limit = M_Neutrino + M_Phi;
			break;
		case _ECHARM:
			Limit = M_Electron + M_Charm;
			break;
		default:
			Limit = 2*GetMass();
			break;
	}

	return (GetMass() >= Limit);
}

//total decay width
double DecayRates::Total()
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
double DecayRates::ExpALL()
{
	return (nEE() + nMUE() + nMUMU() +
		EPI() + MUPI() +
		EKA() + MUKA() +
		ERHO() + MURHO() );
}

//individual decay channels
//all mixing factors are factorised out
double DecayRates::nnn()
{
	if (IsAllowed() && (fnnn < 0 || IsChanged()))
	{
		if (IsAllowed(_nnn))
			fnnn = Const::fGF2 * pow(GetMass(), 5) /
			       (96.0 * Const::fPi3);
		else fnnn = 0.0;
	}

	return fnnn * (Ue*Ue + Um*Um + Ut*Ut);
}

double DecayRates::nGAMMA()
{
	if (fnGAMMA < 0 || IsChanged())
	{
		if (IsAllowed(_nGAMMA))
		{
			double AemPi = Const::fAem / Const::fPi;
			fnGAMMA = Const::fGF2 * pow(GetMass(), 5) *
			       (27.0/32.0 * AemPi) / (192.0 * Const::fPi3);
		}
		else fnGAMMA = 0.0;
	}

	return fnGAMMA * (Ue*Ue + Um*Um + Ut*Ut);
}

//M_Sterile > 2 M_Electron (always)
double DecayRates::nEE()
{
	if (fnEE_e < 0 || fnEE_mt < 0 || IsChanged())
	{
		if (IsAllowed(_nEE))
			NeutrinoLeptonAA(fnEE_e, fnEE_mt, M_Electron);
		else
		{
			fnEE_e  = 0.0;
			fnEE_mt = 0.0;
		}
	}

	return fnEE_e * Ue*Ue + fnEE_mt * (Um*Um + Ut*Ut);
}

//M_Sterile > M_Muon + M_Electron
double DecayRates::nEMU()	//Antiparticle is Elec
{
	if (fnEMU < 0 || IsChanged())
	{
		if (IsAllowed(_nEMU))
			fnEMU = NeutrinoLeptonAB(M_Electron, M_Muon);
		else
			fnEMU = 0.0;
	}

	return fnEMU * Ue*Ue;
}

double DecayRates::nMUE()	//Anti is Muon
{
	if (fnMUE < 0 || IsChanged())
	{
		if (IsAllowed(_nMUE))
			fnMUE = NeutrinoLeptonAB(M_Muon, M_Electron);
		else
			fnMUE = 0.0;
	}

	return fnMUE * Um*Um;
}

//M_Sterile > 2 M_Muon
double DecayRates::nMUMU()
{
	if (fnMUMU_m < 0 || fnMUMU_et < 0 || IsChanged())
	{
		if (IsAllowed(_nMUMU)) 
			NeutrinoLeptonAA(fnMUMU_m, fnMUMU_et, M_Muon);
		else
		{
			fnMUMU_m  = 0.0;
			fnMUMU_et = 0.0;
		}
	}

	return (fnMUMU_m * Um*Um + fnMUMU_et * (Ue*Ue + Ut*Ut));
}

//M_Sterile > M_Tau + M_Electron
double DecayRates::nET()	//Antiparticle is Elec
{
	if (fnET < 0 || IsChanged())
	{
		if (IsAllowed(_nET))
			fnET = NeutrinoLeptonAB(M_Electron, M_Tau);
		else
			fnET = 0.0;
	}

	return fnET * Ue*Ue;
}

double DecayRates::nTE()	//Anti is Tau
{
	if (fnTE < 0 || IsChanged())
	{
		if (IsAllowed(_nTE))
			fnTE = NeutrinoLeptonAB(M_Tau, M_Electron);
		else
			fnTE = 0.0;
	}

	return fnTE * Ut*Ut;
}

//M_Sterile > M_Tau + M_Muon
double DecayRates::nMUT()	//Antiparticle is Muon
{
	if (fnMUT < 0 || IsChanged())
	{
		if (IsAllowed(_nMUT))
			fnMUT = NeutrinoLeptonAB(M_Muon, M_Tau);
		else
			fnMUT = 0.0;
	}

	return fnMUT * Um*Um;
}

double DecayRates::nTMU()	//Anti is Tau
{
	if (fnTMU < 0 || IsChanged())
	{
		if (IsAllowed(_nTMU))
			fnTMU = NeutrinoLeptonAB(M_Tau, M_Muon);
		else
			fnTMU = 0.0;
	}

	return fnTMU * Ut*Ut;
}

//M_Sterile > M_Pion0
double DecayRates::nPI0()
{
	if (fnPI0 < 0 || IsChanged())
	{
		if (IsAllowed(_nPI0))
			fnPI0 = NeutrinoPseudoMeson(M_Pion, Const::fDPion02);	//check
		else
			fnPI0 = 0.0;
	}

	return fnPI0 * (Ue*Ue + Um*Um + Ut*Ut);
}

//M_Sterile > M_Pion
double DecayRates::EPI()
{
	if (fEPI < 0 || IsChanged())
	{
		if (IsAllowed(_EPI))
			fEPI = LeptonPseudoMeson(M_Electron, M_Pion, Const::fU_ud, Const::fDPion2);
		else
			fEPI = 0.0;
	}
	
	return fEPI * Ue*Ue;
}

//M_Sterile > M_Pion + M_Muon
double DecayRates::MUPI()
{
	if (fMUPI < 0 || IsChanged())
	{
		if (IsAllowed(_MUPI))
			fMUPI = LeptonPseudoMeson(M_Muon, M_Pion, Const::fU_ud, Const::fDPion2);
		else
			fMUPI = 0.0;
	}
	
	return fMUPI * Um*Um;
}

//M_Sterile > M_Tau + M_Pion
double DecayRates::TPI()
{
	if (fTPI < 0 || IsChanged())
	{
		if (IsAllowed(_TPI))
			fTPI = LeptonPseudoMeson(M_Tau, M_Pion, Const::fU_ud, Const::fDPion2);
		else
			fTPI = 0.0;
	}
	
	return fTPI * Ut*Ut;
}

//M_Sterile > M_Kaon + M_Electron
double DecayRates::EKA()
{
	if (fEKA < 0 || IsChanged())
	{
		if (IsAllowed(_EKA))
			fEKA = LeptonPseudoMeson(M_Electron, M_Kaon, Const::fU_us, Const::fDKaon2);
		else
			fEKA = 0.0;
	}

	return fEKA * Ue*Ue;
}

//M_Sterile > M_Kaon + M_Muon
double DecayRates::MUKA()
{
	if (fMUKA < 0 || IsChanged())
	{
		if (IsAllowed(_MUKA))
			fMUKA = LeptonPseudoMeson(M_Muon, M_Kaon, Const::fU_us, Const::fDKaon2);
		else
			fMUKA = 0.0;
	}

	return fMUKA * Um*Um;
}

//M_Sterile > M_Rho
double DecayRates::nRHO0()
{
	if (fnRHO0 < 0 || IsChanged())
	{
		if (IsAllowed(_nRHO0))
			fnRHO0 = NeutrinoVectorMeson(M_Rho0, Const::fDRho02, Const::fVLight);	//check
		else
			fnRHO0 = 0.0;
	}

	return fnRHO0 * (Ue*Ue + Um*Um + Ut*Ut);
}

//M_Sterile > M_Rho + M_Electron 
double DecayRates::ERHO()
{
	if (fERHO < 0 || IsChanged())
	{
		if (IsAllowed(_ERHO))
			fERHO = LeptonVectorMeson(M_Electron, M_Rho, Const::fU_ud, Const::fDRho2, Const::fVLight);	//check
		else
			fERHO = 0.0;
	}

	return fERHO * Ue*Ue;
}

//M_Sterile > M_Rho + M_Muon 
double DecayRates::MURHO()
{
	if (fMURHO < 0 || IsChanged())
	{
		if (IsAllowed(_MURHO))
			fMURHO = LeptonVectorMeson(M_Muon, M_Rho, Const::fU_ud, Const::fDRho2, Const::fVLight);	//check
		else
			fMURHO = 0.0;
	}

	return fMURHO * Um*Um;
}

//M_Sterile > M_Kaon* + M_Electron 
double DecayRates::EKAx()
{
	if (fEKAx < 0 || IsChanged())
	{
		if (IsAllowed(_EKAx))
			fEKAx = LeptonVectorMeson(M_Electron, M_Kaonx, Const::fU_us, Const::fDKaonx2, Const::fVStrange);	//check
		else
			fEKAx = 0.0;
	}

	return fEKAx * Ue*Ue;
}

//M_Sterile > M_Kaon* + M_Muon 
double DecayRates::MUKAx()
{
	if (fMUKAx < 0 || IsChanged())
	{
		if (IsAllowed(_MUKAx))
			fMUKAx = LeptonVectorMeson(M_Muon, M_Kaonx, Const::fU_us, Const::fDKaonx2, Const::fVStrange);	//check
		else
			fMUKAx = 0.0;
	}

	return fMUKAx * Um*Um;
}

//M_Sterile > M_Eta
double DecayRates::nETA()
{
	if (fnETA < 0 || IsChanged())
	{
		if (IsAllowed(_nETA))
			fnETA = NeutrinoPseudoMeson(M_Eta, Const::fDEta2);	//check
		else
			fnETA = 0.0;
	}

	return fnETA * (Ue*Ue + Um*Um + Ut*Ut);
}

//M_Sterile > M_Eta'
double DecayRates::nETAi()
{
	if (fnETAi < 0 || IsChanged())
	{
		if (IsAllowed(_nETAi))
			fnETAi = NeutrinoPseudoMeson(M_Etai, Const::fDEtai2);	//check
		else
			fnETAi = 0.0;
	}

	return fnETAi * (Ue*Ue + Um*Um + Ut*Ut);
}

//M_Sterile > M_Omega 
double DecayRates::nOMEGA()
{
	if (fnOMEGA < 0 || IsChanged())
	{
		if (IsAllowed(_nOMEGA))
			fnOMEGA = NeutrinoVectorMeson(M_Omega, Const::fDOmega2, Const::fVLight);	//check
		else
			fnOMEGA = 0.0;
	}

	return fnOMEGA * (Ue*Ue + Um*Um + Ut*Ut);
}

//M_Sterile > M_Phi 
double DecayRates::nPHI()
{
	if (fnPHI < 0 || IsChanged())
	{
		if (IsAllowed(_nPHI))
			fnPHI = NeutrinoVectorMeson(M_Phi, Const::fDPhi2, Const::fVStrange);	//check
		else
			fnPHI = 0.0;
	}

	return fnPHI * (Ue*Ue + Um*Um + Ut*Ut);
}

//M_Sterile > M_Charm + M_Electron
double DecayRates::ECHARM()
{
	if (fECHARM < 0 || IsChanged())
	{
		if (IsAllowed(_ECHARM))
			fECHARM = LeptonPseudoMeson(M_Electron, M_Charm, Const::fU_cd, Const::fDCharm2);
		else
			fECHARM = 0.0;
	}

	return fECHARM * Ue*Ue;
}

/////////////////
//Generic decay//
/////////////////
//
//CC version possible
double DecayRates::LeptonPseudoMeson(double M_Lepton, double M_Meson, double vCKM, double fDecay2)
{
	double dML2 = M_Lepton*M_Lepton/GetMass()/GetMass();
	double dMM2 = M_Meson*M_Meson/GetMass()/GetMass();

	return GetMult() * Const::fGF2 * pow(vCKM, 2.0) * fDecay2 / (16.0 * Const::fPi) *
		pow(GetMass(), 3) * I_PseudoMeson(dML2, dMM2);
}

double DecayRates::NeutrinoPseudoMeson(double M_Meson, double fDecay2)
{
	double dMN2 = M_Neutrino*M_Neutrino/GetMass()/GetMass();
	double dMM2 = M_Meson*M_Meson/GetMass()/GetMass();

	return Const::fGF2 * fDecay2 * / (64.0 * Const::fPi);
		pow(GetMass(), 3) * I_PseudoMeson(dMN2, dMM2);
}

double DecayRates::I_PseudoMeson(double x, double y, double theta)
{
	double cos0 = theta < 0 ? 0.0 : cos(theta);
	double fc = theta < 0.0 ? 2.0 : 1.0;

	double Lambda = sqrt(Kine::Kallen(1, x, y)) * cos0;
	if (GetHelicity() < 0)
		return fc * Lambda * (pow(1 - x, 2) - y*(1+x) + Lambda);
	else if (GetHelicity() > 0)
		return fc * Lambda * (pow(1 - x, 2) - y*(1+x) - Lambda);
}

//no helicity version for this one
//CC version possible
double DecayRates::LeptonVectorMeson(double M_Lepton, double M_Meson, double vCKM, double fDecay2, double fVector)
{
	double dML2 = M_Lepton*M_Lepton/GetMass()/GetMass();
	double dMM2 = M_Meson*M_Meson/GetMass()/GetMass();
	double Lambda = sqrt(Kine::Kallen(1, dML2, dMM2));
	fVector *= fVector;

	return GetMult() * Const::fGF2 * pow(vCKM, 2.0) * fDecay2 * fVector / (16.0 * Const::fPi) *
	       pow(GetMass(), 3) * Lambda * (pow(1-dML2, 2) + dMM2 * (1 + dML2 - 2*dMM2));
}

//no helicity version for this one
double DecayRates::NeutrinoVectorMeson(double M_Meson, double fDecay2, double fVector)
{
	double dMN2 = M_Neutrino*M_Neutrino/GetMass()/GetMass();
	double dMM2 = M_Meson*M_Meson/GetMass()/GetMass();
	double Lambda = sqrt(Kine::Kallen(1, dML2, dMM2));
	fVector *= fVector;

	return  Const::fGF2 * fDecay2 * fVector / (2.0 * Const::fPi) * 
		pow(GetMass(), 3) * (1+2*dMM2) * pow(1-dMM2, 2);
}

double DecayRates::NeutrinoLeptonAA(double &fAAA, double &fABB, double M_Lepton)
{
	double dMN2 = M_Neutrino*M_Neutrino/GetMass()/GetMass();
	double dML2 = M_Lepton*M_Lepton/GetMass()/GetMass();

	double gL = -0.5 + Const::fSin2W;
	double gR = Const::fSin2W;

	double IntWW = I_WW(dMN2, dML2, dML2);	//W or Z mediation
	double IntWZ = I_WZ(dMN2, dML2, dML2);	//cross term for W + Z

	double Amp_AAA = (gL*gL + gR*gR + 1 + 2*gL)*IntWW + (gL*gR + gR) * IntWZ;	//both W and Z
	double Amp_ABB = (gL*gL + gR*gR           )*IntWW + (gL*gR     ) * IntWZ;	//only Z

	fAAA = Const::fGF2 * pow(GetMass(), 5) *  / (16.0 * Const::fPi3);	//nu flavour is the same of leptons
	fABB = Const::fGF2 * pow(GetMass(), 5) *  / (16.0 * Const::fPi3);	//nu flavour is different from leptons

	return 0.0;
}

//CC version also available
double DecayRates::NeutrinoLeptonAB(double M_LeptonA, double M_LeptonB, double M_Neutrino)
{
	double dMN2 = M_Neutrino*M_Neutrino/GetMass()/GetMass();
	double dMA2 = M_LeptonA*M_LeptonB/GetMass()/GetMass();
	double dMB2 = M_LeptonB*M_LeptonB/GetMass()/GetMass();

	return GetMult() * Const::fGF2 / (16.0 * Const::fPi3) * 
	 pow(GetMass(), 5) * I_WW(dMN2, dMA2, dMB2);
}

double DecayRates::I_WW(double x, double y, double z, double theta)
{
	double sInf = x*x + y*y + 2*sqrt(x*y);
	double sSup = 1 - 2*sqrt(z) + z;
	double cos0 = theta < 0 ? 0.0 : cos(theta);
	double fc = theta < 0.0 ? 2.0 : 1.0;

	I_var.clear();
	I_var.push_back(sInf);	//0
	I_var.push_back(sSup);	//1

	I_var.push_back(x);	//2
	I_var.push_back(y);	//3
	I_var.push_back(z);	//4
	I_var.push_back(cos0);	//5

	SetIntegrand(&I_WW_s);
	return fc * (sSup - sInf) * Inte::BooleIntegration(this); 
}

double DecayRates::I_WW_s(double s)
{
	const double &sInf = I_var.at(0);
	const double &sSup = I_var.at(0);

	const double &c2 = I_var.at(2);
	const double &b2 = I_var.at(3);
	const double &a2 = I_var.at(4);
	const double &cos0 = I_var.at(5);

	double S = sInf + (sSup - sInf) * s;

	double Lambda0 = sqrt(Kine::Kallen(1, S, a2));
	double Lambda1 = sqrt(Kine::Kallen(S, b2, c2));

	if (GetHelicity() < 0)
		return (S - b2 - c2) * (1 + a2 - S - Lambda0*cos0) * Lambda0 * Lambda1 / S;
	else if (GetHelicity() > 0)
		return (S - b2 - c2) * (1 + a2 - S + Lambda0*cos0) * Lambda0 * Lambda1 / S;
}

double DecayRates::I_WZ(double x, double y, double z, double theta)
{
	double sInf = x*x + y*y + 2*sqrt(x*y);
	double sSup = 1 - 2*sqrt(z) + z;
	double cos0 = theta < 0 ? 0.0 : cos(theta);
	double fc = theta < 0.0 ? 2.0 : 1.0;

	I_var.clear();
	I_var.push_back(sInf);	//0
	I_var.push_back(sSup);	//1

	I_var.push_back(x);	//2
	I_var.push_back(y);	//3
	I_var.push_back(z);	//4
	I_var.push_back(cos0);	//5

	SetIntegrand(&I_WZ_s);
	return fc * (sSup - sInf) * Inte::BooleIntegration(this); 
}

double DecayRates::I_WZ_s(double s)	//2 -> x, 3 -> y, 4 -> z, 0 -> 5
{
	const double &sInf = I_var.at(0);
	const double &sSup = I_var.at(0);

	const double &c2 = I_var.at(2);
	const double &b2 = I_var.at(3);
	const double &a2 = I_var.at(4);
	const double &cos0 = I_var.at(5);

	double S = sInf + (sSup - sInf) * s;

	double Lambda0 = sqrt(Kine::Kallen(1, S, a2));
	double Lambda1 = sqrt(Kine::Kallen(S, b2, c2));
	double Sum = S - a2 + Lambda0*cos0;

	//not sure here!
	if (GetHelicity() < 0)
		return b2 * c2 * (S - a2 - Lambda0 * (1 + cos0) / 2.0) * Lambda0 * Lambda1 / S;
	else if (GetHelicity() > 0)
		return b2 * c2 * (S - a2 + Lambda0 * (1 + cos0) / 2.0) * Lambda0 * Lambda1 / S;
}

//Integration set up
//
void DecayRates::SetIntegrand(double (DecayRates::*FF)(double))
{
	I_integrand = FF;
}

double DecayRate::Integrand(double x)
{
	return (*I_integrand)(x);
}

/////////////////////////
//end of generic decays//
/////////////////////////


/*
std::vector<std::string> DecayRates::ListChannels()
{
	std::map<std::string, ChannelName>::iterator it;
	std::vector<std::string> vList;
	for (it = mapChannel.begin(); it != mapChannel.end(); ++it)
		vList.push_back(it->first);
	return vList;
}
*/

bool DecayRates::IsChanged()
{
	bool Ret = (fabs(GetMass() - M_Sterile_prev) > 1e-9);

	M_Sterile_prev = GetMass();

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
TLorentzVector *DecayRates::GetNvec()
{
	return N_vec;
}

TLorentzVector DecayRates::GetDecayProduct(int i, int &ID)
//Particle *DecayRates::GetDecayProduct(int i, int &ID)
{
	ID = PdgCode[i];
	TLorentzVector Daughter = *(Event->GetDecay(i));
	Daughter.Boost(N_vec->BoostVector());

	return Daughter;
}

double DecayRates::GetMass()
{
	return M_Sterile;
}

double DecayRates::GetUe()
{
	return Ue;
}

double DecayRates::GetUm()
{
	return Um;
}

double DecayRates::GetUt()
{
	return Ut;
}

bool DecayRates::GetFermion()
{
	return bFermion;
}

int DecayRates::GetMult()
{
	return 2-GetFermion();
}

int DecayRates::GetHelicity()
{
	return iHel;
}

//Set functions
void DecayRates::SetNvec(TLorentzVector &X)
{
	*N_vec = X;
	N_rest->SetPxPyPzE(0, 0, 0, N_vec->M());
}

void DecayRates::SetMass(double X)
{
	M_Sterile = X;
	TheSpace->SetSterileMass(X);
}

void DecayRates::SetUe(double X)
{
	Ue = X;
	TheSpace->SetUe(X);
}

void DecayRates::SetUm(double X)
{
	Um = X;
	TheSpace->SetUm(X);
}

void DecayRates::SetUt(double X)
{
	Ut = X;
	TheSpace->SetUt(X);
}

void DecayRates::SetFermion(bool Fermion)
{
	bFermion = Fermion;	//true for Dirac, false for Majorana
}

void DecayRates::SetHelicity(int Helicity)
{
	iHel = Helicity;	//-1 for Left, +1 for Right, 0 for unpolarised
}

void DecayRates::SetParent(double Mass, double* Mixings, bool Fermion, bool Helicity)
{
	SetMass(Mass);
	SetUe(Mixings[0]);
	SetUm(Mixings[1]);
	SetUt(Mixings[2]);
	SetFermion(Fermion);
	SetHelicity(Helicity);
}
