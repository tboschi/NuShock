#include "Physics/DecayRates.h"

DecayRates::DecayRates()
{
	Channel_prev = _undefined;
	Reset();
}

Amplitude::Channel DecayRates::FindChannel(std::string Name)
{
	if (chMap.size() == 0)
		LoadMap();

	std::map<Amplitude::Channel, std::string>::iterator it = chMap.begin();
	for (std::advance(it, 2); it != chMap.end(); ++it)
		if (it->second == Name)
			return it->first;

	return _undefined;
}

std::string DecayRates::FindChannel(Amplitude::Channel Name)
{
	return chMap[Name];
}

std::vector<Amplitude::Channel> DecayRates::ListChannels()
{
	if (chMap.size() == 0)
		LoadMap();

	std::vector<Amplitude::Channel> vName;

	std::map<Amplitude::Channel, std::string>::iterator itS = chMap.begin();
	std::map<Amplitude::Channel, std::string>::iterator itE = chMap.begin();
	std::map<Amplitude::Channel, std::string>::iterator it;
	std::advance(itS, 2);
	std::advance(itE, 29);

	for (it = itS ; it != itE; ++it)
		vName.push_back(it->first);

	return vName;
}

bool DecayRates::IsAllowed(Channel Name)
{
	if (Channel_prev != Name)
	{
		LoadMass(Name);
		Channel_prev = Name;
	}

	double Limit = 0.0;
	for (unsigned int i = 0; i < vMass.size(); ++i)
		Limit += vMass.at(i);

	return (MassN() >= Limit);
}

//return the decay width (Gamma)
//
double DecayRates::Gamma(Channel Name)
{
	double Result = 0.0;

	//if (IsChanged())
	//	Reset();

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
		case _nEM:
			Result = nEM();
			break;
		case _nME:
			Result = nME();
			break;
		case _nMM:
			Result = nMM();
			break;
		case _nET:
			Result = nET();
			break;
		case _nTE:
			Result = nTE();
			break;
		case _nMT:
			Result = nMT();
			break;
		case _nTM:
			Result = nTM();
			break;
		case _nPI0:
			Result = nPI0();
			break;
		case _EPI:
			Result = EPI();
			break;
		case _MPI:
			Result = MPI();
			break;
		case _TPI:
			Result = TPI();
			break;
		case _EKA:
			Result = EKA();
			break;
		case _MKA:
			Result = MKA();
			break;
		case _EKAx:
			Result = EKAx();
			break;
		case _MKAx:
			Result = MKAx();
			break;
		case _nRHO0:
			Result = nRHO0();
			break;
		case _ERHO:
			Result = ERHO();
			break;
		case _MRHO:
			Result = MRHO();
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
			std::cerr << ShowChannel(Name) << ": channel unknown" << std::endl;
			Result = 0.0;
			break;
	}

	return (1.0 + !GetFermion()) * Result;	//!fermion is majorana, Gamma is twice as much
}

//Return Gamma_tot - Gamma of interest
//
double DecayRates::Other(Channel Name)
{
	if (Gamma(Name) < 0.0)
		return 0.0;
	else return Gamma(_ALL) - Gamma(Name);
}

//Return the branching ration
//
double DecayRates::Branch(Channel Name)
{
	if (Gamma(Name) < 0.0)
		return 0.0;
	else return Gamma(Name)/Gamma(_ALL);
}

//total decay width
double DecayRates::Total()
{
	return (nnn() + nGAMMA() +
		nEE() + nEM() + nME() + nMM() + nET() + nTE() + nMT() + nTM() +
		nPI0() + EPI() + MPI() + TPI() +
		EKA() + MKA() + 
		nRHO0() + ERHO() + MRHO() +
		EKAx() + MKAx() + 
		nETA() + nETAi() + nOMEGA() + nPHI() +
		ECHARM() );
}

//special here
double DecayRates::ExpALL()
{
	return (nEE() + nME() + nMM() +
		EPI() + MPI() +
		EKA() + MKA() +
		ERHO() + MRHO() );
}

//individual decay channels
//all mixing factors are factorised out
double DecayRates::nGAMMA()
{
	if (fnGAMMA < 0 || IsChanged())
		fnGAMMA = IsAllowed(_nGAMMA) ?
			  Const::fGF2 * MassN(5) * (27.0 * Const::fAem) /
			  (192.0 * 64.0 * Const::fPi4) : 0.0;

	return fnGAMMA * (Ue(2) + Um(2) + Ut(2));
}

double DecayRates::nnn()
{
	if (fnnn < 0 || IsChanged())
		fnnn = IsAllowed(_nnn) ?
		       Const::fGF2 * MassN(5) / (192.0 * Const::fPi3) : 0.0;

	return fnnn * (Ue(2) + Um(2) + Ut(2));
}

//M_Sterile > 2 M_Electron (always)
double DecayRates::nEE()
{
	if (fnEE_e < 0 || fnEE_o < 0 || IsChanged())
	{
		if (IsAllowed(_nEE))
			NeutrinoLeptonAA(fnEE_e, fnEE_o, M_Neutrino, M_Electron);
		else
		{
			fnEE_e = 0.0;
			fnEE_o = 0.0;
		}
	}

	return fnEE_e * Ue(2) + fnEE_o * (Um(2) + Ut(2));
}

//M_Sterile > M_Muon + M_Electron
double DecayRates::nEM()	//Antiparticle is Elec
{
	if (fnEM < 0 || IsChanged())
		fnEM = IsAllowed(_nEM) ? 
			NeutrinoLeptonAB(M_Neutrino, M_Electron, M_Muon) : 0.0;

	return fnEM * Ue(2);
}

double DecayRates::nME()	//Anti is Muon
{
	if (fnME < 0 || IsChanged())
		fnME = IsAllowed(_nME) ? 
			NeutrinoLeptonAB(M_Neutrino, M_Muon, M_Electron) : 0.0;

	return fnME * Um(2);
}

//M_Sterile > 2 M_Muon
double DecayRates::nMM()
{
	if (fnMM_m < 0 || fnMM_o < 0 || IsChanged())
	{
		if (IsAllowed(_nMM)) 
			NeutrinoLeptonAA(fnMM_m, fnMM_o, M_Neutrino, M_Muon);
		else
		{
			fnMM_m  = 0.0;
			fnMM_o = 0.0;
		}
	}

	return fnMM_m * Um(2) + fnMM_o * (Ue(2) + Ut(2));
}

//M_Sterile > M_Tau + M_Electron
double DecayRates::nET()	//Antiparticle is Elec
{
	if (fnET < 0 || IsChanged())
		fnET = IsAllowed(_nET) ? 
			NeutrinoLeptonAB(M_Neutrino, M_Electron, M_Tau) : 0.0;

	return fnET * Ue(2);
}

double DecayRates::nTE()	//Anti is Tau
{
	if (fnTE < 0 || IsChanged())
		fnTE = IsAllowed(_nTE) ? 
			NeutrinoLeptonAB(M_Neutrino, M_Tau, M_Electron) : 0.0;

	return fnTE * Ut(2);
}

//M_Sterile > M_Tau + M_Muon
double DecayRates::nMT()	//Antiparticle is Muon
{
	if (fnMT < 0 || IsChanged())
		fnMT = IsAllowed(_nMT) ? 
			NeutrinoLeptonAB(M_Neutrino, M_Muon, M_Tau) : 0.0;

	return fnMT * Um(2);
}

double DecayRates::nTM()	//Anti is Tau
{
	if (fnTM < 0 || IsChanged())
		fnTM = IsAllowed(_nTM) ? 
			NeutrinoLeptonAB(M_Neutrino, M_Tau, M_Muon) : 0.0;

	return fnTM * Ut(2);
}

//M_Sterile > M_Pion0
double DecayRates::nPI0()
{
	if (fnPI0 < 0 || IsChanged())
		fnPI0 = IsAllowed(_nPI0) ?
			pow(Const::fDPion, 2) * NeutrinoPseudoMeson(M_Neutrino, M_Pion0) : 0.0;

	return fnPI0 * (Ue(2) + Um(2) + Ut(2));
}

//M_Sterile > M_Pion
double DecayRates::EPI()
{
	if (fEPI < 0 || IsChanged())
	{
		fEPI = IsAllowed(_EPI) ?
			pow(Const::fU_ud * Const::fDPion, 2) * LeptonPseudoMeson(M_Electron, M_Pion) : 0.0;
	}
	
	return fEPI * Ue(2);
}

//M_Sterile > M_Pion + M_Muon
double DecayRates::MPI()
{
	if (fMPI < 0 || IsChanged())
		fMPI = IsAllowed(_MPI) ?
			pow(Const::fU_ud * Const::fDPion, 2) * LeptonPseudoMeson(M_Muon, M_Pion) : 0.0;
	
	return fMPI * Um(2);
}

//M_Sterile > M_Tau + M_Pion
double DecayRates::TPI()
{
	if (fTPI < 0 || IsChanged())
		fTPI = IsAllowed(_TPI) ?
			pow(Const::fU_ud * Const::fDPion, 2) * LeptonPseudoMeson(M_Tau, M_Pion) : 0.0;
	
	return fTPI * Ut(2);
}

//M_Sterile > M_Kaon + M_Electron
double DecayRates::EKA()
{
	if (fEKA < 0 || IsChanged())
		fEKA = IsAllowed(_EKA) ?
			pow(Const::fU_us * Const::fDKaon, 2) * LeptonPseudoMeson(M_Electron, M_Kaon) : 0.0;

	return fEKA * Ue(2);
}

//M_Sterile > M_Kaon + M_Muon
double DecayRates::MKA()
{
	if (fMKA < 0 || IsChanged())
		fMKA = IsAllowed(_MKA) ?
			pow(Const::fU_us * Const::fDKaon, 2) * LeptonPseudoMeson(M_Muon, M_Kaon) : 0.0;

	return fMKA * Um(2);
}

//M_Sterile > M_Rho
double DecayRates::nRHO0()
{
	if (fnRHO0 < 0 || IsChanged())
		fnRHO0 = IsAllowed(_nRHO0) ? 
			pow(Const::fDRho * Const::fVRho, 2) * NeutrinoVectorMeson(M_Neutrino, M_Rho0) : 0.0;

	return fnRHO0 * (Ue(2) + Um(2) + Ut(2));
}

//M_Sterile > M_Rho + M_Electron 
double DecayRates::ERHO()
{
	if (fERHO < 0 || IsChanged())
		fERHO = IsAllowed(_ERHO) ? 
			pow(Const::fU_ud * Const::fDRho, 2) * LeptonVectorMeson(M_Electron, M_Rho) : 0.0;

	return fERHO * Ue(2);
}

//M_Sterile > M_Rho + M_Muon 
double DecayRates::MRHO()
{
	if (fMRHO < 0 || IsChanged())
		fMRHO = IsAllowed(_MRHO) ? 
			pow(Const::fU_ud * Const::fDRho, 2) * LeptonVectorMeson(M_Muon, M_Rho) : 0.0;

	return fMRHO * Um(2);
}

//M_Sterile > M_Kaon* + M_Electron 
double DecayRates::EKAx()
{
	if (fEKAx < 0 || IsChanged())
		fEKAx = IsAllowed(_EKAx) ? 
			pow(Const::fU_us * Const::fDKaonx, 2) * LeptonVectorMeson(M_Electron, M_Kaonx) : 0.0;

	return fEKAx * Ue(2);
}

//M_Sterile > M_Kaon* + M_Muon 
double DecayRates::MKAx()
{
	if (fMKAx < 0 || IsChanged())
		fMKAx = IsAllowed(_MKAx) ? 
			pow(Const::fU_us * Const::fDKaonx, 2) * LeptonVectorMeson(M_Muon, M_Kaonx) : 0.0;

	return fMKAx * Um(2);
}

//M_Sterile > M_Eta
double DecayRates::nETA()
{
	if (fnETA < 0 || IsChanged())
		fnETA = IsAllowed(_nETA) ? 
			pow(Const::fDEta, 2) * NeutrinoPseudoMeson(M_Neutrino, M_Eta) : 0.0;

	return fnETA * (Ue(2) + Um(2) + Ut(2));
}

//M_Sterile > M_Eta'
double DecayRates::nETAi()
{
	if (fnETAi < 0 || IsChanged())
		fnETAi = IsAllowed(_nETAi) ? 
			pow(Const::fDEtai, 2) * NeutrinoPseudoMeson(M_Neutrino, M_Etai) : 0.0;

	return fnETAi * (Ue(2) + Um(2) + Ut(2));
}

//M_Sterile > M_Omega 
double DecayRates::nOMEGA()
{
	if (fnOMEGA < 0 || IsChanged())
		fnOMEGA = IsAllowed(_nOMEGA) ? 
			pow(Const::fDOmega * Const::fVOmega, 2) * NeutrinoVectorMeson(M_Neutrino, M_Omega) : 0.0;

	return fnOMEGA * (Ue(2) + Um(2) + Ut(2));
}

//M_Sterile > M_Phi 
double DecayRates::nPHI()
{
	if (fnPHI < 0 || IsChanged())
		fnPHI = IsAllowed(_nPHI) ? 
			pow(Const::fDPhi * Const::fVPhi, 2) * NeutrinoVectorMeson(M_Neutrino, M_Phi) : 0.0;

	return fnPHI * (Ue(2) + Um(2) + Ut(2));
}

//M_Sterile > M_Charm + M_Electron
double DecayRates::ECHARM()
{
	if (fECHARM < 0 || IsChanged())
		fECHARM = IsAllowed(_ECHARM) ?
			pow(Const::fU_cd * Const::fDCharm, 2) * LeptonPseudoMeson(M_Electron, M_Charm) : 0.0;

	return fECHARM * Ue(2);
}

/////////////////
//Generic decay//
/////////////////
//
//CC version possible
double DecayRates::LeptonPseudoMeson(double M_Lepton, double M_Meson)
{
	SetMass(MassN());
	double dML2 = M_Lepton*M_Lepton/Mass(2);
	double dMM2 = M_Meson*M_Meson/Mass(2);

	return I_LeptonPseudoMeson(dML2, dMM2);
}

//integrated over angular dep.
double DecayRates::I_LeptonPseudoMeson(double x, double y)
{
	//double cos0 = theta < 0 ? 0.0 : cos(theta);
	//double fc = theta < 0.0 ? 2.0 : 1.0;

	double M2 = M2_LeptonPseudoMeson(0, x, y);
	return dGammad0_2B(M2, x, y);
}

double DecayRates::NeutrinoPseudoMeson(double M_Neut, double M_Meson)
{
	SetMass(MassN());
	double dMn2 = M_Neut*M_Neut/Mass(2);
	double dMM2 = M_Meson*M_Meson/Mass(2);

	return I_NeutrinoPseudoMeson(dMn2, dMM2);
}

//integrated over angular dep.
double DecayRates::I_NeutrinoPseudoMeson(double x, double y)
{
	//double cos0 = theta < 0 ? 0.0 : cos(theta);
	//double fc = theta < 0.0 ? 2.0 : 1.0;

	double M2 = M2_NeutrinoPseudoMeson(0, x, y);
	return dGammad0_2B(M2, x, y);
}

//CC version possible
double DecayRates::LeptonVectorMeson(double M_Lepton, double M_Meson)
{
	SetMass(MassN());
	double dML2 = M_Lepton*M_Lepton/Mass(2);
	double dMM2 = M_Meson*M_Meson/Mass(2);

	return I_LeptonVectorMeson(dML2, dMM2);
}

//integrated over angular dep.
double DecayRates::I_LeptonVectorMeson(double x, double y)
{
	//double cos0 = theta < 0 ? 0.0 : cos(theta);
	//double fc = theta < 0.0 ? 2.0 : 1.0;

	double M2 = M2_LeptonVectorMeson(0, x, y);
	return dGammad0_2B(M2, x, y);
}

double DecayRates::NeutrinoVectorMeson(double M_Neut, double M_Meson)
{
	SetMass(MassN());
	double dMn2 = M_Neut*M_Neut/Mass(2);
	double dMM2 = M_Meson*M_Meson/Mass(2);

	return  I_NeutrinoVectorMeson(dMn2, dMM2);
}

//integrated over angular dep.
double DecayRates::I_NeutrinoVectorMeson(double x, double y)
{
	//double cos0 = theta < 0 ? 0.0 : cos(theta);
	//double fc = theta < 0.0 ? 2.0 : 1.0;

	double M2 = M2_NeutrinoVectorMeson(0, x, y);
	return dGammad0_2B(M2, x, y);
}

double DecayRates::NeutrinoLeptonAA(double &Amp_Ul, double &Amp_Un, double M_Neut, double M_Lepton)
{
	SetMass(MassN());
	double dMn2 = M_Neut*M_Neut/Mass(2);
	double dML2 = M_Lepton*M_Lepton/Mass(2);

	//neutrino flavour is the same as leptons -> both Z and W interaction
	double gL_CC = -0.5 + Const::fSin2W + (2*GetParticle() - 1);	//times U(lepton flavour)
	double gR_CC = Const::fSin2W;

	//neutrino flavour is different from leptons -> Z interaction
	double gL_NC = -0.5 + Const::fSin2W;	//times U(neutrino flavour)
	double gR_NC = Const::fSin2W;

	Amp_Ul = NeutrinoLeptonLepton(dMn2, dML2, dML2, gL_CC, gR_CC);
	Amp_Un = NeutrinoLeptonLepton(dMn2, dML2, dML2, gL_NC, gR_NC);

	return 0.0;
}

//CC version also available
double DecayRates::NeutrinoLeptonAB(double M_Neut, double M_LeptonA, double M_LeptonB)
{
	SetMass(MassN());
	double dMn2 = M_Neut*M_Neut/Mass(2);
	double dMA2 = M_LeptonA*M_LeptonA/Mass(2);
	double dMB2 = M_LeptonB*M_LeptonB/Mass(2);

	//for CC
	double gL = 1.0;
	double gR = 0.0;

	return NeutrinoLeptonLepton(dMn2, dMA2, dMB2, gL, gR);
}

double DecayRates::NeutrinoLeptonLepton(double x, double y, double z, double gL, double gR)
{
	double M2 = I_NeutrinoLeptonLepton(x, y, z, gL, gR);
	return dGammad2_3B(M2);
}

double DecayRates::I_NeutrinoLeptonLepton(double x, double y, double z, double gL, double gR)//, double theta)
{
	F_var.clear();

	F_var.push_back(x);	//0
	F_var.push_back(y);	//1
	F_var.push_back(z);	//2
	F_var.push_back(gL);	//3
	F_var.push_back(gR);	//4

	//F_var.push_back(theta);	//3

	SetFunction(&DecayRates::I_NeutrinoLeptonLepton_s);
	return Inte::BooleIntegration(this); 
}

double DecayRates::I_NeutrinoLeptonLepton_s(double s)
{
	double x = F_var.at(0);
	double y = F_var.at(1);
	double z = F_var.at(2);
	double gL = F_var.at(3);
	double gR = F_var.at(4);

	double s_ = s, t_ = s, u_ = s;
	double fcu = Limit(u_, y, z, x);
	double fct = Limit(t_, z, x, y);
	double fcs = Limit(s_, x, y, z);

	return     gL * gL * fcs * M2_WW(s_, 0, x, y, z) +
	           gR * gR * fct * M2_WW(t_, 0, x, y, z) +
	       2 * gL * gR * fcu * M2_WZ(u_, 0, x, y, z);
}

/////////////////////////
//end of generic decays//
/////////////////////////

void DecayRates::Reset()
{
	fnnn    = -1.0;
	fnGAMMA = -1.0;
	fnEE_e  = -1.0;
	fnEE_o  = -1.0;
	fnEM    = -1.0;
	fnME    = -1.0;
	fnMM_m  = -1.0;
	fnMM_o  = -1.0;
	fnET    = -1.0;
	fnTE    = -1.0;
	fnMT    = -1.0;
	fnTM    = -1.0;
	fnPI0   = -1.0;
	fEPI    = -1.0;
	fMPI    = -1.0;
	fTPI    = -1.0;
	fEKA    = -1.0;
	fMKA    = -1.0;
	fnRHO0  = -1.0;
	fERHO   = -1.0;
	fMRHO   = -1.0;
	fEKAx   = -1.0;
	fMKAx   = -1.0;
	fnETA   = -1.0;
	fnETAi  = -1.0;
	fnOMEGA = -1.0;
	fnPHI   = -1.0;
	fECHARM = -1.0;
}

void DecayRates::SetFunction(double (DecayRates::*FF)(double))
{
	double (Amplitude::*Function)(double) = 
		static_cast<double (Amplitude::*)(double)>(FF); // ok!
	Amplitude::SetFunction(Function);
}
