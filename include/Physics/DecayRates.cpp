#include "DecayRatesRates.h"

DecayRates::DecayRates()	:
{
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

	return (GetMass() >= Limit);
}

//return the decay width (Gamma)
//
double DecayRates::Gamma(Channel Name)
{
	double Result = 0.0;

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
			Result = 0.0;
			break;
	}

	return Result;
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

//total decay width
double DecayRates::Total()
{
	return (nnn() + nGAMMA() +
		nEE() + nEM() + nME() + nMM() + nET() + nTE() + nMT() + nTM() +
		nPI0() + EPI() + + MPI() + TPI() +
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
	if (fnEE_e < 0 || fnEE_a < 0 || IsChanged())
	{
		if (IsAllowed(_nEE))
			NeutrinoLeptonAA(fnEE_e, fnEE_a, M_Neutrino, M_Electron);
		else
		{
			fnEE_e  = 0.0;
			fnEE_a = 0.0;
		}
	}

	return (fnEE_e + fnEE_a) * Ue*Ue + fnEE_a * (Um*Um + Ut*Ut);
}

//M_Sterile > M_Muon + M_Electron
double DecayRates::nEM()	//Antiparticle is Elec
{
	if (fnEM < 0 || IsChanged())
	{
		if (IsAllowed(_nEM))
			fnEM = NeutrinoLeptonAB(M_Neutrino, M_Electron, M_Muon);
		else
			fnEM = 0.0;
	}

	return fnEM * Ue*Ue;
}

double DecayRates::nME()	//Anti is Muon
{
	if (fnME < 0 || IsChanged())
	{
		if (IsAllowed(_nME))
			fnME = NeutrinoLeptonAB(M_Neutrino, M_Muon, M_Electron);
		else
			fnME = 0.0;
	}

	return fnME * Um*Um;
}

//M_Sterile > 2 M_Muon
double DecayRates::nMM(double )
{
	if (fnMM_m < 0 || fnMM_a < 0 || IsChanged())
	{
		if (IsAllowed(_nMM)) 
			NeutrinoLeptonAA(fnMM_m, fnMM_a, M_Neutrino, M_Muon);
		else
		{
			fnMM_m  = 0.0;
			fnMM_a = 0.0;
		}
	}

	return (fnMM_m + fnMM_a) * Um*Um + fnMM_a * (Ue*Ue + Ut*Ut);
}

//M_Sterile > M_Tau + M_Electron
double DecayRates::nET()	//Antiparticle is Elec
{
	if (fnET < 0 || IsChanged())
	{
		if (IsAllowed(_nET))
			fnET = NeutrinoLeptonAB(M_Neutrino, M_Electron, M_Tau);
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
			fnTE = NeutrinoLeptonAB(M_Neutrino, M_Tau, M_Electron);
		else
			fnTE = 0.0;
	}

	return fnTE * Ut*Ut;
}

//M_Sterile > M_Tau + M_Muon
double DecayRates::nMT()	//Antiparticle is Muon
{
	if (fnMT < 0 || IsChanged())
	{
		if (IsAllowed(_nMT))
			fnMT = NeutrinoLeptonAB(M_Neutrino, M_Muon, M_Tau);
		else
			fnMT = 0.0;
	}

	return fnMT * Um*Um;
}

double DecayRates::nTM()	//Anti is Tau
{
	if (fnTM < 0 || IsChanged())
	{
		if (IsAllowed(_nTM))
			fnTM = NeutrinoLeptonAB(M_Neutrino, M_Tau, M_Muon);
		else
			fnTM = 0.0;
	}

	return fnTM * Ut*Ut;
}

//M_Sterile > M_Pion0
double DecayRates::nPI0()
{
	if (fnPI0 < 0 || IsChanged())
	{
		if (IsAllowed(_nPI0))
			fnPI0 = pow(Const::fDPion, 2) * NeutrinoPseudoMeson(M_Neutrino, M_Pion);	//check
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
			fEPI =  pow(Const::fU_ud * Const::fDPion, 2) * LeptonPseudoMeson(M_Electron, M_Pion);
		else
			fEPI = 0.0;
	}
	
	return fEPI * Ue*Ue;
}

//M_Sterile > M_Pion + M_Muon
double DecayRates::MPI()
{
	if (fMPI < 0 || IsChanged())
	{
		if (IsAllowed(_MPI))
			fMPI =  pow(Const::fU_ud * Const::fDPion, 2) * LeptonPseudoMeson(M_Muon, M_Pion);
		else
			fMPI = 0.0;
	}
	
	return fMPI * Um*Um;
}

//M_Sterile > M_Tau + M_Pion
double DecayRates::TPI()
{
	if (fTPI < 0 || IsChanged())
	{
		if (IsAllowed(_TPI))
			fTPI =  pow(Const::fU_ud * Const::fDPion, 2) * LeptonPseudoMeson(M_Tau, M_Pion);
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
			fEKA =  pow(Const::fU_us * Const::fDKaon, 2) * LeptonPseudoMeson(M_Electron, M_Kaon);
		else
			fEKA = 0.0;
	}

	return fEKA * Ue*Ue;
}

//M_Sterile > M_Kaon + M_Muon
double DecayRates::MKA()
{
	if (fMKA < 0 || IsChanged())
	{
		if (IsAllowed(_MKA))
			fMKA =  pow(Const::fU_us * Const::fDKaon, 2) * LeptonPseudoMeson(M_Muon, M_Kaon);
		else
			fMKA = 0.0;
	}

	return fMKA * Um*Um;
}

//M_Sterile > M_Rho
double DecayRates::nRHO0()
{
	if (fnRHO0 < 0 || IsChanged())
	{
		if (IsAllowed(_nRHO0))
			fnRHO0 = pow(Const::fDRho * fVRho, 2) * NeutrinoVectorMeson(M_Neutrino, M_Rho0);	//check
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
			fERHO = pow(Const::fU_ud * Const::fDRho, 2) * LeptonVectorMeson(M_Electron, M_Rho);	//check
		else
			fERHO = 0.0;
	}

	return fERHO * Ue*Ue;
}

//M_Sterile > M_Rho + M_Muon 
double DecayRates::MRHO()
{
	if (fMRHO < 0 || IsChanged())
	{
		if (IsAllowed(_MRHO))
			fMRHO = pow(Const::fU_ud * Const::fDRho, 2) * LeptonVectorMeson(M_Muon, M_Rho);	//check
		else
			fMRHO = 0.0;
	}

	return fMRHO * Um*Um;
}

//M_Sterile > M_Kaon* + M_Electron 
double DecayRates::EKAx()
{
	if (fEKAx < 0 || IsChanged())
	{
		if (IsAllowed(_EKAx))
			fEKAx = pow(Const::fU_us * Const::fDKaonx, 2) * LeptonVectorMeson(M_Electron, M_Kaonx);	//check
		else
			fEKAx = 0.0;
	}

	return fEKAx * Ue*Ue;
}

//M_Sterile > M_Kaon* + M_Muon 
double DecayRates::MKAx()
{
	if (fMKAx < 0 || IsChanged())
	{
		if (IsAllowed(_MKAx))
			fMKAx = pow(Const::fU_us * Const::fDKaonx, 2) * LeptonVectorMeson(M_Muon, M_Kaonx);	//check
		else
			fMKAx = 0.0;
	}

	return fMKAx * Um*Um;
}

//M_Sterile > M_Eta
double DecayRates::nETA()
{
	if (fnETA < 0 || IsChanged())
	{
		if (IsAllowed(_nETA))
			fnETA = pow(Const::fDEta, 2) * NeutrinoPseudoMeson(M_Neutrino, M_Eta);	//check
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
			fnETAi = pow(Const::fDEtai, 2) * NeutrinoPseudoMeson(M_Neutrino, M_Etai);	//check
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
			fnOMEGA = pow(Const::fDOmega * Const::fVOmega, 2) * NeutrinoVectorMeson(M_Neutirno, M_Omega);	//check
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
			fnPHI = pow(Const::fDPhi * Const::fVPhi, 2) * NeutrinoVectorMeson(M_Neutrino, M_Phi);	//check
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
			fECHARM = pow(Const::fU_cd * Const::fDCharm, 2) * LeptonPseudoMeson(M_Electron, M_Charm);
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
double DecayRates::LeptonPseudoMeson(double M_Lepton, double M_Meson)
{
	double dML2 = M_Lepton*M_Lepton/GetMass()/GetMass();
	double dMM2 = M_Meson*M_Meson/GetMass()/GetMass();

	return I_LeptonPseudoMeson(dML2, dMM2);
}

//integrated over angular dep.
double DecayRates::I_LeptonPseudoMeson(double x, double y)
{
	//double cos0 = theta < 0 ? 0.0 : cos(theta);
	//double fc = theta < 0.0 ? 2.0 : 1.0;

	double M2 = M2_LeptonPseudoMeson(x, y, 0);
	return dGammad0_2B(M2, x, y);
}

double DecayRates::NeutrinoPseudoMeson(double M_Neut, double M_Meson)
{
	double dMN2 = M_Neut*M_Neut/GetMass()/GetMass();
	double dMM2 = M_Meson*M_Meson/GetMass()/GetMass();

	return I_NeutrinoPseudoMeson(dMN2, dMM2);
}

//integrated over angular dep.
double DecayRates::I_NeutrinoPseudoMeson(double x, double y)
{
	//double cos0 = theta < 0 ? 0.0 : cos(theta);
	//double fc = theta < 0.0 ? 2.0 : 1.0;

	double M2 = M2_NeutrinoPseudoMeson(x, y, 0);
	return dGammad0_2B(M2, x, y);
}

//CC version possible
double DecayRates::LeptonVectorMeson(double M_Lepton, double M_Meson)
{
	double dML2 = M_Lepton*M_Lepton/GetMass()/GetMass();
	double dMM2 = M_Meson*M_Meson/GetMass()/GetMass();

	return I_LeptonVectorMeson(dMN2, dMM2);
}

//integrated over angular dep.
double DecayRates::I_LeptonVectorMeson(double x, double y)
{
	//double cos0 = theta < 0 ? 0.0 : cos(theta);
	//double fc = theta < 0.0 ? 2.0 : 1.0;

	double M2 = M2_LeptonVectorMeson(x, y, 0);
	return dGammad0_2B(M2, x, y);
}

double DecayRates::NeutrinoVectorMeson(double M_Neut, double M_Meson)
{
	double dMN2 = M_Neut*M_Neut/GetMass()/GetMass();
	double dMM2 = M_Meson*M_Meson/GetMass()/GetMass();

	return  I_NeutrinoVectorMeson(dMN2, dMM2);
}

//integrated over angular dep.
double DecayRates::I_LeptonVectorMeson(double x, double y)
{
	//double cos0 = theta < 0 ? 0.0 : cos(theta);
	//double fc = theta < 0.0 ? 2.0 : 1.0;

	double M2 = M2_NeutrinoVectorMeson(x, y, 0);
	return dGammad0_2B(M2, x, y);
}

double DecayRates::NeutrinoLeptonAA(double &Amp_Ul, double &Amp_Un, double M_Neut, double M_Lepton)
{
	double dMN2 = M_Neut*M_Neut/GetMass()/GetMass();
	double dML2 = M_Lepton*M_Lepton/GetMass()/GetMass();

	//for CCNC
	double gL_CC = 2 - 4*GetParticle();	//times U(lepton flavour)
	double gR_CC = 0;
	double gL_NC = -1 + 2*Const::fSin2W;	//times U(neutrino flavour)
	double gR_NC = 2*Const::fSin2W;

	double IntWW = I_WW(dMN2, dML2, dML2);	//W or Z mediation
	double IntWZ = I_WZ(dMN2, dML2, dML2);	//cross term for W + Z

	Amp_Ul = NeutrinoLeptonLepton(dMN2, dML2, dML2, gL_CC, gL_CC);
	Amp_Un = NeutrinoLeptonLepton(dMN2, dML2, dML2, gL_NC, gL_NC);

	return 0.0;
}

//CC version also available
double DecayRates::NeutrinoLeptonAB(double M_Neut, double M_LeptonA, double M_LeptonB)
{
	double dMN2 = M_Neut*M_Neut/GetMass()/GetMass();
	double dMA2 = M_LeptonA*M_LeptonB/GetMass()/GetMass();
	double dMB2 = M_LeptonB*M_LeptonB/GetMass()/GetMass();

	//for CC
	double gL = 2;
	double gR = 0;

	return NeutrinoLeptonLepton(dMN2, dMA2, dMB2, gL, gR);
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

	SetFunction(&I_NeutrinoLeptonLepton_s);				//Integrand will fix the integration volume
	return Inte::BooleIntegration(this); 
}

double DecayRates::I_NeutrinoLeptonLepton_s(double s)
{
	const double &x = F_var.at(0);
	const double &y = F_var.at(1);
	const double &z = F_var.at(2);
	const double &gL = F_var.at(3);
	const double &gR = F_var.at(4);

	double t = s, u = s;
	double fcu = Limit(u, y, z, x);
	double fct = Limit(t, z, x, y);
	double fcs = Limit(s, x, y, z);

	//double cos0 = theta < 0 ? 0.0 : cos(theta);		//means I am integrating on theta only
	//double fc = theta < 0.0 ? 2.0 : 1.0;			//so I remove cos0 terms and multiply by 2

	return gL * gL * fcs * M2_WW(x, y, z, s) +
	       gR * gR * fct * M2_WW(x, y, z, t) + 
	       2 * gL * gR * fcu * M2_WZ(x, y, z, u);
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

void DecayRates::Reset()
{
	fnnn    = -1.0;
	fnGAMMA = -1.0;
	fnEE_e  = -1.0;
	fnEE_a  = -1.0;
	fnEM    = -1.0;
	fnME    = -1.0;
	fnMM_m  = -1.0;
	fnMM_a  = -1.0;
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
