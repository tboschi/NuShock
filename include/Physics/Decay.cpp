#include "Decay.h"

Decay::Decay()	:
{
}

bool Decay::IsAllowed(Channel Name)
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
double Decay::Gamma(Channel Name)
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
	else return Gamma(Name)/Gamma(_ALL);
}

//total decay width
double Decay::Total()
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
	if (IsAllowed() && (fnnn < 0 || IsChanged()))
	{
		if (IsAllowed(_nnn))
			fnnn = Const::fGF2 * pow(GetMass(), 5) /
			       (96.0 * Const::fPi3);
		else fnnn = 0.0;
	}

	return fnnn * (Ue*Ue + Um*Um + Ut*Ut);
}

double Decay::nGAMMA()
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
double Decay::nEE()
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
double Decay::nEMU()	//Antiparticle is Elec
{
	if (fnEMU < 0 || IsChanged())
	{
		if (IsAllowed(_nEMU))
			fnEMU = NeutrinoLeptonAB(M_Neutrino, M_Electron, M_Muon);
		else
			fnEMU = 0.0;
	}

	return fnEMU * Ue*Ue;
}

double Decay::nMUE()	//Anti is Muon
{
	if (fnMUE < 0 || IsChanged())
	{
		if (IsAllowed(_nMUE))
			fnMUE = NeutrinoLeptonAB(M_Neutrino, M_Muon, M_Electron);
		else
			fnMUE = 0.0;
	}

	return fnMUE * Um*Um;
}

//M_Sterile > 2 M_Muon
double Decay::nMUMU(double )
{
	if (fnMUMU_m < 0 || fnMUMU_a < 0 || IsChanged())
	{
		if (IsAllowed(_nMUMU)) 
			NeutrinoLeptonAA(fnMUMU_m, fnMUMU_a, M_Neutrino, M_Muon);
		else
		{
			fnMUMU_m  = 0.0;
			fnMUMU_a = 0.0;
		}
	}

	return (fnMUMU_m + fnMUMU_a) * Um*Um + fnMUMU_a * (Ue*Ue + Ut*Ut);
}

//M_Sterile > M_Tau + M_Electron
double Decay::nET()	//Antiparticle is Elec
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

double Decay::nTE()	//Anti is Tau
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
double Decay::nMUT()	//Antiparticle is Muon
{
	if (fnMUT < 0 || IsChanged())
	{
		if (IsAllowed(_nMUT))
			fnMUT = NeutrinoLeptonAB(M_Neutrino, M_Muon, M_Tau);
		else
			fnMUT = 0.0;
	}

	return fnMUT * Um*Um;
}

double Decay::nTMU()	//Anti is Tau
{
	if (fnTMU < 0 || IsChanged())
	{
		if (IsAllowed(_nTMU))
			fnTMU = NeutrinoLeptonAB(M_Neutrino, M_Tau, M_Muon);
		else
			fnTMU = 0.0;
	}

	return fnTMU * Ut*Ut;
}

//M_Sterile > M_Pion0
double Decay::nPI0()
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
double Decay::EPI()
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
double Decay::MUPI()
{
	if (fMUPI < 0 || IsChanged())
	{
		if (IsAllowed(_MUPI))
			fMUPI =  pow(Const::fU_ud * Const::fDPion, 2) * LeptonPseudoMeson(M_Muon, M_Pion);
		else
			fMUPI = 0.0;
	}
	
	return fMUPI * Um*Um;
}

//M_Sterile > M_Tau + M_Pion
double Decay::TPI()
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
double Decay::EKA()
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
double Decay::MUKA()
{
	if (fMUKA < 0 || IsChanged())
	{
		if (IsAllowed(_MUKA))
			fMUKA =  pow(Const::fU_us * Const::fDKaon, 2) * LeptonPseudoMeson(M_Muon, M_Kaon);
		else
			fMUKA = 0.0;
	}

	return fMUKA * Um*Um;
}

//M_Sterile > M_Rho
double Decay::nRHO0()
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
double Decay::ERHO()
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
double Decay::MURHO()
{
	if (fMURHO < 0 || IsChanged())
	{
		if (IsAllowed(_MURHO))
			fMURHO = pow(Const::fU_ud * Const::fDRho, 2) * LeptonVectorMeson(M_Muon, M_Rho);	//check
		else
			fMURHO = 0.0;
	}

	return fMURHO * Um*Um;
}

//M_Sterile > M_Kaon* + M_Electron 
double Decay::EKAx()
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
double Decay::MUKAx()
{
	if (fMUKAx < 0 || IsChanged())
	{
		if (IsAllowed(_MUKAx))
			fMUKAx = pow(Const::fU_us * Const::fDKaonx, 2) * LeptonVectorMeson(M_Muon, M_Kaonx);	//check
		else
			fMUKAx = 0.0;
	}

	return fMUKAx * Um*Um;
}

//M_Sterile > M_Eta
double Decay::nETA()
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
double Decay::nETAi()
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
double Decay::nOMEGA()
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
double Decay::nPHI()
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
double Decay::ECHARM()
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
double Decay::LeptonPseudoMeson(double M_Lepton, double M_Meson)
{
	double dML2 = M_Lepton*M_Lepton/GetMass()/GetMass();
	double dMM2 = M_Meson*M_Meson/GetMass()/GetMass();

	return I_LeptonPseudoMeson(dML2, dMM2);
}

//integrated over angular dep.
double Decay::I_LeptonPseudoMeson(double x, double y)
{
	//double cos0 = theta < 0 ? 0.0 : cos(theta);
	//double fc = theta < 0.0 ? 2.0 : 1.0;

	double M2 = M2_LeptonPseudoMeson(x, y, 0);
	return dGammad0_2B(M2, x, y);
}

double Decay::NeutrinoPseudoMeson(double M_Meson, double fDecay2)
{
	double dMN2 = M_Neutrino*M_Neutrino/GetMass()/GetMass();
	double dMM2 = M_Meson*M_Meson/GetMass()/GetMass();

	return fDecay2 * I_NeutrinoPseudoMeson(dMN2, dMM2);
}

//integrated over angular dep.
double Decay::I_NeutrinoPseudoMeson(double x, double y)
{
	//double cos0 = theta < 0 ? 0.0 : cos(theta);
	//double fc = theta < 0.0 ? 2.0 : 1.0;

	double M2 = M2_NeutrinoPseudoMeson(x, y, 0);
	return dGammad0_2B(M2, x, y);
}

//CC version possible
double Decay::LeptonVectorMeson(double M_Lepton, double M_Meson)
{
	double dML2 = M_Lepton*M_Lepton/GetMass()/GetMass();
	double dMM2 = M_Meson*M_Meson/GetMass()/GetMass();

	return I_LeptonVectorMeson(dMN2, dMM2);
}

//integrated over angular dep.
double Decay::I_LeptonVectorMeson(double x, double y)
{
	//double cos0 = theta < 0 ? 0.0 : cos(theta);
	//double fc = theta < 0.0 ? 2.0 : 1.0;

	double M2 = M2_LeptonVectorMeson(x, y, 0);
	return dGammad0_2B(M2, x, y);
}

double Decay::NeutrinoVectorMeson(double M_Neut, double M_Meson)
{
	double dMN2 = M_Neut*M_Neut/GetMass()/GetMass();
	double dMM2 = M_Meson*M_Meson/GetMass()/GetMass();

	return  I_NeutrinoVectorMeson(dMN2, dMM2);
}

//integrated over angular dep.
double Decay::I_LeptonVectorMeson(double x, double y)
{
	//double cos0 = theta < 0 ? 0.0 : cos(theta);
	//double fc = theta < 0.0 ? 2.0 : 1.0;

	double M2 = M2_NeutrinoVectorMeson(x, y, 0);
	return dGammad0_2B(M2, x, y);
}

double Decay::NeutrinoLeptonAA(double &Amp_Ul, double &Amp_Un, double M_Neut, double M_Lepton)
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
double Decay::NeutrinoLeptonAB(double M_Neut, double M_LeptonA, double M_LeptonB)
{
	double dMN2 = M_Neut*M_Neut/GetMass()/GetMass();
	double dMA2 = M_LeptonA*M_LeptonB/GetMass()/GetMass();
	double dMB2 = M_LeptonB*M_LeptonB/GetMass()/GetMass();

	//for CC
	double gL = 2;
	double gR = 0;

	return NeutrinoLeptonLepton(dMN2, dMA2, dMB2, gL, gR);
}

double Decay::NeutrinoLeptonLepton(double x, double y, double z, double gL, double gR)
{
	double M2 = (gL * gL + gR * gR) * I_WW(x, y, z) +
		    (2 * gL * gR) * I_WZ(x, y, z);

	return dGammad2_3B(M2, x, y);
}

double Decay::I_WW(double x, double y, double z)//, double theta)
{
	F_var.clear();

	F_var.push_back(x);	//0
	F_var.push_back(y);	//1
	F_var.push_back(z);	//2
	//F_var.push_back(theta);	//3

	SetFunction(&I_WW_s);				//Integrand will fix the integration volume
	return Inte::BooleIntegration(this); 
}

double Decay::I_WW_s(double s)
{
	const double &x = F_var.at(0);
	const double &y = F_var.at(1);
	const double &z = F_var.at(2);

	double t = 0;
	double fc = Limit(s, t, x, y, z);

	//double cos0 = theta < 0 ? 0.0 : cos(theta);		//means I am integrating on theta only
	//double fc = theta < 0.0 ? 2.0 : 1.0;			//so I remove cos0 terms and multiply by 2

	return fc * M2_WW(x, y, z, s, 0);
}

double Decay::I_WZ(double x, double y, double z)//, double theta)
{
	F_var.clear();

	F_var.push_back(x);	//0
	F_var.push_back(y);	//1
	F_var.push_back(z);	//2
	//F_var.push_back(theta);	//3

	SetFunction(&I_WZ_s);
	return Inte::BooleIntegration(this); 
}

double Decay::I_WZ_s(double s)
{
	F_var.push_back(s);

	SetFunction(&I_WZ_t);
	double Ret = Inte::BooleIntegration(this);

	F_var.pop_back();
	SetFunction(&I_WZ_s);
	return Ret;
}

double Decay::I_WZ_t(double t)
{
	const double &x = F_var.at(0);
	const double &y = F_var.at(1);
	const double &z = F_var.at(2);
	//const double &theta = F_var.at(3);
	double &s = F_var.at(3);

	double fc = Limit(s, t, x, y, z);
	//double cos0 = theta < 0 ? 0.0 : cos(theta);
	//double fc = theta < 0.0 ? 2.0 : 1.0;

	return fc * M2_WZ(x, y, z, s, t, 0, 0);
}

/////////////////////////
//end of generic decays//
/////////////////////////


/*
std::vector<std::string> Decay::ListChannels()
{
	std::map<std::string, ChannelName>::iterator it;
	std::vector<std::string> vList;
	for (it = mapChannel.begin(); it != mapChannel.end(); ++it)
		vList.push_back(it->first);
	return vList;
}
*/

vecbool Decay::IsChanged()
{
	bool Ret = (fabs(GetMass() - M_Sterile_prev) > 1e-9);

	M_Sterile_prev = GetMass();

	//Reset decay widths if changed
	if (Ret)
	{
		fnnn      = -1.0;
                fnGAMMA   = -1.0;
                fnEE_e    = -1.0;
                fnEE_a   = -1.0;
                fnEMU     = -1.0;
                fnMUE     = -1.0;
                fnMUMU_m  = -1.0;
                fnMUMU_a = -1.0;
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

