#include "Production.h"

Production::Production()
{
	Reset();
}

Amplitude::Channel Production::FindChannel(std::string Name)
{
	if (chMap.size() == 0)
		LoadMap();

	std::map<Amplitude::Channel, std::string>::iterator it = chMap.begin();
	for (std::advance(it, 29); it != chMap.end(); ++it)
	{
		if (it->second == Name)
			return it->first;
	}

	return _undefined;
}

std::vector<std::string> Production::ListChannel()
{
	if (chMap.size() == 0)
		LoadMap();

	std::vector<std::string> vName;

	std::map<Amplitude::Channel, std::string>::iterator it = chMap.begin();
	for (std::advance(it, 29); it != chMap.end(); ++it)
		vName.push_back(it->second);

	return vName;
}

bool Production::IsAllowed(Channel Name)
{
	if (Channel_prev != Name)
	{
		LoadMass(Name);
		Channel_prev = Name;
	}

	double Limit = vMass.at(0);
	for (unsigned int i = 1; i < vMass.size(); ++i)
		Limit -= vMass.at(i);

	return (Limit >= MassN());
}

double Production::Gamma(Channel Name)
{
	double Result = 0.0;

	switch(Name)
	{
		case _ALL:
			Result = Total();
			break;
		case _MuonE:
			Result = MuonE();
			break;
		case _MuonM:
			Result = MuonM();
			break;
		case _TauEE:
			Result = TauEE();
			break;
		case _TauET:
			Result = TauET();
			break;
		case _TauMM:
			Result = TauMM();
			break;
		case _TauMT:
			Result = TauMT();
			break;
		case _TauPion:
			Result = TauPion();
			break;
		case _PionE:
			Result = PionE();
			break;
		case _PionM:
			Result = PionM();
			break;
		case _KaonE:
			Result = KaonE();
			break;
		case _KaonM:
			Result = KaonM();
			break;
		case _CharmE:
			Result = CharmE();
			break;
		case _CharmM:
			Result = CharmM();
			break;
		case _CharmT:
			Result = CharmT();
			break;
		case _Kaon0E:
			Result = Kaon0E();
			break;
		case _Kaon0M:
			Result = Kaon0M();
			break;
		case _KaonCE:
			Result = KaonCE();
			break;
		case _KaonCM:
			Result = KaonCM();
			break;
		default:
			std::cerr << ShowChannel(Name) << ": channel unknown" << std::endl;
			Result = 0.0;
			break;
	}

	return Result;
}

double Production::Scale(Channel Name)
{
	double mass = MassN();
	int hel = Helicity();
	double GammaN = Gamma(Name);

	if (GammaN > 0)
	{
		SetMassN(0);
		SetHelicity(0);
		double Gamma0 = 2*Gamma(Name);

		SetMassN(mass);
		SetHelicity(hel);
		return GammaN/Gamma0;
	}
	else
		return 0.0;
}

double Production::Total()
{
	return (MuonE() + MuonM() + TauEE() + TauET() + TauMM() + TauMT() +
		PionE() + PionM() + KaonE() + KaonM() + CharmE() + CharmM() + CharmT() +
		Kaon0E() + Kaon0M() + KaonCE() + KaonCM());
}

double Production::MuonE()
{
	double *fAmp = NULL;
	bool Cond = false;

	if (MassN() > 0)
	{
		fAmp = &fMuonE;
		if (*fAmp < 0 || IsChanged())
			Cond = true;
	}
	else
	{
		fAmp = &lMuonE;
		if (*fAmp < 0)
			Cond = true;
	}

	if (Cond)
	{
		if (IsAllowed(_MuonE))
			*fAmp = AntiLeptonNeutrinoDecay(M_Muon, M_Electron, M_Neutrino);
		else
			*fAmp  = 0.0;
	}

	return *fAmp * Ue(2);
}

double Production::MuonM()
{
	double *fAmp = NULL;

	bool Cond = false;
	if (MassN() > 0)
	{
		fAmp = &fMuonM;
		if (*fAmp < 0 || IsChanged())
			Cond = true;
	}
	else
	{
		fAmp = &lMuonM;
		if (*fAmp < 0)
			Cond = true;
	}

	if (Cond)
	{
		if (IsAllowed(_MuonM))
			*fAmp = LeptonNeutrinoDecay(M_Muon, M_Electron, M_Neutrino);
		else
			*fAmp  = 0.0;
	}

	return *fAmp * Um(2);
}

double Production::TauEE()
{
	double *fAmp = NULL;
	bool Cond = false;

	if (MassN() > 0)
	{
		fAmp = &fTauEE;
		if (*fAmp < 0 || IsChanged())
			Cond = true;
	}
	else
	{
		fAmp = &lTauEE;
		if (*fAmp < 0)
			Cond = true;
	}

	if (Cond)
	{
		if (IsAllowed(_TauEE))
			*fAmp = AntiLeptonNeutrinoDecay(M_Tau, M_Electron, M_Neutrino);
		else
			*fAmp  = 0.0;
	}

	return *fAmp * Ue(2);
}

double Production::TauET()
{
	double *fAmp = NULL;
	bool Cond = false;

	if (MassN() > 0)
	{
		fAmp = &fTauET;
		if (*fAmp < 0 || IsChanged())
			Cond = true;
	}
	else
	{
		fAmp = &lTauET;
		if (*fAmp < 0)
			Cond = true;
	}

	if (Cond)
	{
		if (IsAllowed(_TauET))
			*fAmp = LeptonNeutrinoDecay(M_Tau, M_Electron, M_Neutrino);
		else
			*fAmp  = 0.0;
	}

	return *fAmp * Ut(2);
}

double Production::TauMM()
{
	double *fAmp = NULL;
	bool Cond = false;

	if (MassN() > 0)
	{
		fAmp = &fTauMM;
		if (*fAmp < 0 || IsChanged())
			Cond = true;
	}
	else
	{
		fAmp = &lTauMM;
		if (*fAmp < 0)
			Cond = true;
	}

	if (Cond)
	{
		if (IsAllowed(_TauMM))
			*fAmp = AntiLeptonNeutrinoDecay(M_Tau, M_Muon, M_Neutrino);
		else
			*fAmp  = 0.0;
	}

	return *fAmp * Um(2);
}

double Production::TauMT()
{
	double *fAmp = NULL;
	bool Cond = false;

	if (MassN() > 0)
	{
		fAmp = &fTauMT;
		if (*fAmp < 0 || IsChanged())
			Cond = true;
	}
	else
	{
		fAmp = &lTauMT;
		if (*fAmp < 0)
			Cond = true;
	}

	if (Cond)
	{
		if (IsAllowed(_TauMT))
			*fAmp = LeptonNeutrinoDecay(M_Tau, M_Muon, M_Neutrino);
		else
			*fAmp  = 0.0;
	}

	return *fAmp * Ut(2);
}

double Production::TauPion()
{
	double *fAmp = NULL;
	bool Cond = false;

	if (MassN() > 0)
	{
		fAmp = &fTauPion;
		if (*fAmp < 0 || IsChanged())
			Cond = true;
	}
	else
	{
		fAmp = &lTauPion;
		if (*fAmp < 0)
			Cond = true;
	}

	if (Cond)
	{
		if (IsAllowed(_TauPion))
			*fAmp = pow(Const::fU_ud * Const::fDPion, 2) * LeptonMesonDecay(M_Tau, M_Pion);
		else
			*fAmp  = 0.0;
	}

	return *fAmp * Ut(2);
}

double Production::PionE()
{
	double *fAmp = NULL;
	bool Cond = false;

	if (MassN() > 0)
	{
		fAmp = &fPionE;
		if (*fAmp < 0 || IsChanged())
			Cond = true;
	}
	else
	{
		fAmp = &lPionE;
		if (*fAmp < 0)
			Cond = true;
	}

	if (Cond)
	{
		if (IsAllowed(_PionE))
			*fAmp = pow(Const::fU_ud * Const::fDPion, 2) * MesonTwoDecay(M_Pion, M_Electron);
		else
			*fAmp  = 0.0;
	}

	return *fAmp * Ue(2);
}

double Production::PionM()
{
	double *fAmp = NULL;
	bool Cond = false;

	if (MassN() > 0)
	{
		fAmp = &fPionM;
		if (*fAmp < 0 || IsChanged())
			Cond = true;
	}
	else
	{
		fAmp = &lPionM;
		if (*fAmp < 0)
			Cond = true;
	}

	if (Cond)
	{
		if (IsAllowed(_PionM))
			*fAmp = pow(Const::fU_ud * Const::fDPion, 2) * MesonTwoDecay(M_Pion, M_Muon);
		else
			*fAmp  = 0.0;
	}

	return *fAmp * Um(2);
}

double Production::KaonE()
{
	double *fAmp = NULL;
	bool Cond = false;

	if (MassN() > 0)
	{
		fAmp = &fKaonE;
		if (*fAmp < 0 || IsChanged())
			Cond = true;
	}
	else
	{
		fAmp = &lKaonE;
		if (*fAmp < 0)
			Cond = true;
	}

	if (Cond)
	{
		if (IsAllowed(_KaonE))
			*fAmp = pow(Const::fU_us * Const::fDKaon, 2) * MesonTwoDecay(M_Kaon, M_Electron);
		else
			*fAmp  = 0.0;
	}

	else
	return *fAmp * Ue(2);
}

double Production::KaonM()
{
	double *fAmp = NULL;
	bool Cond = false;

	if (MassN() > 0)
	{
		fAmp = &fKaonM;
		if (*fAmp < 0 || IsChanged())
			Cond = true;
	}
	else
	{
		fAmp = &lKaonM;
		if (*fAmp < 0)
			Cond = true;
	}

	if (Cond)
	{
		if (IsAllowed(_KaonM))
			*fAmp = pow(Const::fU_us * Const::fDKaon, 2) * MesonTwoDecay(M_Kaon, M_Muon);
		else
			*fAmp  = 0.0;
	}

	return *fAmp * Um(2);
}

double Production::CharmE()
{
	double *fAmp = NULL;
	bool Cond = false;

	if (MassN() > 0)
	{
		fAmp = &fCharmE;
		if (*fAmp < 0 || IsChanged())
			Cond = true;
	}
	else
	{
		fAmp = &lCharmE;
		if (*fAmp < 0)
			Cond = true;
	}

	if (Cond)
	{
		if (IsAllowed(_CharmE))
			*fAmp = pow(Const::fU_cs * Const::fDCharm, 2) * MesonTwoDecay(M_CharmS, M_Electron);
		else
			*fAmp  = 0.0;
	}

	return *fAmp * Ue(2);
}

double Production::CharmM()
{
	double *fAmp = NULL;
	bool Cond = false;

	if (MassN() > 0)
	{
		fAmp = &fCharmM;
		if (*fAmp < 0 || IsChanged())
			Cond = true;
	}
	else
	{
		fAmp = &lCharmM;
		if (*fAmp < 0)
			Cond = true;
	}

	if (Cond)
	{
		if (IsAllowed(_CharmM))
			*fAmp = pow(Const::fU_cs * Const::fDCharm, 2) * MesonTwoDecay(M_CharmS, M_Muon);
		else
			*fAmp  = 0.0;
	}

	return *fAmp * Um(2);
}

double Production::CharmT()
{
	double *fAmp = NULL;
	bool Cond = false;

	if (MassN() > 0)
	{
		fAmp = &fCharmT;
		if (*fAmp < 0 || IsChanged())
			Cond = true;
	}
	else
	{
		fAmp = &lCharmT;
		if (*fAmp < 0)
			Cond = true;
	}

	if (Cond)
	{
		if (IsAllowed(_CharmT))
			*fAmp = pow(Const::fU_cs * Const::fDCharm, 2) * MesonTwoDecay(M_CharmS, M_Tau);
		else
			*fAmp  = 0.0;
	}

	return *fAmp * Ut(2);
}

double Production::Kaon0E()
{
	double *fAmp = NULL;
	bool Cond = false;

	if (MassN() > 0)
	{
		fAmp = &fKaon0E;
		if (*fAmp < 0 || IsChanged())
			Cond = true;
	}
	else
	{
		fAmp = &lKaon0E;
		if (*fAmp < 0)
			Cond = true;
	}

	if (Cond)
	{
		if (IsAllowed(_Kaon0E))
			*fAmp = pow(Const::fU_us * Const::fKaPi, 2) * 
				MesonThreeDecay(M_Kaon0, M_Pion, M_Electron, Const::fK0L_, Const::fK0L0);
		else
			*fAmp  = 0.0;
	}

	return *fAmp * Ue(2);
}

double Production::Kaon0M()
{
	double *fAmp = NULL;
	bool Cond = false;

	if (MassN() > 0)
	{
		fAmp = &fKaon0M;
		if (*fAmp < 0 || IsChanged())
			Cond = true;
	}
	else
	{
		fAmp = &lKaon0M;
		if (*fAmp < 0)
			Cond = true;
	}

	if (Cond)
	{
		if (IsAllowed(_Kaon0M))
			*fAmp = pow(Const::fU_us * Const::fKaPi, 2) * 
				MesonThreeDecay(M_Kaon0, M_Pion, M_Muon, Const::fK0L_, Const::fK0L0);
		else
			*fAmp  = 0.0;
	}

	return *fAmp * Um(2);
}

double Production::KaonCE()
{
	double *fAmp = NULL;
	bool Cond = false;

	if (MassN() > 0)
	{
		fAmp = &fKaonCE;
		if (*fAmp < 0 || IsChanged())
			Cond = true;
	}
	else
	{
		fAmp = &lKaonCE;
		if (*fAmp < 0)
			Cond = true;
	}

	if (Cond)
	{
		if (IsAllowed(_KaonCE))
			*fAmp = pow(Const::fU_us * Const::fKaPi, 2) / 2.0 * 
				MesonThreeDecay(M_Kaon, M_Pion0, M_Electron, Const::fKCL_, Const::fKCL0);
		else
			*fAmp  = 0.0;
	}

	return *fAmp * Ue(2);
}

double Production::KaonCM()
{
	double *fAmp = NULL;
	bool Cond = false;

	if (MassN() > 0)
	{
		fAmp = &fKaonCM;
		if (*fAmp < 0 || IsChanged())
			Cond = true;
	}
	else
	{
		fAmp = &lKaonCM;
		if (*fAmp < 0)
			Cond = true;
	}

	if (Cond)
	{
		if (IsAllowed(_KaonCM))
			*fAmp = pow(Const::fU_us * Const::fKaPi, 2) / 2.0 * 
				MesonThreeDecay(M_Kaon, M_Pion0, M_Muon, Const::fKCL_, Const::fKCL0);
		else
			*fAmp  = 0.0;
	}

	return *fAmp * Um(2);
}

/////////////////
//Generic decay//
/////////////////
//
//production from lepton -> neutrino is from same leptonic line
double Production::LeptonNeutrinoDecay(double M_Lepton0, double M_Lepton, double M_Neut)
{
	SetMass(M_Lepton0);
	double dML2 = M_Lepton*M_Lepton/Mass(2);
	double dMn2 = M_Neut*M_Neut/Mass(2);
	double dMN2 = MassN(2)/Mass(2);

	double M2 = I_LeptonNeutrino(dMn2, dML2, dMN2);
	return dGammad2_3B(M2);
}

//						c	  b	    a
double Production::I_LeptonNeutrino(double x, double y, double z)//, double theta)	//integrate first in y and then in x
{
	F_var.clear();

	F_var.push_back(x);	//0	//x
	F_var.push_back(y);	//1	//y
	F_var.push_back(z);	//2	//z
	//F_var.push_back(cos0);	//3	//theta

	SetFunction(&Production::I_LeptonNeutrino_u);
	return Inte::BooleIntegration(this); 
}

double Production::I_LeptonNeutrino_u(const double u)	//fixing one variable
{
	//aliases for clarity
	double &x = F_var.at(0);	//light neutrino	u
	double &y = F_var.at(1);	//lepton
	double &z = F_var.at(2);	//heavy neutrino
	//const double &cos0 = F_var.at(3);

	double u_ = u;
	double fc = Limit(u_, y, z, x);

	//double cos0 = theta < 0 ? 0.0 : cos(theta);	//theta < 0 means integration over theta
	//double fc = theta < 0 ? 2.0 : 1.0;		//so an overall factor of 2

	return fc * M2_LeptonNeutrino(x, y, z, u_);
}


//production from antilepton -> neutrino is from opposite leptonic line
double Production::AntiLeptonNeutrinoDecay(double M_Lepton0, double M_Lepton, double M_Neut)
{
	SetMass(M_Lepton0);
	double dML2 = M_Lepton*M_Lepton/Mass(2);
	double dMn2 = M_Neut*M_Neut/Mass(2);
	double dMN2 = MassN(2)/Mass(2);

	double M2 = I_AntiLeptonNeutrino(dMn2, dML2, dMN2);
	return dGammad2_3B(M2);
}
						//  c	      b		a
double Production::I_AntiLeptonNeutrino(double x, double y, double z)//, double theta)	//integrate first in y and then in x
{
	F_var.clear();

	F_var.push_back(x);	//0	//x
	F_var.push_back(y);	//1	//y
	F_var.push_back(z);	//2	//z
	//F_var.push_back(cos0);	//3	//theta

	SetFunction(&Production::I_AntiLeptonNeutrino_s);
	return Inte::BooleIntegration(this);
}

double Production::I_AntiLeptonNeutrino_s(double s)	//the term is written for a neutrino production
{								//therefore with heliciies inverted
	//aliases for clarity
	double &x = F_var.at(0);
	double &y = F_var.at(1);
	double &z = F_var.at(2);
	//const double &cos0 = F_var.at(5);

	//create S var
	double s_ = s;
	double fc = Limit(s_, x, y, z);

	//double cos0 = theta < 0 ? 0.0 : cos(theta);	//theta < 0 means integration over theta
	//double fc = theta < 0 ? 2.0 : 1.0;		//so an overall factor of 2

	return fc * M2_AntiLeptonNeutrino(x, y, z, s_);
}

double Production::LeptonMesonDecay(double M_Lepton, double M_Meson)
{
	SetMass(M_Lepton);
	double dMN2 = MassN(2)/Mass(2);
	double dMM2 = M_Meson*M_Meson/Mass(2);

	return I_LeptonMeson(dMN2, dMM2);
}

double Production::I_LeptonMeson(double x, double y)
{
	double M2 = M2_LeptonMeson(x, y);
	return dGammad0_2B(M2, x, y);
}

double Production::MesonTwoDecay(double M_Meson, double M_Lepton)
{
	SetMass(M_Meson);
	double dMN2 = MassN(2)/Mass(2);
	double dML2 = M_Lepton*M_Lepton/Mass(2);

	return I_MesonTwo(dMN2, dML2);
}

double Production::I_MesonTwo(double x, double y)	//symetric in x and y
{
	double M2 = M2_MesonTwo(x, y);
	return dGammad0_2B(M2, x, y);
}

double Production::MesonThreeDecay(double M_Meson0, double M_Meson, double M_Lepton, double L_, double L0)	//decay constant not important
{
	SetMass(M_Meson0);
	double dMM2 = M_Meson*M_Meson/Mass(2);
	double dML2 = M_Lepton*M_Lepton/Mass(2);
	double dMN2 = MassN(2)/Mass(2);	

	return I_MesonThree(dMM2, dML2, dMN2, L_, L0);
}

double Production::I_MesonThree(double x, double y, double z, double L_, double L0)	//no angle dependence
{
	F_var.clear();

	F_var.push_back(x);	//0	//c2
	F_var.push_back(y);	//1	//b2
	F_var.push_back(z);	//2	//a2
	//F_var.push_back(cos0);	//3	//theta
	F_var.push_back(L_);	//4	//linear dep for decay constant
	F_var.push_back(L0);	//5	//linear dep for decay constant


	int Trial, Fail;
	double Error, Chi2Prob;

	SetFunction(&Production::I_MesonThree_s);
	double Boole = Inte::BooleIntegration(this); 		//switch to Vega

	//SetFunction_D(&Production::I_MesonThree_D);
	//double Vegas = Inte::VegasIntegration(this, 2, Trial, Fail, Error, Chi2Prob);
	//std::cout << MassN() << "\t" << Boole << "\t" << Vegas << "\t" << Boole-Vegas << std::endl;

	return Boole;
}

double Production::I_MesonThree_s(const double s)	//fixing one variable
{
	F_var.push_back(s);	//6	//x var

	SetFunction(&Production::I_MesonThree_t);
	double Ret = Inte::BooleIntegration(this);

	F_var.pop_back();
	SetFunction(&Production::I_MesonThree_s);
	std::cout << "\n\n";

	return Ret;
}

double Production::I_MesonThree_t(const double t)
{
	//aliases for clarity
	double &x = F_var.at(0);
	double &y = F_var.at(1);
	double &z = F_var.at(2);
	//const double &cos0 = F_var.at(3);
	double &L0 = F_var.at(3);
	double &L_ = F_var.at(4);

	double &s_ = F_var.at(5);
	double t_ = t;
	double fc = Limit(s_, t_, x, y, z);

	//double cos0 = theta < 0 ? 0.0 : cos(theta);	//theta < 0 means integration over theta
	//double fc = theta < 0 ? 2.0 : 1.0;              //so an overall factor of 2

	std::cout << s_ << "\t" << t_ << "\t" << fc * M2_MesonThree(s_, t_, x, y, z, L_, L0) << std::endl;
	return fc * M2_MesonThree(s_, t_, x, y, z, L_, L0);
}

double Production::I_MesonThree_D(const double *v)
{
	//aliases for clarity
	double &x = F_var.at(0);
	double &y = F_var.at(1);
	double &z = F_var.at(2);
	//const double &cos0 = F_var.at(3);
	double &L0 = F_var.at(3);
	double &L_ = F_var.at(4);

	double s_ = v[0];
	double t_ = v[1];
	double fc = Limit(s_, t_, x, y, z);

	//double cos0 = theta < 0 ? 0.0 : cos(theta);	//theta < 0 means integration over theta
	//double fc = theta < 0 ? 2.0 : 1.0;              //so an overall factor of 2

	return fc * M2_MesonThree(s_, t_, x, y, z, L_, L0);
}

void Production::Reset()
{
	//heavy neutrino f
	fMuonE	= -1.0;
	fMuonM	= -1.0;
	fTauEE	= -1.0;
	fTauET	= -1.0;
	fTauMM	= -1.0;
	fTauMT	= -1.0;
	fPionE	= -1.0;
	fPionM	= -1.0;
	fKaonE	= -1.0;
	fKaonM	= -1.0;
	fCharmE	= -1.0;
	fCharmM	= -1.0;
	fCharmT	= -1.0;
	fKaon0E	= -1.0;
	fKaon0M	= -1.0;
	fKaonCE	= -1.0;
	fKaonCM	= -1.0;

	//light neutrino l
	lMuonE	= -1.0;
	lMuonM	= -1.0;
	lTauEE	= -1.0;
	lTauET	= -1.0;
	lTauMM	= -1.0;
	lTauMT	= -1.0;
	lPionE	= -1.0;
	lPionM	= -1.0;
	lKaonE	= -1.0;
	lKaonM	= -1.0;
	lCharmE	= -1.0;
	lCharmM	= -1.0;
	lCharmT	= -1.0;
	lKaon0E	= -1.0;
	lKaon0M	= -1.0;
	lKaonCE	= -1.0;
	lKaonCM	= -1.0;
}

void Production::SetFunction(double (Production::*FF)(const double))
{
	double (Amplitude::*Function)(const double) = 
		static_cast<double (Amplitude::*)(const double)>(FF); // ok!
	Amplitude::SetFunction(Function);
}

void Production::SetFunction_D(double (Production::*FF)(const double*))
{
	double (Amplitude::*Function)(const double*) = 
		static_cast<double (Amplitude::*)(const double*)>(FF); // ok!
	Amplitude::SetFunction_D(Function);
}
