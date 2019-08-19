#include "Production.h"

Production::Production()
{
	Channel_prev = _undefined;
	Reset();
}

Amplitude::Channel Production::FindChannel(std::string Name)
{
	if (chMap.size() == 0)
		LoadMap();

	std::map<Amplitude::Channel, std::string>::iterator it = chMap.begin();
	for (std::advance(it, 26); it != chMap.end(); ++it)
		if (it->second == Name)
			return it->first;

	return _undefined;
}

std::string Production::FindChannel(Amplitude::Channel Name)
{
	return chMap[Name];
}

//std::vector<std::string> Production::ListChannels()
std::vector<Amplitude::Channel> Production::ListChannels()
{
	if (chMap.size() == 0)
		LoadMap();

	std::vector<Amplitude::Channel> vName;

	std::map<Amplitude::Channel, std::string>::iterator it = chMap.begin();
	for (std::advance(it, 26); it != chMap.end(); ++it)
		vName.push_back(it->first);

	return vName;
}

bool Production::IsAllowed(Channel Name)
{
	return (MassThreshold(Name) >= MassN());
}

double Production::MassThreshold(Channel Name)
{
	if (Channel_prev != Name)
	{
		LoadMass(Name);
		Channel_prev = Name;
	}

	double Limit = vMass.at(0);
	for (int i = 1; i < vMass.size(); ++i)
		Limit -= vMass.at(i);

	return Limit;
}

double Production::Gamma(Channel Name, bool Unitary)
{
	double Result = 0.0;

	double Mix[3] = {Ue(), Um(), Ut()};
	if (Unitary)
	{
		SetUe(1.0);
		SetUm(1.0);
		SetUt(1.0);
	}

	IsChanged();
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
		case _TauPI:
			Result = TauPI();
			break;
		case _Tau2PI:
			Result = Tau2PI();
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

	if (Unitary)
	{
		SetUe(Mix[0]);
		SetUm(Mix[1]);
		SetUt(Mix[2]);
	}

	if (Result < 1e-27)
		Result = 0.0;
	return (Helicity() ? 1.0 : 2.0) * Result;
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
		double Gamma0 = Gamma(Name);

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
	if (fMuonE < 0)// || IsChanged())
		fMuonE = IsAllowed(_MuonE) ? 
			AntiLeptonNeutrinoDecay(M_Muon, M_Electron, M_Neutrino) : 0.0;

	return fMuonE * Ue(2);
}

double Production::MuonM()
{
	if (fMuonM < 0)// || IsChanged())
		fMuonM = IsAllowed(_MuonM) ? 
			LeptonNeutrinoDecay(M_Muon, M_Electron, M_Neutrino) : 0.0;

	return fMuonM * Um(2);
}

double Production::TauEE()
{

	if (fTauEE < 0)// || IsChanged())
		fTauEE = IsAllowed(_TauEE) ?
			AntiLeptonNeutrinoDecay(M_Tau, M_Electron, M_Neutrino) : 0.0;

	return fTauEE * Ue(2);
}

double Production::TauET()
{

	if (fTauET < 0)// || IsChanged())
		fTauET = IsAllowed(_TauET) ? 
			LeptonNeutrinoDecay(M_Tau, M_Electron, M_Neutrino) : 0.0;

	return fTauET * Ut(2);
}

double Production::TauMM()
{

	if (fTauMM < 0)// || IsChanged())
		fTauMM = IsAllowed(_TauMM) ? 
			AntiLeptonNeutrinoDecay(M_Tau, M_Muon, M_Neutrino) : 0.0;

	return fTauMM * Um(2);
}

double Production::TauMT()
{

	if (fTauMT < 0)// || IsChanged())
		fTauMT = IsAllowed(_TauMT) ? 
			LeptonNeutrinoDecay(M_Tau, M_Muon, M_Neutrino) : 0.0;

	return fTauMT * Ut(2);
}

double Production::TauPI()
{

	if (fTauPI < 0)// || IsChanged())
		fTauPI = IsAllowed(_TauPI) ? 
			pow(Const::U_ud * Const::DPion, 2) * LeptonTwoDecay(M_Tau, M_Pion) : 0.0;

	return fTauPI * Ut(2);
}

double Production::Tau2PI()
{
	if (fTau2PI < 0)// || IsChanged())
		fTau2PI = IsAllowed(_Tau2PI) ? 
			pow(Const::U_ud * Const::DPion, 2) * LeptonThreeDecay(M_Tau, M_Pion, M_Pion0) : 0.0;

	return fTau2PI * Ut(2);
}

double Production::PionE()
{
	if (fPionE < 0 || IsChanged())
		fPionE = IsAllowed(_PionE) ? 
			pow(Const::U_ud * Const::DPion, 2) * MesonTwoDecay(M_Pion, M_Electron) : 0.0;

	return fPionE * Ue(2);
}

double Production::PionM()
{

	if (fPionM < 0)// || IsChanged())
		fPionM = IsAllowed(_PionM) ? 
			pow(Const::U_ud * Const::DPion, 2) * MesonTwoDecay(M_Pion, M_Muon) : 0.0;

	return fPionM * Um(2);
}

double Production::KaonE()
{

	if (fKaonE < 0)// || IsChanged())
		fKaonE = IsAllowed(_KaonE) ? 
			pow(Const::U_us * Const::DKaon, 2) * MesonTwoDecay(M_Kaon, M_Electron) : 0.0;

	return fKaonE * Ue(2);
}

double Production::KaonM()
{

	if (fKaonM < 0)// || IsChanged())
		fKaonM = IsAllowed(_KaonM) ? 
			pow(Const::U_us * Const::DKaon, 2) * MesonTwoDecay(M_Kaon, M_Muon) : 0.0;

	return fKaonM * Um(2);
}

double Production::CharmE()
{

	if (fCharmE < 0)// || IsChanged())
		fCharmE = IsAllowed(_CharmE) ? 
			pow(Const::U_cs * Const::DCharm, 2) * MesonTwoDecay(M_CharmS, M_Electron) : 0.0;

	return fCharmE * Ue(2);
}

double Production::CharmM()
{

	if (fCharmM < 0)// || IsChanged())
		fCharmM = IsAllowed(_CharmM) ? 
			pow(Const::U_cs * Const::DCharm, 2) * MesonTwoDecay(M_CharmS, M_Muon) : 0.0;

	return fCharmM * Um(2);
}

double Production::CharmT()
{

	if (fCharmT < 0)// || IsChanged())
		fCharmT = IsAllowed(_CharmT) ? 
			pow(Const::U_cs * Const::DCharm, 2) * MesonTwoDecay(M_CharmS, M_Tau) : 0.0;

	return fCharmT * Ut(2);
}

double Production::Kaon0E()
{

	if (fKaon0E < 0)// || IsChanged())
		fKaon0E = IsAllowed(_Kaon0E) ?
			pow(Const::U_us * Const::KaPi, 2) * 
			//MesonThreeDecay(M_Kaon0, M_Pion, M_Electron, Const::K0L_, Const::K0L0) : 0.0;
			MesonThreeDecay(M_Kaon0, M_Pion, M_Electron, 0, 0) : 0.0;

	return fKaon0E * Ue(2);
}

double Production::Kaon0M()
{

	if (fKaon0M < 0)// || IsChanged())
		fKaon0M = IsAllowed(_Kaon0M) ? 
			pow(Const::U_us * Const::KaPi, 2) * 
			MesonThreeDecay(M_Kaon0, M_Pion, M_Muon, Const::K0L_, Const::K0L0) : 0.0;

	return fKaon0M * Um(2);
}

double Production::KaonCE()
{
	if (fKaonCE < 0)// || IsChanged())
		fKaonCE = IsAllowed(_KaonCE) ? 
			pow(Const::U_us * Const::KaPi, 2) / 2.0 * 
			MesonThreeDecay(M_Kaon, M_Pion0, M_Electron, Const::KCL_, Const::KCL0) : 0.0;

	return fKaonCE * Ue(2);
}

double Production::KaonCM()
{
	if (fKaonCM < 0)// || IsChanged())
		fKaonCM = IsAllowed(_KaonCM) ? 
			pow(Const::U_us * Const::KaPi, 2) / 2.0 * 
			MesonThreeDecay(M_Kaon, M_Pion0, M_Muon, Const::KCL_, Const::KCL0) : 0.0;

	return fKaonCM * Um(2);
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

	return I_LeptonNeutrino(dMn2, dML2, dMN2);
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
	return num::BooleIntegration(this); 
}

double Production::I_LeptonNeutrino_u(double u)	//fixing one variable
{
	//shouldnt use alias for clarity
	double x = F_var.at(0);	//light neutrino	u
	double y = F_var.at(1);	//lepton
	double z = F_var.at(2);	//heavy neutrino
	//const double &cos0 = F_var.at(3);

	double u_ = u;
	double fc = Limit(u_, y, z, x);

	//double cos0 = theta < 0 ? 0.0 : cos(theta);	//theta < 0 means integration over theta
	//double fc = theta < 0 ? 2.0 : 1.0;		//so an overall factor of 2

	double M2 = fc * M2_LeptonNeutrino(u_, x, y, z);
	return dGammad2_3B(M2);
}


//production from antilepton -> neutrino is from opposite leptonic line
double Production::AntiLeptonNeutrinoDecay(double M_Lepton0, double M_Lepton, double M_Neut)
{
	SetMass(M_Lepton0);
	double dML2 = M_Lepton*M_Lepton/Mass(2);
	double dMn2 = M_Neut*M_Neut/Mass(2);
	double dMN2 = MassN(2)/Mass(2);

	return I_AntiLeptonNeutrino(dMn2, dML2, dMN2);
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
	return num::BooleIntegration(this);
}

double Production::I_AntiLeptonNeutrino_s(double s)	//the term is written for a neutrino production
{								//therefore with heliciies inverted
	//aliases for clarity
	double x = F_var.at(0);
	double y = F_var.at(1);
	double z = F_var.at(2);
	//const double &cos0 = F_var.at(5);

	//create S var
	double s_ = s;
	double fc = Limit(s_, x, y, z);

	//double cos0 = theta < 0 ? 0.0 : cos(theta);	//theta < 0 means integration over theta
	//double fc = theta < 0 ? 2.0 : 1.0;		//so an overall factor of 2

	double M2 = fc * M2_AntiLeptonNeutrino(s_, x, y, z);
	return dGammad2_3B(M2);
}

double Production::LeptonTwoDecay(double M_Lepton, double M_Meson)
{
	SetMass(M_Lepton);
	double dMN2 = MassN(2)/Mass(2);
	double dMM2 = M_Meson*M_Meson/Mass(2);

	return I_LeptonTwo(dMN2, dMM2);
}

double Production::I_LeptonTwo(double x, double y)
{
	double M2 = M2_LeptonTwo(x, y);
	return dGammad0_2B(M2, x, y);
}

double Production::LeptonThreeDecay(double M_Lepton, double M_Meson0, double M_Meson)
{
	SetMass(M_Lepton);
	double dMN2 = MassN(2)/Mass(2);
	double dMM2 = M_Meson*M_Meson/Mass(2);
	double dM02 = M_Meson0*M_Meson0/Mass(2);

	return I_LeptonThree(dMN2, dMM2, dM02);
}

double Production::I_LeptonThree(double x, double y, double z)
{
	F_var.clear();

	F_var.push_back(x);	//0	//c2
	F_var.push_back(y);	//1	//b2
	F_var.push_back(z);	//2	//a2

	SetFunction(&Production::I_LeptonThree_s);
	return num::BooleIntegration(this); 		//switch to Vega
}

double Production::I_LeptonThree_s(double s)
{
	double x = F_var.at(0);
	double y = F_var.at(1);
	double z = F_var.at(2);
	double s_ = s;
	double fc = Limit(s_, x, y, z);

	double M2 = fc * M2_LeptonThree(x, y, z);
	return dGammad2_3B(M2);
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

/////from 1805.08567
double Production::MesonThreeDecay2(double M_Meson0, double M_Meson, double M_Lepton, double L_, double L0)	//decay constant not important
{
	SetMass(M_Meson0);
	double x = M_Meson*M_Meson/Mass(2);
	double y = M_Lepton*M_Lepton/Mass(2);
	double z = MassN(2)/Mass(2);	

	F_var.clear();

	F_var.push_back(x);	//0	//c2
	F_var.push_back(y);	//1	//b2
	F_var.push_back(z);	//2	//a2
	//F_var.push_back(cos0);	//3	//theta
	F_var.push_back(L_);	//4	//linear dep for decay constant
	F_var.push_back(L0);	//5	//linear dep for decay constant

	//int Trial, Fail;
	//double Error, Chi2Prob;

	SetFunction(&Production::I_MesonThree_2);
	double Int = num::BooleIntegration(this); 		//switch to Vega
	return Const::GF2 * Mass(5) / (64.0 * Const::pi3) * pow(Const::U_us, 2) * Int;
}
double Production::I_MesonThree_2(double s)
{
	//aliases for clarity
	double z = F_var.at(0);	//h'
	double y = F_var.at(1);	//l
	double x = F_var.at(2);	//N
	//const double &cos0 = F_var.at(3);
	double L0 = F_var.at(3);
	double L_ = F_var.at(4);

	double sInf = x + y + 2*sqrt(x*y);
	double sSup = 1 + z - 2*sqrt(z);
	double fc = sSup - sInf;
	double s_ = sInf + fc*s;

	double G = s_ * (x+y) - pow((x-y), 2);
	double L = SqrtKallen(1, z, s_)*SqrtKallen(s_, x, y);

	//double fp = Const::KaPi * (1 - L_ * s_ / x);
	//double f0 = Const::KaPi * (1 - L0 * s_ / x);
	double fp = Const::KaPi;
	double f0 = Const::KaPi;

	double I1 = fp*fp * L * L * L			/ (3 * s_*s_*s_);
	double I2 = fp*fp * L * G * Kallen(1, z, s_)	/ (2 * s_*s_*s_);
	double I3 = f0*f0 * L * G * (1-z)*(1-z)		/ (2 * s_*s_*s_);
	if (s_ == 0)
		return 0.0;

	return fc * (I1 + I2 + I3);
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

	//int Trial, Fail;
	//double Error, Chi2Prob;

	SetFunction(&Production::I_MesonThree_s);
	return num::BooleIntegration(this); 		//switch to Vega
}

double Production::I_MesonThree_s(double s)	//fixing one variable
{
	F_var.push_back(s);	//6	//x var

	SetFunction(&Production::I_MesonThree_t);
	double Ret = num::BooleIntegration(this);

	F_var.pop_back();
	SetFunction(&Production::I_MesonThree_s);

	return Ret;
}

double Production::I_MesonThree_t(double t)
{
	//aliases for clarity
	double x = F_var.at(0);
	double y = F_var.at(1);
	double z = F_var.at(2);
	//const double &cos0 = F_var.at(3);
	double L0 = F_var.at(3);
	double L_ = F_var.at(4);

	double s_ = F_var.at(5);
	double t_ = t;
	double fc = Limit(s_, t_, x, y, z);

	//double cos0 = theta < 0 ? 0.0 : cos(theta);	//theta < 0 means integration over theta
	//double fc = theta < 0 ? 2.0 : 1.0;              //so an overall factor of 2

	double M2 = fc * M2_MesonThree(s_, t_, x, y, z, L_, L0);
	return dGammad2_3B(M2);
}

double Production::I_MesonThree_D(double *v)
{
	//aliases for clarity
	double x = F_var.at(0);
	double y = F_var.at(1);
	double z = F_var.at(2);
	//const double &cos0 = F_var.at(3);
	double L0 = F_var.at(3);
	double L_ = F_var.at(4);

	double s_ = v[0];
	double t_ = v[1];
	double fc = Limit(s_, t_, x, y, z);

	//double cos0 = theta < 0 ? 0.0 : cos(theta);	//theta < 0 means integration over theta
	//double fc = theta < 0 ? 2.0 : 1.0;              //so an overall factor of 2

	return fc * M2_MesonThree(s_, t_, x, y, z, L_, L0);
}

void Production::Reset()
{
	fMuonE	= -1.0;
	fMuonM	= -1.0;
	fTauEE	= -1.0;
	fTauET	= -1.0;
	fTauMM	= -1.0;
	fTauMT	= -1.0;
	fTauPI  = -1.0;
	fTau2PI = -1.0;
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
}

void Production::SetFunction(double (Production::*FF)(double))
{
	double (Amplitude::*Function)(double) = 
		static_cast<double (Amplitude::*)(double)>(FF); // ok!
	Amplitude::SetFunction(Function);
}

void Production::SetFunction_D(double (Production::*FF)(double*))
{
	double (Amplitude::*Function)(double*) = 
		static_cast<double (Amplitude::*)(double*)>(FF); // ok!
	Amplitude::SetFunction_D(Function);
}
