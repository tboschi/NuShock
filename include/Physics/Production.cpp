#include "Production.h"

Production::Production() :
{
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

	return (Limit >= GetMass());
}

double Production::dGamma(Channel Name)
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
				Result = 0.0;
				break;
		}
		return Result;
	}
}

double Production::Scale(Channel Name)
{
	int ih = GetHelicity();
	double mass = GetMass();

	SetHelicity(0);
	SetMass(0);
	double Gamma0 = 2*Gamma(Name);

	SetHelicity(ih);
	SetMass(mass);

	return Gamma(Name)/Gamma0;
}

double Production::Total()
{
	return (MuonE() + MuonM() + TauEE() + TauET() + TauME() + TauMT() +
		PionE() + PionM() + KaonE() + KaonM() + CharmE() + CharmM() + CharmT() +
		Kaon0E() + Kaon0M() + KaonCE() + KaonCM());
}

double Production::MuonE()
{
	if (fMuonE < 0 || IsChanged())
	{
		if (IsAllowed(_MuonE))
			fMuonE = LeptonAntineutrinoDecay(M_Muon, M_Electron, M_Neutrino);
		else
			fMuonE  = 0.0;
	}

	return fMuonE * Ue*Ue;
}

double Production::MuonM()
{
	if (fMuonM < 0 || IsChanged())
	{
		if (IsAllowed(_MuonM))
			fMuonM = LeptonNeutrinoDecay(M_Muon, M_Electron, M_Neutrino);
		else
			fMuonM  = 0.0;
	}

	return fMuonM * Um*Um;
}

double Production::TauEE()
{
	if (fTauEE < 0 || IsChanged())
	{
		if (IsAllowed(_TauEE))
			fTauEE = LeptonAntineutrinoDecay(M_Tau, M_Electron, M_Neutrino);
		else
			fTauEE  = 0.0;
	}

	return fTauEE * Ue*Ue;
}

double Production::TauET()
{
	if (fTauET < 0 || IsChanged())
	{
		if (IsAllowed(_TauET))
			fTauET = LeptonNeutrinoDecay(M_Tau, M_Electron, M_Neutrino);
		else
			fTauET  = 0.0;
	}

	return fTauET * Ut*Ut;
}

double Production::TauMM()
{
	if (fTauME < 0 || IsChanged())
	{
		if (IsAllowed(_TauME))
			fTauME = LeptonAntineutrinoDecay(M_Tau, M_Muon, M_Neutrino);
		else
			fTauME  = 0.0;
	}

	return fTauME * Um*Um;
}

double Production::TauMT()
{
	if (fTauMT < 0 || IsChanged())
	{
		if (IsAllowed(_TauMT))
			fTauMT = LeptonNeutrinoDecay(M_Tau, M_Muon, M_Neutrino);
		else
			fTauMT  = 0.0;
	}

	return fTauMT * Ut*Ut;
}

double Production::PionE()
{
	if (fPionE < 0 || IsChanged())
	{
		if (IsAllowed(_PionE))
			fPionE = MesonTwoDecay(M_Pion, M_Electron, Const::fPion);
		else
			fPionE  = 0.0;
	}

	return fPionE * Ue*Ue;
}

double Production::PionM()
{
	if (fPionM < 0 || IsChanged())
	{
		if (IsAllowed(_PionM))
			fPionM = MesonTwoDecay(M_Pion, M_Muon, Const::fPion);
		else
			fPionM  = 0.0;
	}

	return fPionM * Um*Um;
}

double Production::KaonE()
{
	if (fKaonE < 0 || IsChanged())
	{
		if (IsAllowed(_KaonE))
			fKaonE = MesonTwoDecay(M_Kaon, M_Electron, Const::fKaon);
		else
			fKaonE  = 0.0;
	}

	return fKaonE * Ue*Ue;
}

double Production::KaonM()
{
	if (fKaonM < 0 || IsChanged())
	{
		if (IsAllowed(_KaonM))
			fKaonM = MesonTwoDecay(M_Kaon, M_Muon, Const::fKaon);
		else
			fKaonM  = 0.0;
	}

	return fKaonM * Um*Um;
}

double Production::CharmE()
{
	if (fCharmE < 0 || IsChanged())
	{
		if (IsAllowed(_CharmE))
			fCharmE = MesonTwoDecay(M_CharmS, M_Electron, Const::fCharm);
		else
			fCharmE  = 0.0;
	}

	return fCharmE * Ue*Ue;
}

double Production::CharmM()
{
	if (fCharmM < 0 || IsChanged())
	{
		if (IsAllowed(_CharmM))
			fCharmM = MesonTwoDecay(M_CharmS, M_Muon, Const::fCharm);
		else
			fCharmM  = 0.0;
	}

	return fCharmM * Um*Um;
}

double Production::CharmT()
{
	if (fCharmT < 0 || IsChanged())
	{
		if (IsAllowed(_CharmT))
			fCharmT = MesonTwoDecay(M_CharmS, M_Tau, Const::fCharm);
		else
			fCharmT  = 0.0;
	}

	return fCharmT * Ut*Ut;
}

double Production::Kaon0E()
{
	if (fKaon0E < 0 || IsChanged())
	{
		if (IsAllowed(_Kaon0E))
			fKaon0E = MesonThreeDecay(M_Kaon0, M_Pion, M_Electron, Const::fUus, Const::fK0L_, Const::fK0L0);
		else
			fKaon0E  = 0.0;
	}

	return fKaon0E * Ue*Ue;
}

double Production::Kaon0M()
{
	if (fKaon0M < 0 || IsChanged())
	{
		if (IsAllowed(_Kaon0M))
			fKaon0M = MesonThreeDecay(M_Kaon0, M_Pion, M_Muon, Const::fUus, Const::fK0L_, Const::fK0L0);
		else
			fKaon0M  = 0.0;
	}

	return fKaon0M * Ue*Ue;
}

double Production::KaonCE()
{
	if (fKaonCE < 0 || IsChanged())
	{
		if (IsAllowed(_KaonCE))
			fKaonCE = MesonThreeDecay(M_Kaon0, M_Pion, M_Electron, Const::fUus, Const::fK_L_, Const::fK_L0);
		else
			fKaonCE  = 0.0;
	}

	return fKaonCE * Ue*Ue;
}

double Production::KaonCM()
{
	if (fKaonCM < 0 || IsChanged())
	{
		if (IsAllowed(_KaonCM))
			fKaonCM = MesonThreeDecay(M_Kaon, M_Pion0, M_Electron, Const::fUus, Const::fK_L_, Const::fK_L0);
		else
			fKaonCM  = 0.0;
	}

	return fKaonCM * Um*Um;
}

/////////////////
//Generic decay//
/////////////////
//

double Production::LeptonNeutrinoDecay(double M_LeptonA, double M_LeptonB, double M_Neut)
{
	double dMB2 = M_LeptonB*M_LeptonB/M_LeptonA/M_LeptonA;
	double dMn2 = M_Neut*M_Neut/M_LeptonA/M_LeptonA;
	double dMN2 = GetMass()*GetMass()/M_LeptonA/M_LeptonA;

	double M2 = I_LeptonNeutrino(dMB2, dMn2, dMN2);
	return gGamma2_3B(M2);
	//return Const::fGF2 / (16.0 * Const::fPi3) *
	//	d_LeptonNeutrino(dMB2, dMn2, dMN2);
}

//						c	  b	    a
double PFoductionRates::I_LeptonNeutrino(double x, double y, double z)//, double theta)	//integrate first in y and then in x
{
	F_var.clear();

	F_var.push_back(x);	//0	//x
	F_var.push_back(y);	//1	//y
	F_var.push_back(z);	//2	//z
	//F_var.push_back(cos0);	//3	//theta

	SetFunction(&I_LeptonNeutrino_s);
	return Inte::BooleIntegration(this); 
}

double Production::I_LeptonNeutrino_s(double s)	//fixing one variable
{
	F_var.push_back(s)	//4
	SetFunction(&I_LeptonNeutrino_t);
	double Ret = Inte::BooleIntegration(this);

	F_var.pop_back();
	SetFunction(&I_LeptonNeutrino_s);

	return Ret;
}

double Production::I_LeptonNeutrino_t(double t)
{
	//aliases for clarity
	const double &x = F_var.at(0);
	const double &y = F_var.at(1);
	const double &z = F_var.at(2);
	//const double &cos0 = F_var.at(3);

	double &s = F_var.at(3);

	double fc = Limit(s, t, x, y, z);

	//double cos0 = theta < 0 ? 0.0 : cos(theta);	//theta < 0 means integration over theta
	//double fc = theta < 0 ? 2.0 : 1.0;		//so an overall factor of 2

	return fc * M2_LeptonNeutrino(x, y, z, s, t, 0);
}

double Production::LeptonAntineutrinoDecay(double M_LeptonA, double M_LeptonB, double M_Neut)
{
	double dMB2 = M_LeptonB*M_LeptonB/M_LeptonA/M_LeptonA;
	double dMn2 = M_Neut*M_Neut/M_LeptonA/M_LeptonA;
	double dMN2 = GetMass()*GetMass()/M_LeptonA/M_LeptonA;

	double M2 = I_LeptonAntineutrino(dMB2, dMn2, dMN2);
	return dGammad2_3B(M2);
}
						//  c	      b		a
double Production::I_LeptonAntineutrino(double x, double y, double z)//, double theta)	//integrate first in y and then in x
{
	F_var.clear();

	F_var.push_back(x);	//0	//x
	F_var.push_back(y);	//1	//y
	F_var.push_back(z);	//2	//z
	//F_var.push_back(cos0);	//3	//theta

	SetFunction(&I_LeptonAntineutrino_s);
	return Inte::BooleIntegration(this); 
}

double Production::I_LeptonAntineutrino_s(double s)	//the term is written for a neutrino production
{								//therefore with heliciies inverted
	//aliases for clarity
	const double &x = F_var.at(0);
	const double &y = F_var.at(1);
	const double &z = F_var.at(2);
	//const double &cos0 = F_var.at(5);

	//create S var
	double t = 0;
	Limit(s, t, x, y, z);

	//double cos0 = theta < 0 ? 0.0 : cos(theta);	//theta < 0 means integration over theta
	//double fc = theta < 0 ? 2.0 : 1.0;		//so an overall factor of 2

	return fc * M2_LeptonAntineutrino(x, y, z, s, t, 0);
}

double Production::LeptonMesonDecay(double M_Lepton, double M_Meson)
{
	double dMN2 = GetMass()*GetMass()/M_Lepton/M_Lepton;
	double dMM2 = M_Meson*M_Meson/M_Lepton/M_Lepton;

	return I_LeptonMeson(dMN2, dMM2);
}

double PRoductionRates::I_LeptonMeson(double x, double y)
{
	double M2 = M2_LeptonMeson(x, y, cos0);
	return dGamma0_2B(M2, x, y);
}

double Production::MesonTwoDecay(double M_Meson, double M_Lepton)
{
	double dMN2 = GetMass()*GetMass()/M_Meson/M_Meson;
	double dML2 = M_Lepton*M_Lepton/M_Meson/M_Meson;

	return I_MesonTwo(dMN2, dML2)
}

double Production::I_MesonTwo(double x, double y)	//symetric in x and y
{
	double M2 = M2_MesonTwo(x, y);
	return dGammad0_2B(M2, x, y);
}

double Production::MesonThreeDecay(double M_Meson0, double M_Meson1, double M_Lepton, double vCKM, double L_, double L0)	//decay constant not important
{
	double dMM2 = M_Meson1*M_Meson1/M_Meson0/M_Meson0;
	double dML2 = M_Lepton*M_Lepton/M_Meson0/M_Meson0;
	double dMN2 = GetMass()*GetMass()/M_Meson0/M_Meson0;	

	return Const::fGF2 * pow(vCKM, 2) / (64.0 * Const::fPi3) *
		pow(M_Meson0, 5) * d_MesonThree(dMM2, dML2, dMN2, L_, L0);
}

double Production::d_MesonThree(double x, double y, double z, double L_, double L0)	//no angle dependence
{
	F_var.clear();

	F_var.push_back(x);	//0	//c2
	F_var.push_back(y);	//1	//b2
	F_var.push_back(z);	//2	//a2
	F_var.push_back(cos0);	//3	//theta
	F_var.push_back(L_);	//4	//linear dep for decay constant
	F_var.push_back(L0);	//5	//linear dep for decay constant

	SetFunction(&d_MesonThree_s);
	return Inte::BooleIntegration(this); 
}

double Production::d_MesonThree_s(double s)	//fixing one variable
{
	F_var.push_back(s);	//6	//x var

	SetFunction(&d_MesonThree_t);
	double Ret = Inte::BooleIntegration(this);

	F_var.pop_back();
	SetFunction(&d_MesonThree_s);

	return Ret;
}
double Production::d_MesonThree_t(double t)
{
	//aliases for clarity
	const double &x = F_var.at(0);
	const double &y = F_var.at(1);
	const double &z = F_var.at(2);
	const double &cos0 = F_var.at(3);
	const double &L0 = F_var.at(4);
	const double &L_ = F_var.at(5);

	const double &s = F_var.at(6);

	double sInf = x*x + y*y + 2*sqrt(x*y);
	double sSup = 1 - 2*sqrt(z) + z;
	double S = sInf + (sSup - sInf) * s;

	double tInf = ( (2 - S)*(S + y - x) - sqrt(Kine::Kallen(S, y, x)*Kine::Kallen(1, S, z)) ) / (2*S);
	double tSup = ( (2 - S)*(S + y - x) + sqrt(Kine::Kallen(S, y, x)*Kine::Kallen(1, S, z)) ) / (2*S);
	double T = yInf + (ySup - yInf) * t;

	double cos0 = theta < 0 ? 0.0 : cos(theta);	//theta < 0 means integration over theta
	double fc = theta < 0 ? 2.0 : 1.0;              //so an overall factor of 2

	return fc * (sSup - sInf) * (tSup - tInf) * d_MesonThree(1 + z - S, T, x, y, z, cos0, L_, L0);
}

double Production::d_MesonThree(double Z, double Y, double x, double y, double cos0, double L_, double L0)
{
	double F = 2 * - 2*L_ * (X + Y - 1 + x) / x;
	double G = F/2.0 - (L0 - L_) * (1 - x) / x;
	double CompZ = Z - sqrt(Z*Z - 4*z);
	double CompY = Y - sqrt(Y*Y - 4*y)*cos0;

	if (GetHelicity() < 0)
	{
		double A = X*Y + Z + Y - z - y - 1 + x - CompZ*CompY / 4.0;
		double B = z * (1 - Z + z - Y + y - x) + 2*z*y - (z - y) * CompZ*CompY;
		double C = z*Y + 2*y*Z + CompZ / 2.0 * (1 - Z + z -Y - y - x + Z * CompY / 2.0);

		return F*F * A + G*G * B - F*G * C;
	}
	else if (GetHelicity() > 0)
	{
		double A = CompS*CompY / 4.0;
		double B = y * (1 - Z + z - Y + y - x) + 2*z*y + (z - y) * CompZ*CompY;
		double C = z*Y + - CompZ / 2.0 * (1 - Z + z -Y - y - x + Z * CompY / 2.0);

		return F*F * A + G*G * B - F*G * C;
	}
}

void Production::Reset()
{
	fMuonE	= -1.0;
	fMuonM	= -1.0;
	fTauEE	= -1.0;
	fTauET	= -1.0;
	fTauME	= -1.0;
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

//Set functions
void DecayRates::SetNvec(TLorentzVector &X)
{
	*N_vec = X;
	N_rest->SetPxPyPzE(0, 0, 0, N_vec->M());
}
