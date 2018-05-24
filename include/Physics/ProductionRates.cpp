#include "ProductionRates.h"

ProductionRates::ProductionRates(double Mass, double Ue, double Um, double Ut) :
	M_Neutrino(0.0),
	M_Photon(0.0),
	M_Electron(Const::fMElectron),
	M_Muon(Const::fMMuon),
        M_Tau(Const::fMTau),
	M_Pion(Const::fMPion),
	M_Pion0(Const::fMPion0),
	M_Kaon(Const::fMKaon),
	M_Kaon0(Const::fMKaon0),
        M_CharmS(Const::fMDs),
	fKaon(Const::fVusFKaon),
	fLambda1(Const::fLambdaPlus),
	fLambda0(Const::fLambdaZero)
{
	M_Sterile_prev = -1.0;

	SetMass(Mass);
	SetUe(Ue);
	SetUm(Um);
	SetUt(Ut);

	IsElectron = false;
	IsMuon = false;
	IsTau = false;

	InitConst();
}

bool ProductionRates::IsAllowed(Channel Name)
{
	double Limit = 0.0;

	switch(Name)
	{
		case _ALL:
		case _MuonE:
		case _MuonM:
			Limit = M_Muon - M_Electron - M_Neutrino;
			break;
		case _TauEE:
		case _TauET:
			Limit = M_Tau - M_Electron - M_Neutrino;
			break;
		case _TauMM:
		case _TauMT:
			Limit = M_Tau - M_Muon;
			break;
		case _PionE:
			Limit = M_Pion - M_Electron;
			break;
		case _PionM:
			Limit = M_Pion - M_Muon;
			break;
		case _KaonE:
			Limit = M_Kaon - M_Electron;
			break;
		case _KaonM:
			Limit = M_KaonM - M_Muon;
			break;
		case _CharmE:
			Limit = M_CharmS - M_Electron;
			break;
		case _CharmM:
			Limit = M_CharmS - M_Muon;
			break;
		case _CharmT:
			Limit = M_CharmS - M_Tau;
			break;
		case _Kaon0E:
			Limit = M_Kaon0 - M_Pion - M_Electron;
			break;
		case _Kaon0M:
			Limit = M_Kaon0 - M_Pion - M_Muon;
			break;
		case _KaonCE:
			Limit = M_Kaon - M_Pion0 - M_Electron;
			break;
		case _KaonCM:
			Limit = M_Kaon - M_Pion0 - M_Muon;
			break;
		default:
			Limit = -1.0;
			break;
	}

	return (Limit >= GetMass());
}

double ProductionRates::Gamma(Channel Name)
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

double ProductionRates::Scale(Channel Name)
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

double ProductionRates::Total()
{
	return (MuonE() + MuonM() + TauEE() + TauET() + TauME() + TauMT() +
		PionE() + PionM() + KaonE() + KaonM() + CharmE() + CharmM() + CharmT() +
		Kaon0E() + Kaon0M() + KaonCE() + KaonCM());
}

double ProductionRates::MuonE()
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

double ProductionRates::MuonM()
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

double ProductionRates::TauEE()
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

double ProductionRates::TauET()
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

double ProductionRates::TauMM()
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

double ProductionRates::TauMT()
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

double ProductionRates::PionE()
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

double ProductionRates::PionM()
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

double ProductionRates::KaonE()
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

double ProductionRates::KaonM()
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

double ProductionRates::CharmE()
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

double ProductionRates::CharmM()
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

double ProductionRates::CharmT()
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

/////////////////
//Generic decay//
/////////////////
//

double ProductionRates::LeptonNeutrinoDecay(double M_LeptonA, double M_LeptonB, double M_Neutrino)
{
	double dMB2 = M_LeptonB*M_LeptonB/M_LeptonA/M_LeptonA;
	double dMn2 = M_Neutrino*M_Neutrino/M_LeptonA/M_LeptonA;
	double dMN2 = GetMass()*GetMass()/M_LeptonA/M_LeptonA;

	return Const::fGF2 / (16.0 * Const::fPi3) *
		I_LeptonNeutrino(dMB2, dMn2, dMN2);
}
//						c	  b	    a
double ProductionRates::I_LeptonNeutrino(double x, double y, double z, double theta)	//integrate first in y and then in x
{
	I_var.clear();

	I_var.push_back(x);	//0	//x
	I_var.push_back(y);	//1	//y
	I_var.push_back(z);	//2	//z
	I_var.push_back(cos0);	//3	//theta

	SetFunction(&I_LeptonNeutrino_s);
	return Inte::BooleIntegration(this); 
}

double ProductionRates::I_LeptonNeutrino_s(double s)	//fixing one variable
{
	I_var.push_back(s)	//4
	SetFunction(&I_LeptonNeutrino_t);
	double Ret = Inte::BooleIntegration(this);

	I_var.pop_back();
	SetFunction(&I_LeptonNeutrino_s);

	return Ret;
}

double ProductionRates::I_LeptonNeutrino_t(double t)
{
	//aliases for clarity
	const double &x = I_var.at(0);
	const double &y = I_var.at(1);
	const double &z = I_var.at(2);
	const double &cos0 = I_var.at(3);

	const double &s = I_var(4);

	const double &X = I_var.at(8);
	const double &yInf = I_var.at(9);
	const double &ySup = I_var.at(10);

	//create X & Y var
	double xInf = 2*sqrt(z);
	double xSup = 1 + z - x - y - 2*sqrt(x*y);
	double X = 1 + z - (xInf + (xSup - xInf) * s);

	double Y = yInf + (ySup - yInf) * t;
	double yInf = ( (2 - X)*(X + y - x) - sqrt(Kine::Kallen(X, y, x)*Kine::Kallen(1, X, z)) ) / (2*X);
	double ySup = ( (2 - X)*(X + y - x) + sqrt(Kine::Kallen(X, y, x)*Kine::Kallen(1, X, z)) ) / (2*X);

	double cos0 = theta < 0 ? 0.0 : cos(theta);	//theta < 0 means integration over theta
	double fc = theta < 0 ? 2.0 : 1.0;		//so an overall factor of 2

	return fc * (xSup - xInf) * (ySup - yInf) * d_LeptonNeutrino(X, Y, x, y, z, cos0);
}

double ProductionRates::d_LeptonNeutrino(double X, double Y, double x, double y, double z, double cos0);
{
	double kX = sqrt(X*X - 4*z);
	double kY = sqrt(Y*Y - 4*y)*cos0;

	if (GetHelicity() < 0)
		return 1 + y - z - x - Y - (2 - X - Y - kX - kY)*(X - kX)/4.0;
	else if (GetHelicity() > 0)
		return  (2 - X - Y - kX - kY)*(X - kX)/4.0;
}

double ProductionRates::LeptonAntineutrinoDecay(double M_LeptonA, double M_LeptonB, double M_Neutrino)
{
	double dMB2 = M_LeptonB*M_LeptonB/M_LeptonA/M_LeptonA;
	double dMn2 = M_Neutrino*M_Neutrino/M_LeptonA/M_LeptonA;
	double dMN2 = GetMass()*GetMass()/M_LeptonA/M_LeptonA;

	return Const::fGF2 / (16.0 * Const::fPi3) *
		I_LeptonNeutrino(dMB2, dMn2, dMN2);
}
						//  c	      b		a
double ProductionRates::I_LeptonAntineutrino(double x, double y, double z, double theta)	//integrate first in y and then in x
{
	I_var.clear();

	I_var.push_back(x);	//0	//x
	I_var.push_back(y);	//1	//y
	I_var.push_back(z);	//2	//z
	I_var.push_back(cos0);	//3	//theta

	SetFunction(&I_LeptonAntineutrino_s);
	return Inte::BooleIntegration(this); 
}

double ProductionRates::I_LeptonAntineutrino_s(double s)	//the term is written for a neutrino production
{								//therefore with heliciies inverted
	//aliases for clarity
	const double &x = I_var.at(2);
	const double &y = I_var.at(3);
	const double &z = I_var.at(4);
	const double &cos0 = I_var.at(5);

	//create S var
	double S = sInf + (sSup - sInf) * s;
	double sInf = x*x + y*y + 2*sqrt(x*y);
	double sSup = 1 - 2*sqrt(z) + z;

	double cos0 = theta < 0 ? 0.0 : cos(theta);	//theta < 0 means integration over theta
	double fc = theta < 0 ? 2.0 : 1.0;		//so an overall factor of 2

	return fc * (sSup - sInf) * d_LeptonAntineutrino_s(S, x, y, z, cos0);
}

double ProductionRates::d_LeptonAntineutrino_s(double S, double x, double y, double z, double cos0)
{
	double Lambda0 = sqrt(Kine::Kallen(1, S, z));
	double Lambda1 = sqrt(Kine::Kallen(S, y, x));

	if (GetHelicity() < 0)
		return  (S - y - x) * (1 + z - S + Lambda0*cos0) * Lambda0 * Lambda1 / S;
	else if (GetHelicity() > 0)
		return (S - y - x) * (1 + z - S - Lambda0*cos0) * Lambda0 * Lambda1 / S;
}

double ProductionRates::LeptonMesonDecay(double M_Lepton, double M_Meson, double fDecay2)
{
	double dMN2 = GetMass()*GetMass()/M_Lepton/M_Lepton;
	double dMM2 = M_Meson*M_Meson/M_Lepton/M_Lepton;

	return Const::fGF2 * fDecay2 / (16.0 * Const::fPi) *
		pow(M_Lepton, 3) * I_LeptonMeson(dMN2, dMM2);
}

double ProductionRates::I_LeptonMeson(double x, double y)	//y is the meson
{
	double Lambda = sqrt(Kine::Kallen(1, x, y));
	if (GetHelicity() < 0)
		return Lambda * (pow(1 - x, 2) - y * (1 + x) - (x - 1) * Lambda);
	else if (GetHelicity() > 0)
		return Lambda * (pow(1 - x, 2) - y * (1 + x) + (x - 1) * Lambda);
}

double ProductionRates::MesonTwoDecay(double M_Meson, double M_Lepton, double fDecay2)
{
	double dMN2 = GetMass()*GetMass()/M_Meson/M_Meson;
	double dML2 = M_Lepton*M_Lepton/M_Meson/M_Meson;

	return Const::fGF2 * fDecay2 / (16.0 * Const::fPi)
		pow(M_Meson, 3) * I_MesonTwo(dMN2, dML2);
}

double ProductionRates::I_MesonTwo(double x, double y)	//symetric in x and y
{
	double Lambda = sqrt(Kine::Kallen(1, x, y));
	if (GetHelicity() < 0)
		return Lambda * (x + y - pow(x - y, 2) - (x - y) * Lambda);
	else if (GetHelicity() > 0)
		return Lambda * (x + y - pow(x - y, 2) + (x - y) * Lambda);
}

double ProductionRates::MesonThreeDecay(double M_Meson0, double M_Meson1, double M_Lepton, double vCKM, double Decay2, double L_, double L0)
{
	double dMM2 = M_Meson1*M_Meson1/M_Meson0/M_Meson0;	// c
	double dML2 = M_Lepton*M_Lepton/M_Meson0/M_Meson0;	// b
	double dMN2 = GetMass()*GetMass()/M_Meson0/M_Meson0;	// a

	return Const::fGF2 * pow(vCKM, 2) * fDecay2 / (64.0 * Const::fPi3) *
		pow(M_Meson0, 5) * I_MesonThree(dMM2, dML2, dMN2, L_, L0);
}

double ProductionRates::I_MesonThree(double x, double y, double z, double L_, double L0)	//no angle dependence
{
	I_var.clear();

	I_var.push_back(x);	//0	//c2
	I_var.push_back(y);	//1	//b2
	I_var.push_back(z);	//2	//a2
	I_var.push_back(cos0);	//3	//theta
	I_var.push_back(L_);	//4	//linear dep for decay constant
	I_var.push_back(L0);	//5	//linear dep for decay constant

	SetFunction(&I_MesonThree_s);
	return Inte::BooleIntegration(this); 
}

double ProductionRates::I_MesonThree_s(double s)	//fixing one variable
{
	I_var.push_back(s);	//6	//x var

	SetFunction(&I_MesonThree_t);
	double Ret = Inte::BooleIntegration(this);

	I_var.pop_back();
	SetFunction(&I_MesonThree_s);

	return Ret;
}
double ProductionRates::I_MesonThree_t(double t)
{
	//aliases for clarity
	const double &x = I_var.at(0);
	const double &y = I_var.at(1);
	const double &z = I_var.at(2);
	const double &cos0 = I_var.at(3);
	const double &L0 = I_var.at(4);
	const double &L_ = I_var.at(5);

	const double &s = I_var.at(6);

	double xInf = 2*sqrt(z);
	double xSup = 1 + z - x - y - 2*sqrt(x*y);
	double X = 1 + z - (xInf + (xSup - xInf) * s);

	double yInf = ( (2 - X)*(X + y - x) - sqrt(Kine::Kallen(X, y, x)*Kine::Kallen(1, X, z)) ) / (2*X);
	double ySup = ( (2 - X)*(X + y - x) + sqrt(Kine::Kallen(X, y, x)*Kine::Kallen(1, X, z)) ) / (2*X);
	double Y = yInf + (ySup - yInf) * t;

	double cos0 = theta < 0 ? 0.0 : cos(theta);	//theta < 0 means integration over theta
	double fc = theta < 0 ? 2.0 : 1.0;              //so an overall factor of 2

	return fc * (xSup - xInf) * (ySup - yInf) * d_MesonThree(X, Y, x, y, z, cos0, L_, L0);
}

double ProductionRates::d_MesonThree(double S, double Y, double x, double y, double cos0, double L_, double L0)
{
	double F = 2 * - 2*L_ * (X + Y - 1 + x) / x;
	double G = F/2.0 - (L0 - L_) * (1 - x) / x;
	double CompS = 1 + z - S - sqrt(Kine::Kallen(1, z, S));
	double CompY = Y - sqrt(Y*Y - 4*y)*cos0;

	if (GetHelicity() < 0)
	{
		double A = (1 + z - S)*Y - S + Y - y + x - CompS*CompY / 4.0;
		double B = z *(S -Y + y - x) + 2*z*y - (z - y) * CompS*CompY;
		double C = z*Y + 2*y*(1 + z - S) + CompS / 2.0 * (S - Y - y - x + (1 + z - S)*CompY / 2.0);

		return F*F * A + G*G * B - F*G * C;
	}
	else if (GetHelicity() > 0)
	{
		double A = CompS*CompY / 4.0;
		double B = y *(S -Y + y - x) + 2*z*y + (z - y) * CompS*CompY;
		double C = z*Y - CompX / 2.0 * (1 + z - y - x - X - Y + X*CompY / 2.0 );
		double C = z * (1 + z - S) - CompS / 2.0 * (S - Y -y -x + (1 + z - S)*CompY / 2.0);

		return F*F * A + G*G * B - F*G * C;
	}
}

//Integration set up
//
void ProductionRates::SetFunction(double (ProductionRates::*FF)(double))
{
	fFunction = FF;
}

double ProductionRates::Function(double x)
{
	return (*fFunction)(x);
}

//Getter
//
double ProductionRates::GetMass()
{
	return fMass;
}

int ProductionRates::GetHelicity()
{
	return iHel;
}

//Setter
void ProductionRates::SetMass(double Mass)
{
	fMass = Mass;
}

void ProductionRates::SetHelicity(int Helicity)
{
	iHel = Helicity;
}

void ProductionRates::SetUe(double X)
{
	Ue = X;
}

void ProductionRates::SetUm(double X)
{
	Um = X;
}

void ProductionRates::SetUt(double X)
{
	Ut = X;
}

bool ProductionRates::IsChanged()
{
	bool Ret = (fabs(GetMass() - M_Sterile_prev) > 1e-9);

	M_Sterile_prev = GetMass();

	//Reset decay widths if changed
	if (Ret)
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
                fKaon0	= -1.0;
                fKaonCE	= -1.0;
                fKaonCM	= -1.0;
	}

	return Ret;
}
