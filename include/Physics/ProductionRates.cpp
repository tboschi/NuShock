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

bool DecayRates::IsAllowed(Channel Name)
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
	SetMas(0);
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
			fMuonE = LeptonDecay(M_Muon, M_Electron, GetMass());
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
			fMuonM = LeptonDecay(M_Muon, GetMass(), M_Electron);
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
			fTauEE = LeptonDecay(M_Tau, M_Electron, GetMass());
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
			fTauET = LeptonDecay(M_Tau, GetMass(), M_Electron);
		else
			fTauET  = 0.0;
	}

	return fTauET * Ut*Ut;
}

double ProductionRates::TauME()
{
	if (fTauME < 0 || IsChanged())
	{
		if (IsAllowed(_TauME))
			fTauME = LeptonDecay(M_Tau, M_Muon, GetMass());
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
			fTauMT = LeptonDecay(M_Tau, GetMass(), M_Muon);
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

double ProductionRates::LeptonThree(double M_Lepton0, double M_Lepton1, double M_Lepton3)
{
	double dM12 = 
}

double ProductionRates::LeptonMesonDecay(double M_Lepton, double M_Meson, double fDecay2)
{
	double dMN2 = GetMass()*GetMass()/M_Lepton/M_Lepton;
	double dMM2 = M_Meson*M_Meson/M_Lepton/M_Lepton;

	return Const::fGF2 * fDecay2 / (16.0 * Const::fPi) *
		pow(M_Lepton, 3) * I_Lepton(dMN2, dMM2);
}

double ProductionRates::I_Lepton(double x, double y)	//y is the meson
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
		pow(M_Meson, 3) * I_Meson(dMN2, dML2);
}

double ProductionRates::I_Meson(double x, double y)	//symetric in x and y
{
	double Lambda = sqrt(Kine::Kallen(1, x, y));
	if (GetHelicity() < 0)
		return Lambda * (x + y - pow(x - y, 2) - (x - y) * Lambda);
	else if (GetHelicity() > 0)
		return Lambda * (x + y - pow(x - y, 2) + (x - y) * Lambda);
}

double ProductionRates::MesonThreeDecay(double M_Meson0, double M_Meson1, double M_Lepton, double vCKM, double fDecay2)
{
	double dMM2 = M_Meson1*M_Meson1/M_Meson0/M_Meson0;	// c
	double dML2 = M_Lepton*M_Lepton/M_Meson0/M_Meson0;	// b
	double dMN2 = GetMass()*GetMass()/M_Meson0/M_Meson0;	// a



	return Const::fGF2 * pow(vCKM, 2) * fDecay2 / (64.0 * Const::fPi3) *
		pow(M_Meson0, 5) * I_MesonThree(dMM2);
}

double ProductionRates::I_MesonThree(double x, double y, double z, double theta)	//no angle dependence
{
	double Inf = 2*sqrt(z);
	double Sup = 1 + z - x - y - 2*sqrt(x*y);

	I_var.clear();
	I_var.push_back(Inf);	//0
	I_var.push_back(Sup);	//1

	I_var.push_back(x);	//2	//c2
	I_var.push_back(y);	//3	//b2
	I_var.push_back(z);	//4	//a2
	I_var.push_back(theta);	//5	//theta

	SetIntegrand2(&I_MesonThree_s);
	return (Sup - Inf) * Inte::BooleIntegration(this); 
}

double ProductionRates::I_MesonThree_s(double s)	//fixing one variable
{
	//aliases for clarity
	const double &xInf = I_var.at(0);
	const double &xSup = I_var.at(1);

	const double &c2 = I_var.at(2);
	const double &b2 = I_var.at(3);
	const double &a2 = I_var.at(4);

	//create X var and Y lim
	double X = 1 + a2 - (xInf + (xSup - xInf) * s);
	double yInf = ( (2 - X)*(X + b2 - c2) - sqrt(Kine::Kallen(X, b2, c2)*Kine::Kallen(1, X, a2)) ) / (2*X);
	double ySup = ( (2 - X)*(X + b2 - c2) + sqrt(Kine::Kallen(X, b2, c2)*Kine::Kallen(1, X, a2)) ) / (2*X);

	I_var.push_back(X);	//6	//x var
	I_var.push_back(yInf);	//7	//lim inf of y
	I_var.push_back(ySup);	//8	//lim sup of y

	SetIntegrand(&I_MesonThree_t);
	double Ret = Inte::BooleIntegration(this);

	I_var.pop_back();
	I_var.pop_back();
	I_var.pop_back();
	SetIntegrand(&I_MesonThree_s);

	return (ySup - yInf) * Ret;
}
double ProductionRates::I_MesonThree_t(double t)
{
	//aliases for clarity
	const double &xInf = I_var.at(0);
	const double &xSup = I_var.at(1);

	const double &c2 = I_var.at(2);
	const double &b2 = I_var.at(3);
	const double &a2 = I_var.at(4);
	const double &theta = I_var.at(5);

	const double &L0 = I_var.at(6);
	const double &L_ = I_var.at(7);

	const double &X = I_var.at(8);
	const double &yInf = I_var.at(9);
	const double &ySup = I_var.at(10);

	//create Y var
	double Y = yInf + (ySup - yInf) * t;

	double F = 2 * - 2*L_ * (X + Y - 1 + c2) / c2;
	double G = F/2.0 - (L0 - L_) * (1 - c2) / c2;
	double CompX = X - sqrt(X*X - 4a2);
	double CompY = Y - sqrt(Y*Y - 4b2)*cos(theta);

	if (GetHelicity() < 0)
	{
		double A = X*Y - CompX*CompY/2.0 + 1 +a2 + b2 -c1 - X - Y;
		double B = b2 *( X + Y - 1 - a2 - b2 + c2 ) + 2*a2*b2 + (a2 - b2)*CompX*CompY/4;
		double C = b2 * X - CompX * ( X + Y - 1 - a2 + b2 + c2 - X*CompY / 2.0 ) /2.0;
		return F*F * A + G*G * B - F*G * C;
	}
	else if (GetHelicity() > 0)
	{
		double A = CompX*CompY/2.0;
		double B = a2 *( X + Y - 1 - a2 - b2 + c2 ) + 2*a2*b2 - (a2 - b2)*CompX*CompY/4;
		double C = a2 * X + CompX * ( X + Y - 1 - a2 + b2 + c2 - X*CompY / 2.0 ) /2.0;
		return F*F * A + G*G * B - F*G * C;
	}
}

double ProductionRates::I_MesonThree_s(double s)	//fixing one variable
{
	double fPlus = 1 - GetLambda1() * (Ex + Ey - 1 + c2) / c2;
	double fMinus = (GetLambda0()-GetLambda1()) * (1-c2) / c2;
	double F = 2*fp;
	double G = fp - fm;

	double Lambda = sqrt(Kine::Kallen(1, s, I_var.at(4)));

	double dY2 = (SupY*SupY - InfY*Inf* ) / 2.0;
	double dY1 = SupY - InfY;

	double Comb = 1+I_var.at(5) - s - Lambda;

	return F*F


	if (GetHelicity() < 0)
	{
		double A = (1+a2-s) * dY2 - Comb * dY2 / 2.0 + (s + b2 - c2)*dY1 - dY2;
		double B = b2 * (dY2 - (s +b2 -c2)*dY1) + 2*a2*b2 + (a2-b2) * Comb * dY2 / 4.0;
		double C = b2*(1+a2-s)*dY1 - Comp / 2.0 * ( (1 - a2 + s) * dY2 / 2.0 - (s - b2 - c2) * dY1 );
		return I_var.at(0) + (I_var.at(1) - I_var.at(0)) * 
			(F*F*A + G*G*B - F*G*C)
	else if (GetHelicity() > 0)
	double A_p = Comb * dY2 / 2.0;
	double B_p = a2 * (dY2 - (s +b2 -c2)*dY1) + 2*a2*b2 - (a2-b2) * Comb * dY2 / 4.0;
	double C_p = a2*dY2 + Comp / 2.0 * ( (1 - a2 + s) * dY2 / 2.0 - (s - b2 - c2) * dY1 );
		return I_var.at(0) + (I_var.at(1) - I_var.at(0)) * 
			I_var.at(2)*I_var.at(3) * Sum * Lambda0 * Lambda1 / s;

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


void ProductionRates::InitConst()
{
	M_Sterile_prev = -1.0;
	M_Parent_prev = -1.0;
	U_e_prev = -1.0;
	U_m_prev = -1.0;
	U_t_prev = -1.0;
	fMax = -1.0;

	switch(mapParent[GetParent()])
	{
		case _Muon:	//Muon decays into nu_mu (c), nu_e (a,x), e(b,y)
			M_Parent = M_Muon;

			if (IsElectron && !IsMuon) //nu_e is the sterile
			{
				fA = M_Sterile/M_Parent;
				fC = M_Neutrino/M_Parent;
			}
			if (!IsElectron && IsMuon) //nu_mu is the sterile
			{
				fA = M_Neutrino/M_Parent;
				fC = M_Sterile/M_Parent;
			}

			fB = M_Electron/M_Parent;

			break;

		case _TauE:	//Tau decays into nu_tau (c), nu_e (a,x), e(b,y)
			M_Parent = M_Tau;

			if (IsElectron && !IsTau) //nu_e is the sterile
			{
				fA = M_Sterile/M_Parent;
				fC = M_Neutrino/M_Parent;
			}
			if (!IsElectron && IsTau) //nu_tau is the sterile
			{
				fA = M_Neutrino/M_Parent;
				fC = M_Sterile/M_Parent;
			}

			fB = M_Electron/M_Parent;

			break;

		case _TauM:	//Tau decays into nu_tau (c), nu_mu (a,x), mu(b,y)
			M_Parent = M_Tau;

			if (IsMuon && !IsTau) //nu_e is the sterile
			{
				fA = M_Sterile/M_Parent;
				fC = M_Neutrino/M_Parent;
			}
			if (!IsMuon && IsTau) //nu_tau is the sterile
			{
				fA = M_Neutrino/M_Parent;
				fC = M_Sterile/M_Parent;
			}

			fB = M_Muon/M_Parent;

			break;

		case _Kaon:	//Kaon decays into pi0 (c), lepton (a,x), neutrino (b,y)
			M_Parent = M_Kaon;

			if (IsElectron)
				fA = M_Electron/M_Kaon;
			else if (IsMuon)
				fA = M_Muon/M_Kaon;
			else fA = 0.0;

			fB = M_Sterile/M_Kaon;
			fC = M_Pion0/M_Kaon;

			break;

		case _Kaon0:	//Kaon0 decays into pi (c), lepton (a,x), neutrino (b,y)
			M_Parent = M_Kaon0;

			if (IsElectron)
				fA = M_Electron/M_Kaon;
			else if (IsMuon)
				fA = M_Muon/M_Kaon;
			else fA = 0.0;

			fB = M_Sterile/M_Kaon0;
			fC = M_Pion/M_Kaon0;

			break;

		case _nEE:
			M_Parent = M_Sterile;

			fA = M_Electron/M_Sterile;
			fB = M_Electron/M_Sterile;
			fC = M_Neutrino/M_Sterile;

			break;

		case _nMUMU:
			M_Parent = M_Sterile;

			fA = M_Muon/M_Sterile;
			fB = M_Muon/M_Sterile;
			fC = M_Neutrino/M_Sterile;

			break;

		case _nEMU:	//e+ mu-
			M_Parent = M_Sterile;

			fA = M_Electron/M_Sterile;
			fB = M_Muon/M_Sterile;
			fC = M_Neutrino/M_Sterile;

			break;

		case _nMUE:	//mu+ e- 
			M_Parent = M_Sterile;
			fA = M_Muon/M_Sterile;
			fB = M_Electron/M_Sterile;
			fC = M_Neutrino/M_Sterile;

			break;

		default:
			M_Parent = 0.0;
			fA = 0.0;
			fB = 0.0;
			fC = 0.0;
			break;
	}
}

double ProductionRates::ddGamma()	//double differential decay width (dG/dExdEy)
{
	if (IsEnergyConserved() && InLimX() && InLimY())
		return M2() / (64.0 * Const::fPi3 * GetParentMass());
	else return 0.0;
}

double ProductionRates::dGamma()	//differential decay width (dG/dEx)
{
	if (IsEnergyConserved() && InLimX())
		return M2IntY() / (64.0 * Const::fPi3 * GetParentMass());
		//return M2IntY() * GetParentMass() / (256 * Const::fPi3);
	else return 0.0;
}

double ProductionRates::Gamma()	//fully integrated decay width (G)
{
	if (IsEnergyConserved())
		return M2IntXY() / (64.0 * Const::fPi3 * GetParentMass());
	else return 0.0;
}

double ProductionRates::ddPhaseSpace()	//double differential phase space (dPS/dExdEy)
{
	if (IsEnergyConserved() && InLimX() && InLimY())
		return 1.0;
	else return 0.0;
}

double ProductionRates::dPhaseSpace()	//differential phase space (dPS/dEx)
{
	if (IsEnergyConserved() && InLimX())
	{
		double ymin, ymax;
		return yLim(ymin, ymax) * GetParentMass() / 2.0;
	}
	else return 0.0;
}

double ProductionRates::PhaseSpace()	//fully integrated phase space (PS)
{
	if (IsEnergyConserved())
	{
		double xmin, xmax;
		double dx = xLim(xmin, xmax);

		double (ProductionRates::*pPS)() = &ProductionRates::dPhaseSpace;
		return Integrate(pPS, xmin, xmax) * GetParentMass() / 2.0;
	}
	else return 0.0;
}

double ProductionRates::M2()		//Unpolarised amplitude
{
	if (IsEnergyConserved() && InLimX() && InLimY())
	{
		double M2;
		switch(mapParent[GetParent()])
		{
			case _Muon:
			case _TauE:
			case _TauM:
				M2 = M2Lept();
				break;
			case _Kaon:
				M2 = M2Kaon();
				break;
			case _Kaon0:
				M2 = M2Kaon0();
				break;
			case _nEE:
				M2 = M2nEE(); 
				break;
			case _nMUMU:
				M2 = M2nMUMU(); 
				break;
			case _nEMU:
				M2 = M2nEMU(); 
				break;
			case _nMUE:
				M2 = M2nMUE(); 
				break;
			default:
				M2 = 0.0;
				break;
		}

		return M2;
	}
	else return 0.0;
}

double ProductionRates::M2IntY()	//Unpolarised amplitude, integrated over Ey
{
	if (IsEnergyConserved() && InLimX())
	{
		double M2;
		switch(mapParent[GetParent()])
		{
			case _Muon:
			case _TauE:
			case _TauM:
				M2 = M2LeptIntY();
				break;
			case _Kaon:
				M2 = M2KaonIntY();
				break;
			case _Kaon0:
				M2 = M2Kaon0IntY();
				break;
			case _nEE:
				M2 = M2nEEIntY();
				break;
			default:
				M2 = 0.0;
				break;
		}

		return M2 * GetParentMass() / 2.0;
	}
	else return 0.0;
}

double ProductionRates::M2IntXY()	//Unpolarised amplitude, integrated over Ex and Ey
{
	if (IsEnergyConserved())
	{
		double xmin, xmax;
		double dx = xLim(xmin, xmax);

		double (ProductionRates::*pM2)() = &ProductionRates::M2IntY;
		return Integrate(pM2, xmin, xmax) * GetParentMass() / 2.0;
	}
	else return 0.0;
}

//Unpolarised amplitudes here after

double ProductionRates::M2_Z() //NC N to n l l, Z propagator
{
	double gV = -0.5 + 2 * Const::fSin2W;
	double gA = -0.5;
	return 16 * Const::fGF2 * pow(M_Sterile, 4) * 
		(   pow((gV+gA), 2) * x() * (1 + a(2) - b(2) - c(2) - x())
		  + pow((gV-gA), 2) * y() * (1 + b(2) - a(2) - c(2) - y())
		  + (gV*gV - gA*gA) * a()*b() * (2 - x() - y()) );
}

double ProductionRates::M2_ZIntY() //NC N to n l l, Z propagator
{
	double ymin, ymax;
	double dy = yLim(ymin, ymax);

	double gV = -0.5 + 2 * Const::fSin2W;
	double gA = -0.5;
	return 16 * Const::fGF2 * pow(M_Sterile, 4) * 
		( dy * ( pow((gV+gA), 2) * x() * (1 + a(2) - b(2) - c(2) - x()) +
		         (gV*gV - gA*gA) * a()*b() * (2 - x()) ) +
		  (ymax*ymax - ymin*ymin) / 2.0 * ( pow((gV-gA), 2) * (1 + b(2) - a(2) - c(2)) -
						    (gV*gV - gA*gA) * a()*b() ) -
		  (ymax*ymax*ymax - ymin*ymin*ymin) / 3.0 * pow((gV-gA), 2) );
}

double ProductionRates::M2_WZ() //Interference term between Z and W propagator
{
	double gV = -0.5 + 2 * Const::fSin2W;
	double gA = -0.5;
	return  8.0 * Const::fGF2 * pow(M_Sterile, 4) * 
	       ( - (gV + gA) * x() * (1 + a(2) - b(2) - c(2) - x())
		 + (gV - gA) * a()*b() * (2 - x() - y()) );
}

double ProductionRates::M2_WZIntY() //Interference term between Z and W propagator
{
	double ymin, ymax;
	double dy = yLim(ymin, ymax);

	double gV = -0.5 + 2 * Const::fSin2W;
	double gA = -0.5;

	return  8.0 * Const::fGF2 * pow(M_Sterile, 4) * 
	       ( dy * (- (gV + gA) * x() * (1 + a(2) - b(2) - c(2) - x()) +
		       (gV - gA) * a()*b() * (2 - x()) ) -
		 (ymax*ymax - ymin*ymin) / 2.0 * (gV - gA) * a()*b() );
}

double ProductionRates::M2Lept()	//Pure leptonic decays (like muon or tau) 	//W propagator
{
	return 16 * Const::fGF2 * GetUu()*GetUu() *
	       	pow(GetParentMass(), 4) * x() * (1 + a(2) - b(2) - c(2) - x());
}

double ProductionRates::M2LeptIntY()	//Leptonic decay, integrated analytically over y
{
	double ymin, ymax;
	return yLim(ymin, ymax) * M2Lept();
}

double ProductionRates::M2Kaon()	//Kaon decay
{
	return Const::fGF2 * GetUu()*GetUu() * pow(M_Kaon,4) * pow(GetDecayConst(),2) * 
		( 4 * fPlus()*fPlus() * ( 1 + a(2) + b(2) - c(2) - x() - y() + x()*y() ) -
	      	  pow(fPlus() - fMinus(), 2) * ( pow(a(2) - b(2),2) + (a(2)+b(2)) * (1 - c(2) - x() - y()) ) );
}

double ProductionRates::M2KaonIntY(double Y)	//Kaon decay primitive
{
	double X = 1 - c(2) - x();
	double AB = a(2)+b(2);
	double A2B = pow(a(2)-b(2),2);
	double L1 = GetLambda1();
	double L0 = GetLambda0();
	double fP = 1 - L1 * (X - Y) / c(2);
	double fM = fMinus();

	double Ret1 = 4 * fP*fP * Y * (X + AB - Y/2 + x()*Y/2) -
		      8 * L1 / c(2) * (1 - L1*X/c(2)) * Y*Y/2 * (X + AB - Y/3 + x()*Y/3) - 
		      8 * L1*L1 / c(4) * Y*Y*Y/3 * (X + AB - Y/4 + x()*Y/4);
	double Ret2 = (fP-fM)*(fP-fM) * Y * (A2B + AB * (X-Y/2)) -
		      2 * L1 / c(2) * (1 - L1*X/c(2) - fM) * Y*Y/2 * (A2B + AB * (X-Y/3)) -
		      2 * L1*L1 / c(4) * Y*Y*Y/3 * (A2B + AB * (X-Y/4));

	return Const::fGF2 * GetUu()*GetUu() * pow(M_Kaon,4) * pow(GetDecayConst(),2) * (Ret1 - Ret2);
}

double ProductionRates::M2KaonIntY()	//Kaon decay, integrated analytically over y
{
	double ymin, ymax;
	double dy = yLim(ymin, ymax);

	double Ret = M2KaonIntY(ymax) - M2KaonIntY(ymin);
	if (Ret > 0)
		return Ret;
	else return 0.0;
}

double ProductionRates::M2Kaon0()	//Kaon 0 decay
{
	return 2*M2Kaon();
}

double ProductionRates::M2Kaon0IntY(double Y)	//Kaon0 decay primitive
{
	return 2*M2KaonIntY(Y);
}

double ProductionRates::M2Kaon0IntY() //Kaon0 decay, integrated analytically over y
{
	double ymin, ymax;
	double dy = yLim(ymin, ymax);

	double Ret = M2Kaon0IntY(ymax) - M2Kaon0IntY(ymin);
	if (Ret > 0)
		return Ret;
	else return 0.0;
}

//Very boring stuff
/*M2 of visible processes with non constant M2 
 * x is for the antilepton, of mass a
 * y is for the lepton, of mass b
 * n is integrated out, with mass c
 * this convention is taken for all M2
 */

double ProductionRates::M2nEE()
{
	return GetUe()*GetUe() * (M2Lept() + M2_WZ() + M2_Z()) + (GetUm()*GetUm() + GetUt()*GetUt()) * M2_Z();
}

double ProductionRates::M2nEEIntY()
{
	return GetUe()*GetUe() * (M2LeptIntY() + M2_WZIntY() + M2_ZIntY()) + (GetUm()*GetUm() + GetUt()*GetUt()) * M2_ZIntY();
}

double ProductionRates::M2nMUMU()
{
	return GetUm()*GetUm() * (M2Lept() + M2_WZ() + M2_Z()) + (GetUe()*GetUe() + GetUt()*GetUt()) * M2_Z();
}

double ProductionRates::M2nEMU()
{
	return GetUm()*GetUm() * M2Lept();
}

double ProductionRates::M2nMUE()
{
	return GetUe()*GetUm() * M2Lept();
}

//Maximum values of ddG for MC purposes
double ProductionRates::MaxGamma()
{
	if (IsChanged() && fMax < 0)
	{
		double temx = x();
		double temy = y();

		fMax = -1.0;

		//Phasespace coordinates are already checked by ddGamma
		//for (double ix = 0.0; ix <= 2.0; ix += 2.0/Kine::Loop)
		for (double ix = 0.0; ix <= 2.0; ix += 2.0/100)
		{
			SetX(ix);
			//for (double iy = 0.0; iy <= 2.0; iy += 2.0/Kine::Loop)
			for (double iy = 0.0; iy <= 2.0; iy += 2.0/100)
			{
				SetY(iy);
				double Gam = ddGamma();
				if (fMax < Gam)
				       fMax = Gam;	
			}
		}

		SetX(temx);
		SetY(temy);
	}

	return fMax;
}

//boundaries of phase space
double ProductionRates::yLim(double &Min, double &Max)	//y integration limits
{
	double X = 1 + a(2) - x();
	double A = (2 - x())*(X + b(2) - c(2));
	double P = x(2) - 4*a(2);
	double L = Kine::Lambda(X, b(2), c(2));

	Min = (A - sqrt(P) * sqrt(L)) / (2*X);
	Max = (A + sqrt(P) * sqrt(L)) / (2*X);

	return Max-Min;
}

double ProductionRates::xLim(double &Min, double &Max)	//x integration limits
{
	Min = 2*a();
	Max = 1 + a(2) - pow((b() + c()),2);

	return Max-Min;
}

double ProductionRates::Integrate(double (ProductionRates::*FF)(), double A, double B)
{
	if (A > B)
	{
		double tmp = B;
		B = A;
		A = tmp;
	}
	double a, b;
	double h = (B-A)/100.0;	//Step
	double Integral = 0;	//Boole's method for integration
	for (a = A; b <= B; a = b)
	{
		b = a + h;
		SetX(a);
		Integral += 7*(this->*FF)();

		SetX((3*a+b)/4.0);
		Integral += 32*(this->*FF)();

		SetX((a+b)/2.0);
		Integral += 12*(this->*FF)();

		SetX((a+3*b)/4.0);
		Integral += 32*(this->*FF)();

		SetX(b);
		Integral += 7*(this->*FF)();
	}	

	return Integral * h/90.0;
}

bool ProductionRates::InLimX()
{
	double xmin, xmax;
	double dx = xLim(xmin, xmax);
	return (xmin <= x() && x() <= xmax);
}

bool ProductionRates::InLimY()
{
	double ymin, ymax;
	double dy = yLim(ymin, ymax);
	return (ymin <= y() && y() <= ymax);
}

//The most important question, after all
bool ProductionRates::IsEnergyConserved()
{
	if (a()+b()+c() <= 1.0)
		return true;
	else return false;
}

//decay constants
double ProductionRates::fPlus()
{
	return 1 - GetLambda1() * (1 + c(2) - x() - y()) / c(2);
}

double ProductionRates::fMinus()
{
	return (GetLambda0()-GetLambda1()) * (1-c(2)) / c(2);
}

//variables
double ProductionRates::a(double p)	//mass ratio
{
	return pow(fA,p);
}

double ProductionRates::b(double p)	//mass ratio
{
	return pow(fB,p);
}

double ProductionRates::c(double p)	//mass ratio
{
	return pow(fC,p);
}

double ProductionRates::x(double p)	//energy over mass
{
	double Ret = 2*GetEnergyX()/GetParentMass();
	return pow(Ret,p);
}

double ProductionRates::y(double p)	//energy over mass
{
	double Ret = 2*GetEnergyY()/GetParentMass();
	return pow(Ret,p);
}

//Channel mode
void ProductionRates::ElectronChannel()
{
	IsElectron = true;
	IsMuon = false;
	IsTau = false;
	InitConst();
}

void ProductionRates::MuonChannel()
{
	IsElectron = false;
	IsMuon = true;
	IsTau = false;
	InitConst();
}

void ProductionRates::TauChannel()
{
	IsElectron = false;
	IsMuon = false;
	IsTau = true;
	InitConst();
}

//Getter
std::string ProductionRates::GetParent()
{
	return sParent;
}

double ProductionRates::GetEnergyX()
{
	return fEX;
}

double ProductionRates::GetEnergyY()
{
	return fEY;
}

double ProductionRates::GetParentMass()
{
	return M_Parent;
}

double ProductionRates::GetSterileMass()
{
	return M_Sterile;
}

double ProductionRates::GetUe()
{
	return U_e;
}

double ProductionRates::GetUm()
{
	return U_m;
}

double ProductionRates::GetUt()
{
	return U_t;
}

double ProductionRates::GetUu()
{
	if (IsElectron)
		return GetUe();
	else if (IsMuon)
		return GetUm();
	else if (IsTau)
		return GetUt();
	else return 1.0;
}

double ProductionRates::GetDecayConst()
{
	return fKaon;
}

double ProductionRates::GetLambda1()
{
	return fLambda1;
}

double ProductionRates::GetLambda0()
{
	return fLambda0;
}

//Setter
void ProductionRates::SetParent(std::string Name)
{
	sParent.assign(Name);
	InitConst();
}

void ProductionRates::SetEnergyX(double X)
{
	fEX = X;
}

void ProductionRates::SetEnergyY(double X)
{
	fEY = X;
}

void ProductionRates::SetX(double X)
{
	fEX = X*GetParentMass()/2;
}

void ProductionRates::SetY(double X)
{
	fEY = X*GetParentMass()/2;
}

void ProductionRates::SetSterileMass(double X)
{
	M_Sterile = X;
	InitConst();
}

void ProductionRates::SetUe(double X)
{
	U_e = X;
}

void ProductionRates::SetUm(double X)
{
	U_m = X;
}

void ProductionRates::SetUt(double X)
{
	U_t = X;
}

bool ProductionRates::IsChanged()
{
	bool Ret = ( M_Sterile != M_Sterile_prev || 
		     M_Parent != M_Parent_prev || 
  		     U_e != U_e_prev ||
  		     U_m != U_m_prev ||
  		     U_t != U_t_prev );

	M_Sterile_prev = M_Sterile;
	M_Parent_prev = M_Parent; 
	U_e_prev = U_e;
	U_m_prev = U_m;
	U_t_prev = U_t;

	return Ret;
}
