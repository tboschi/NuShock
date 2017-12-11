#include "3Body.h"

ThreeBody::ThreeBody(std::string Parent, double MSterile, double Ue, double Um, double Ut)	: //ThreeBody rates calculator
	M_Neutrino(0.0),
	M_Photon(0.0),
	M_Electron(Const::fMElectron),
	//M_Electron(0.0),
	M_Muon(Const::fMMuon),
	M_Pion(Const::fMPion),
	M_Pion0(Const::fMPion0),
	M_Kaon(Const::fMKaon),
	M_Kaon0(Const::fMKaon0),
	fKaon(Const::fVusFKaon),
	fLambda1(Const::fLambdaPlus),
	fLambda0(Const::fLambdaZero)
{
	InitMap();
	SetParent(Parent);

	SetSterileMass(MSterile);
	SetUe(Ue);
	SetUm(Um);
	SetUt(Ut);

	IsElectron = false;
	IsMuon = false;

	InitConst();
}

//Initialisation of map
void ThreeBody::InitMap()
{
	mapParent[""] = _undefined;
	mapParent["Muon"] = _Muon;
	mapParent["Kaon"] = _Kaon;
	mapParent["Kaon0"] = _Kaon0;
	mapParent["nEE"] = _nEE;
	mapParent["nMUMU"] = _nMUMU;
	mapParent["nEMU"] = _nEMU;
	mapParent["nMUE"] = _nMUE;
}

void ThreeBody::InitConst()
{
	switch(mapParent[GetParent()])
	{
		case _Muon:	//Muon decays into nu_mu (c), nu_e (a,x), e(b,y)
			M_Parent = M_Muon;

			if (IsElectron && !IsMuon) //nu_e is the sterile
			{
				fA = M_Sterile/M_Muon;
				fC = M_Neutrino/M_Muon;
			}
			if (!IsElectron && IsMuon) //nu_mu is the sterile
			{
				fA = M_Neutrino/M_Muon;
				fC = M_Sterile/M_Muon;
			}

			fB = M_Electron/M_Muon;

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

double ThreeBody::ddGamma()	//double differential decay width (dG/dExdEy)
{
	if (IsEnergyConserved() && InLimX() && InLimY())
		return M2() / (64.0 * Const::fPi3 * GetParentMass());
	else return 0.0;
}

double ThreeBody::dGamma()	//differential decay width (dG/dEx)
{
	if (IsEnergyConserved() && InLimX())
		return M2IntY() / (64.0 * Const::fPi3 * GetParentMass());
		//return M2IntY() * GetParentMass() / (256 * Const::fPi3);
	else return 0.0;
}

double ThreeBody::Gamma()	//fully integrated decay width (G)
{
	if (IsEnergyConserved())
		return M2IntXY() / (64.0 * Const::fPi3 * GetParentMass());
	else return 0.0;
}

double ThreeBody::ddPhaseSpace()	//double differential phase space (dPS/dExdEy)
{
	if (IsEnergyConserved() && InLimX() && InLimY())
		return 1.0;
	else return 0.0;
}

double ThreeBody::dPhaseSpace()	//differential phase space (dPS/dEx)
{
	if (IsEnergyConserved() && InLimX())
	{
		double ymin, ymax;
		return yLim(ymin, ymax) * GetParentMass() / 2.0;
	}
	else return 0.0;
}

double ThreeBody::PhaseSpace()	//fully integrated phase space (PS)
{
	if (IsEnergyConserved())
	{
		double xmin, xmax;
		double dx = xLim(xmin, xmax);

		double (ThreeBody::*pPS)() = &ThreeBody::dPhaseSpace;
		return Integrate(pPS, xmin, xmax) * GetParentMass() / 2.0;
	}
	else return 0.0;
}

double ThreeBody::M2()		//Unpolarised amplitude
{
	if (IsEnergyConserved() && InLimX() && InLimY())
	{
		double M2;
		switch(mapParent[GetParent()])
		{
			case _Muon:
				M2 = M2Muon();
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

double ThreeBody::M2IntY()	//Unpolarised amplitude, integrated over Ey
{
	if (IsEnergyConserved() && InLimX())
	{
		double M2;
		switch(mapParent[GetParent()])
		{
			case _Muon:
				M2 = M2MuonIntY();
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

double ThreeBody::M2IntXY()	//Unpolarised amplitude, integrated over Ex and Ey
{
	if (IsEnergyConserved())
	{
		double xmin, xmax;
		double dx = xLim(xmin, xmax);

		double (ThreeBody::*pM2)() = &ThreeBody::M2IntY;
		return Integrate(pM2, xmin, xmax) * GetParentMass() / 2.0;
	}
	else return 0.0;
}

//Unpolarised amplitudes here after

double ThreeBody::M2_Z() //NC N to n l l, Z propagator
{
	double gV = -0.5 + 2 * Const::fSin2W;
	double gA = -0.5;
	return 16 * Const::fGF2 * pow(M_Sterile, 4) * 
		(   pow((gV+gA), 2) * x() * (1 + a(2) - b(2) - c(2) - x())
		  + pow((gV-gA), 2) * y() * (1 + b(2) - a(2) - c(2) - y())
		  + (gV*gV - gA*gA) * a()*b() * (2 - x() - y()) );
}

double ThreeBody::M2_ZIntY() //NC N to n l l, Z propagator
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

double ThreeBody::M2_WZ() //Interference term between Z and W propagator
{
	double gV = -0.5 + 2 * Const::fSin2W;
	double gA = -0.5;
	return  8.0 * Const::fGF2 * pow(M_Sterile, 4) * 
	       ( - (gV + gA) * x() * (1 + a(2) - b(2) - c(2) - x())
		 + (gV - gA) * a()*b() * (2 - x() - y()) );
}

double ThreeBody::M2_WZIntY() //Interference term between Z and W propagator
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

double ThreeBody::M2Muon()	//Muon decay 	//W propagator
{
	return 16 * Const::fGF2 * GetUu()*GetUu() *
	       	pow(M_Muon, 4) * x() * (1 + a(2) - b(2) - c(2) - x());
}

double ThreeBody::M2MuonIntY()	//Muon decay, integrated analytically over y
{
	double ymin, ymax;
	return yLim(ymin, ymax) * M2Muon();
}

double ThreeBody::M2Kaon()	//Kaon decay
{
	return Const::fGF2 * GetUu()*GetUu() * pow(M_Kaon,4) * pow(GetDecayConst(),2) * 
		( 4 * fPlus()*fPlus() * ( 1 + a(2) + b(2) - c(2) - x() - y() + x()*y() ) -
	      	  pow(fPlus() - fMinus(), 2) * ( pow(a(2) - b(2),2) + (a(2)+b(2)) * (1 - c(2) - x() - y()) ) );
}

double ThreeBody::M2KaonIntY(double Y)	//Kaon decay primitive
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

double ThreeBody::M2KaonIntY()	//Kaon decay, itnegrated analytically over y
{
	double ymin, ymax;
	double dy = yLim(ymin, ymax);

	double Ret = M2KaonIntY(ymax) - M2KaonIntY(ymin);
	if (Ret > 0)
		return Ret;
	else return 0.0;
}

double ThreeBody::M2Kaon0()	//Kaon 0 decay
{
	return 2*M2Kaon();
}

double ThreeBody::M2Kaon0IntY(double Y)	//Kaon0 decay primitive
{
	return 2*M2KaonIntY(Y);
}

double ThreeBody::M2Kaon0IntY() //Kaon0 decay, integrated analytically over y
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

double ThreeBody::M2nEE()
{
	return GetUe()*GetUe() * (M2Muon() + M2_WZ() + M2_Z()) + (GetUm()*GetUm() + GetUt()*GetUt()) * M2_Z();
}

double ThreeBody::M2nEEIntY()
{
	return GetUe()*GetUe() * (M2MuonIntY() + M2_WZIntY() + M2_ZIntY()) + (GetUm()*GetUm() + GetUt()*GetUt()) * M2_ZIntY();
}

double ThreeBody::M2nMUMU()
{
	return GetUm()*GetUm() * (M2Muon() + M2_WZ() + M2_Z()) + (GetUe()*GetUe() + GetUt()*GetUt()) * M2_Z();
}

double ThreeBody::M2nEMU()
{
	return GetUm()*GetUm() * M2Muon();
}

double ThreeBody::M2nMUE()
{
	return GetUe()*GetUm() * M2Muon();
}

//Maximum values of ddG for MC purposes
double ThreeBody::MaxGamma()
{
	double temx = x();
	double temy = y();

	if (IsChanged())
	{
		fMax = 0.0;

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
	}

	SetX(temx);
	SetY(temy);

	return fMax;
}

//boundaries of phase space
double ThreeBody::yLim(double &Min, double &Max)	//y integration limits
{
	double X = 1 + a(2) - x();
	double A = (2 - x())*(X + b(2) - c(2));
	double P = x(2) - 4*a(2);
	double L = Kine::Lambda(X, b(2), c(2));

	Min = (A - sqrt(P) * sqrt(L)) / (2*X);
	Max = (A + sqrt(P) * sqrt(L)) / (2*X);

	return Max-Min;
}

double ThreeBody::xLim(double &Min, double &Max)	//x integration limits
{
	Min = 2*a();
	Max = 1 + a(2) - pow((b() + c()),2);

	return Max-Min;
}

double ThreeBody::Integrate(double (ThreeBody::*FF)(), double A, double B)
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

bool ThreeBody::InLimX()
{
	double xmin, xmax;
	double dx = xLim(xmin, xmax);
	return (xmin <= x() && x() <= xmax);
}

bool ThreeBody::InLimY()
{
	double ymin, ymax;
	double dy = yLim(ymin, ymax);
	return (ymin <= y() && y() <= ymax);
}

//The most important question, after all
bool ThreeBody::IsEnergyConserved()
{
	if (a()+b()+c() <= 1.0)
		return true;
	else return false;
}

//decay constants
double ThreeBody::fPlus()
{
	return 1 - GetLambda1() * (1 + c(2) - x() - y()) / c(2);
}

double ThreeBody::fMinus()
{
	return (GetLambda0()-GetLambda1()) * (1-c(2)) / c(2);
}

//variables
double ThreeBody::a(double p)	//mass ratio
{
	return pow(fA,p);
}

double ThreeBody::b(double p)	//mass ratio
{
	return pow(fB,p);
}

double ThreeBody::c(double p)	//mass ratio
{
	return pow(fC,p);
}

double ThreeBody::x(double p)	//energy over mass
{
	double Ret = 2*GetEnergyX()/GetParentMass();
	return pow(Ret,p);
}

double ThreeBody::y(double p)	//energy over mass
{
	double Ret = 2*GetEnergyY()/GetParentMass();
	return pow(Ret,p);
}

//Channel mode
void ThreeBody::ElectronChannel()
{
	IsElectron = true;
	IsMuon = false;
	InitConst();
}

void ThreeBody::MuonChannel()
{
	IsElectron = false;
	IsMuon = true;
	InitConst();
}

//Getter
std::string ThreeBody::GetParent()
{
	return sParent;
}

double ThreeBody::GetEnergyX()
{
	return fEX;
}

double ThreeBody::GetEnergyY()
{
	return fEY;
}

double ThreeBody::GetParentMass()
{
	return M_Parent;
}

double ThreeBody::GetSterileMass()
{
	return M_Sterile;
}

double ThreeBody::GetUe()
{
	return U_e;
}

double ThreeBody::GetUm()
{
	return U_m;
}

double ThreeBody::GetUt()
{
	return U_t;
}

double ThreeBody::GetUu()
{
	if (IsElectron)
		return GetUe();
	else if (IsMuon)
		return GetUm();
	else return 1.0;
}

double ThreeBody::GetDecayConst()
{
	return fKaon;
}

double ThreeBody::GetLambda1()
{
	return fLambda1;
}

double ThreeBody::GetLambda0()
{
	return fLambda0;
}

//Setter
void ThreeBody::SetParent(std::string Name)
{
	sParent.assign(Name);
	InitConst();
}

void ThreeBody::SetEnergyX(double X)
{
	fEX = X;
}

void ThreeBody::SetEnergyY(double X)
{
	fEY = X;
}

void ThreeBody::SetX(double X)
{
	fEX = X*GetParentMass()/2;
}

void ThreeBody::SetY(double X)
{
	fEY = X*GetParentMass()/2;
}

void ThreeBody::SetSterileMass(double X)
{
	M_Sterile = X;
	InitConst();
}

void ThreeBody::SetUe(double X)
{
	U_e = X;
}

void ThreeBody::SetUm(double X)
{
	U_m = X;
}

void ThreeBody::SetUt(double X)
{
	U_t = X;
}

bool ThreeBody::IsChanged()
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
