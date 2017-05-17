#include "3Body.h"

ThreeBody::ThreeBody(std::string Parent, double MSterile, double Ue, double Um, double Ut)	: //ThreeBody rates calculator
	M_Neutrino(0.0),
	M_Photon(0.0),
	//M_Electron(Const::fMElectron),
	M_Electron(0.0),
	M_Muon(Const::fMMuon),
	M_Pion(Const::fMPion),
	M_Pion0(Const::fMPion0),
	M_Kaon(Const::fMKaon),
	M_Kaon0(Const::fMKaon0),
	fKaon(Const::fVusFKaon),
	fLambda1(Const::fLambda1m),
	fLambda0(Const::fLambda0m)
{
	SetSterileMass(MSterile);
	SetUe(Ue);
	SetUm(Um);
	SetUt(Ut);

	SetParent(Parent);

	IsElectron = false;
	IsMuon = false;

	InitMap();
	InitConst();
}

//Initialisation of map
void ThreeBody::InitMap()
{
	mapParent["Muon"] = _Muon;
	mapParent["Kaon"] = _Kaon;
	mapParent["Kaon0"] = _Kaon0;
}

void ThreeBody::InitConst()
{
	switch(mapParent[GetParent()])
	{
		case _Muon:	//Will need to change this including heavy neutrino
			M_Parent = M_Muon;
			fA = M_Neutrino/M_Muon;
			fB = M_Electron/M_Muon;
			fC = M_Neutrino/M_Muon;
			break;
		case _Kaon:	//Will need to change this including heavy neutrino
			M_Parent = M_Kaon;
			if (IsElectron)
				fA = M_Electron/M_Kaon;
			else if (IsMuon)
				fA = M_Muon/M_Kaon;
			fB = M_Neutrino/M_Kaon;
			fC = M_Pion0/M_Kaon;
			break;
		case _Kaon0:	//Will need to change this including heavy neutrino
			M_Parent = M_Kaon0;
			if (IsElectron)
				fA = M_Electron/M_Kaon;
			else if (IsMuon)
				fA = M_Muon/M_Kaon;
			fB = M_Neutrino/M_Kaon0;
			fC = M_Pion0/M_Kaon0;
			break;
		default:
			M_Parent = 0.0;
			fA = 0.0;
			fB = 0.0;
			fC = 0.0;
			std::cout << "Unknown Parent particle!" << std::endl;
			break;
	}

	if (!IsEnergyConserved())
	{
		M_Parent = 0.0;
		fA = 0.0;
		fB = 0.0;
		fC = 0.0;
		std::cout << "Warning! Respect the rules, check the masses." << std::endl;
	}
}

double ThreeBody::ddGamma()	//double differential decay width (dG/dxdy)
{
	if (InLimX() && InLimY())
		return M2() * GetParentMass() / (256 * genie::constants::kPi3);
	else return 0.0;
}

double ThreeBody::dGamma()	//differential decay width (dG/dx)
{
	if (InLimX())
		return M2IntY() * GetParentMass() / (256 * genie::constants::kPi3);
	else return 0.0;
}

double ThreeBody::Gamma()	//fully integrated decay width (G)
{
	double xmin, xmax;
	double dx = xLim(xmin, xmax);
	std::cout << "xmin " << xmin << " xmax " << xmax << std::endl;

	double (ThreeBody::*pGamma)() = &ThreeBody::dGamma;
	return Integrate(pGamma, xmin, xmax); 
}

double ThreeBody::M2()		//Unpolarised amplitude
{
	double M2 = 0.0;
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
		default:
			break;
	}

	if (InLimX() && InLimY())
		return M2;
	else return 0.0;
}

double ThreeBody::M2IntY()	//Unpolarised amplitude, integrated over y
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
		default:
			M2 = 0.0;
			break;
	}

	if (InLimX())
		return M2;
	else return 0.0;
}

double ThreeBody::M2IntXY()	//Unpolarised amplitude, integrated over y and x
{
	double xmin, xmax;
	double dx = xLim(xmin, xmax);

	double (ThreeBody::*pM2)() = &ThreeBody::M2IntY;
	return Integrate(pM2, xmin, xmax);
}

//Unpolarised amplitudes here after
double ThreeBody::M2Muon()	//Muon decay
{
	return 16 * Const::fGF2 * GetUu()*GetUu() *
	       	pow(M_Muon, 4) * x() * (1 + a2() - b2() - c2() - x());
}

double ThreeBody::M2MuonIntY()	//Muon decay, integrated analytically over y
{
	double ymin, ymax;
	double dy = yLim(ymin, ymax);
	return dy * M2Muon();
}

double ThreeBody::M2Kaon()	//Kaon decay
{
	double Zeta = fZeta()-1;

	return 2 * Const::fGF2 * GetUu()*GetUu() *
		pow(M_Kaon,4) * pow(GetDecayConst(),2) * fPlus()*fPlus() *
		(1 - c2() +x()*y() - x() - y() + a2()*(1 + y()*Zeta) + b2()*(1 + x()*Zeta) -
		 4*Zeta*Zeta*(a2() + b2()) * (1 + a2() - b2() - c2() - x() - y()));
}

double ThreeBody::M2KaonIntY(double Y)	//Kaon decay primitive
{
	double Zeta = fZeta()-1;

	return 2 * Const::fGF2 * GetUu()*GetUu() *
		pow(M_Kaon,4) * pow(GetDecayConst(),2) * fPlus()*fPlus() * Y *
		(1 - c2() + x()*Y/2.0 - x() - Y/2.0 + a2()*(1 + Y*Zeta/2.0) + b2()*(1 + x()*Zeta) -
		 4*Zeta*Zeta*(a2() + b2()) * (1 + a2() - b2() - c2() - x() - Y/2.0));
}

double ThreeBody::M2KaonIntY()	//Kaon decay, itnegrated analytically over y
{
	double ymin, ymax;
	double dy = yLim(ymin, ymax);

	return M2KaonIntY(ymax) - M2KaonIntY(ymin);
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

	return M2Kaon0IntY(ymax) - M2Kaon0IntY(ymin);
}

double ThreeBody::yLim(double &Min, double &Max)	//y integration limits
{
	double X = 1 + a2() - x();
	double A = (2 - x())*(X + b2() - c2());
	double P = x2() - 4*a2();
	double L = Kine::Lambda(X, b2(), c2());

	Min = (A - sqrt(P) * sqrt(L)) / (2*X);
	Max = (A + sqrt(P) * sqrt(L)) / (2*X);

	return Max-Min;
}

double ThreeBody::xLim(double &Min, double &Max)	//x integration limits
{
	Min = 2*a();
	Max = 1 + a2() - pow((b() + c()),2);

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
	double h = (B-A)/Kine::Sample;	//Step
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
	if (a()+b()+c() <= 1)
		return true;
	else false;
}

//decay constants
double ThreeBody::fPlus()
{
	return 1+(c2() + x() + y() - 1)*GetLambda1()/c2();
}

double ThreeBody::fMinus()
{
	return (GetLambda0()-GetLambda1())*(1-c2())/(c2());
}

double ThreeBody::fZeta()
{
	return fMinus()/fPlus();
}


//variables
double ThreeBody::a()	//mass ratio
{
	return fA;
}

double ThreeBody::b()	//mass ratio
{
	return fB;
}

double ThreeBody::c()	//mass ratio
{
	return fC;
}

double ThreeBody::a2()
{
	return a()*a();
}

double ThreeBody::b2()
{
	return b()*b();
}

double ThreeBody::c2()
{
	return c()*c();
}

double ThreeBody::x()	//energy over mass
{
	return 2*GetEnergyX()/GetParentMass();
}

double ThreeBody::y()	//energy over mass
{
	return 2*GetEnergyY()/GetParentMass();
}

double ThreeBody::x2()
{
	return x()*x();
}

double ThreeBody::y2()
{
	return y()*y();
}

//Channel mode
void ThreeBody::ElectronChannel()
{
	IsElectron = true;
	IsMuon = false;
}

void ThreeBody::MuonChannel()
{
	IsElectron = false;
	IsMuon = true;
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
		return GetUu();
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
