#include "Hadron.h"

Hadron::Hadron(std::string ProbeName, std::string TargetName, int OutQuark)
{
	ProbPdf = LHAPDF::mkPDF(ProbeName);
	TargPdf = LHAPDF::mkPDF(TargetName);

	//Process   : A1 + A2 -> q + q_ + X
	//Hard proc : parton1 + parton2 -> quark3 + quark_4
	Prob = new TLorentzVector();	//known from the start
	Targ = new TLorentzVector();	//known from the start

	//massless partons
	Part1 = new TLorentzVector();	//parton 1 (quark or gluon)
	Part2 = new TLorentzVector();	//parton 2 (quark or gluon)

	Quark3 = new TLorentzVector();	//quark, specific flavour
	Quark4 = new TLorentzVector();	//anti quark, same flavour

	PartonPS = new TGenPhaseSpace;

	ValQuark = OutQuark;
	switch(ValQuark)
	{
		case 1:
			ValMass = Const::fMQuarkU;
			break;
		case 2:
			ValMass = Const::fMQuarkD;
			break;
		case 3:
			ValMass = Const::fMQuarkS;
			break;
		case 4:
			//ValMass = Const::fMQuarkC;
			ValMass = 1.27;
			break;
		case 5:
			ValMass = Const::fMQuarkB;
			break;
		case 6:
			ValMass = Const::fMQuarkT;
			break;
		default:
			break;
	}
}

//Phase Space setting
void Hadron::SetProb(double Px, double Py, double Pz, double E)
{
	Prob->SetPxPyPzE(Px, Py, Pz, E);
}

void Hadron::SetTarg(double Px, double Py, double Pz, double E)
{
	Targ->SetPxPyPzE(Px, Py, Pz, E);
}

void Hadron::SetParton(double s)	//Probe on Target, fixed by experiment and constant
{
	double E = sqrt(s)/2.0;

	//in cm
	Part1->SetPxPyPzE(0, 0,  E, E);
	Part2->SetPxPyPzE(0, 0, -E, E);

	double Mass[2] = {M(), M()};

	TLorentzVector S = *Part1 + *Part2;
	if (PartonPS->SetDecay(S, 2, Mass))
	{
		PartonPS->Generate();

		Quark3->SetVect(PartonPS->GetDecay(0)->Vect());
		Quark3->SetE(PartonPS->GetDecay(0)->E());

		Quark4->SetVect(PartonPS->GetDecay(1)->Vect());
		Quark4->SetE(PartonPS->GetDecay(1)->E());
	}
}

void Hadron::SetFromPS()
{
	PartonPS->Generate();

	Quark3->SetVect(PartonPS->GetDecay(0)->Vect());
	Quark3->SetE(PartonPS->GetDecay(0)->E());

	Quark4->SetVect(PartonPS->GetDecay(1)->Vect());
	Quark4->SetE(PartonPS->GetDecay(1)->E());

	TLorentzVector S1 = *Part1 + *Part2;
	TLorentzVector S2 = *Quark3 + *Quark4;
	std::cout << "S " << S1.M2() << "\t" << S2.M2() << std::endl;
	SetCMEnergy(S2.M2());

	TLorentzVector T1 = *Part1 - *Quark3;
	TLorentzVector T2 = *Part2 - *Quark4;
	std::cout << "T " << T1.M2() << "\t" << T2.M2() << std::endl;
	SetQ2(-T1.M2());
}

double Hadron::Pt(int e)
{
	return pow(Quark3->Pt(), e);
}

bool Hadron::PSCut()
{
	if (fabs(CosT_()) < 0.8)
		return true;
	else
		return false;
}

//Mandelstam variables
double Hadron::M(int e)	//valence mass
{
	return pow(ValMass, e);
}

double Hadron::s(int e)	//Probe on Target, fixed by experiment and constant
{
	TLorentzVector S = *Prob + *Targ;
	return pow(S.M2(), e);
}

double Hadron::s_0(int e)	//threshold s in partonic process
{
	return pow(4*M(2), e);
}

double Hadron::s_(int e)	//s in partonic process
{
	//TLorentzVector S = *Part1 + *Part2;
	//return pow(S.M2(), e);
	return pow(fs_, e);
}

double Hadron::t_(int e)	//t in partonic process
{
	//TLorentzVector T = *Part1 - *Quark3;
	//return pow(T.M2(), e);
	return pow (-Q2_(), e);
}

double Hadron::u_(int e)	//u in partonic process
{
	//TLorentzVector U = *Part1 - *Quark4;
	//return pow(U.M2(), e);
	return 2*M(2) - s_() - t_();
}

double Hadron::Q2_(int e)	//momentum transfer in partonic process)
{
	//return -t_(e);
	return pow(fQ2_, e);
}	

double Hadron::CosT_(int e)	//momentum transfer in partonic process)
{
	return pow(fcosT_, e);
}	

void Hadron::SetOmega(double X)	//X is cosTheta
{
	fcosT_ = X;

	double E = sqrt(s_())/2.0;
	double t = M(2) - 2*E*( E - X*sqrt(E*E-M(2)) );
	fQ2_ = -t;
}

void Hadron::SetCMEnergy(double X)
{
	fs_ = X;
}

void Hadron::SetQ2(double X)
{
	fQ2_ = X;

	double E = sqrt(s_())/2.0;
	double cosT = (2*E*E - M(2) - X) / (2*E * sqrt(E*E - M(2)));
	fcosT_ = cosT;
}

double Hadron::Q2Lim(double &Min, double &Max)	//y integration limits
{
	double E = sqrt(s_())/2.0;
	double P = sqrt(E*E - M(2));

	Min = 2*E*( E - P ) - M(2);
	Max = 2*E*( E + P ) - M(2);

	return Max - Min;
}

double Hadron::aS(int e)
{
	//double muR2 = s_();
	double muR2 = M(2) * 1.6*1.6;
	return pow(ProbPdf->alphasQ2(muR2*muR2), e);
}

//XSections
double Hadron::dXSdQ2_qq()	//differential cross section (dXS/dQ2) for qq_ into qq_
{
	if (PSCut())
	{
		double t1 = t_() - M(2);
		double u1 = u_() - M(2);
		return  4 * aS(2) * Const::fPi / (9 * s_(4)) *
			(t1*t1 + u1*u1 + 2*M(2)*s_()) ;
	}
	else
		return 0.0;
}

double Hadron::dXSdQ2_gg()	//differential cross section (dXS/dQ2) for gg into qq_
{
	if (PSCut())
	{
		double t1 = t_() - M(2);
		double u1 = u_() - M(2);
		return  aS(2) * Const::fPi / s_(2) *
			(4.0/3.0 - 3*t1*u1 / s_(2)) / 8.0 *
			(t1/u1 + u1/t1 + 4*M(2)*s_() / (t1*u1) * ( 1 - M(2)*s_() / (t1*u1) ) );
	}
	else
		return 0.0;
}

double Hadron::dXSdOmega_qq()	//differential cross section (dXS/dOmega) for qq_ into qq_
{
	if (PSCut())
	{
		double t1 = t_() - M(2);
		double u1 = u_() - M(2);
		return aS(2) / (9 * s_(3)) *
		       sqrt(1 - 4*M(2)/s_()) *
		       (t1*t1 + u1*u1 + 2*M(2)*s_());
	}
	else
		return 0.0;
}

double Hadron::dXSdOmega_gg()	//differential cross section (dXS/dOmega) for gg into qq_
{
	if (PSCut())
	{
		double t1 = t_() - M(2);
		double u1 = u_() - M(2);
		return aS(2) / (32 * s_()) *
		       sqrt(1 - 4*M(2)/s_()) *
		       (6*t1*u1 / s_(2) - M(2) * (s_()-4*M(2)) / (3.0 * t1*u1) + 
			4.0 * (t1*u1 - 2*M(2)*(2*M(2) + t1)) / (3.0 * t1*t1) +
			4.0 * (t1*u1 - 2*M(2)*(2*M(2) + u1)) / (3.0 * u1*u1) -
			3.0 * (t1*u1 - M(2)*(u1 - t1)) / (s_() * t1) -
			3.0 * (t1*u1 - M(2)*(t1 - u1)) / (s_() * u1) );
	}
	else
		return 0.0;
}

//for vegas integeration
void Hadron::Integrand(const int *nDim, const double *x, const int *nComp, double f[])
{
	double Tau0, y1, y2, Jac;
	switch(IntegrandType)
	{
		case 1:
			Tau0 = s_0()/s();	//\hat{s} / s
			y1 = (1 - Tau0)*x[0] + Tau0;
			y2 = (1 - Tau0/y1)*x[1] + Tau0/y1;
		
			SetCMEnergy(s() * y1 *y2);
		
			SetOmega(2*x[2] - 1);
		
			//jacobian	pdf integration		* phasespace
			Jac = (1 - Tau0) * (1 - Tau0/y1) * 4.0*Const::fPi;
			f[0] = Jac * Variable(y1, y2);
			break;
		case 2:
			Tau0 = s_0()/s();	//\hat{s} / s
			y1 = (1 - Tau0)*x[0] + Tau0;
			y2 = (1 - Tau0/y1)*x[1] + Tau0/y1;
		
			SetCMEnergy(s() * y1 *y2);
		
			//jacobian	pdf integration
			Jac = (1 - Tau0) * (1 - Tau0/y1);
			f[0] = Jac * Variable(y1, y2);
			break;
		case 3:
			SetOmega(2*x[0] - 1);
			f[0] = 4.0*Const::fPi * dXSdQ2_gg();
			break;
		default:
			f[0] = 0.0;
	}
}

double Hadron::Variable(const double x1, const double x2)
{
	//double Mf2 = s_();
	double Mf2 = M(2) * 2.1*2.1;	//SHiP

	if (ProbPdf->inRangeQ2(Mf2) && TargPdf->inRangeQ2(Mf2) && 
	    ProbPdf->inRangeX(x1) && TargPdf->inRangeX(x1) && 
	    ProbPdf->inRangeX(x2) && TargPdf->inRangeX(x2) ) 
	{
		double pdf_gg = ProbPdf->xfxQ2(0, x1, Mf2) * TargPdf->xfxQ2(0, x2, Mf2);
		      pdf_gg += ProbPdf->xfxQ2(0, x2, Mf2) * TargPdf->xfxQ2(0, x1, Mf2);

		double pdf_qq = 0.0;
		for (int i = 1; i < ValQuark; ++i)	//Charm = 4
		{
			pdf_qq += ProbPdf->xfxQ2( i, x1, Mf2) * TargPdf->xfxQ2(-i, x2, Mf2);	//q/p and q_/A
			pdf_qq += ProbPdf->xfxQ2( i, x2, Mf2) * TargPdf->xfxQ2(-i, x1, Mf2);
			pdf_qq += ProbPdf->xfxQ2(-i, x1, Mf2) * TargPdf->xfxQ2( i, x2, Mf2);	//q_/p and q/A
			pdf_qq += ProbPdf->xfxQ2(-i, x2, Mf2) * TargPdf->xfxQ2( i, x1, Mf2);
		}

		if (IntegrandType == 2)
			return (pdf_gg * dXSdQ2_gg() + pdf_qq * dXSdQ2_qq())/(x1*x2);
		else 
			return (pdf_gg * dXSdOmega_gg() + pdf_qq * dXSdOmega_qq())/(x1*x2);
	}
	else
		return 0.0;
}

double Hadron::Total(double &Error, double &Chi2Prob)	//Total xsec
{
	IntegrandType = 1;
	int Trial, Fail;
	return Kine::VegasIntegration(this, 3, Trial, Fail, Error, Chi2Prob);
}

double Hadron::Spectrum(double &Error, double &Chi2Prob)	//Total xsec
{
	IntegrandType = 2;
	int Trial, Fail;
	return Kine::VegasIntegration(this, 2, Trial, Fail, Error, Chi2Prob);
}

double Hadron::XSec(double &Error, double &Chi2Prob)	//Total xsec
{
	IntegrandType = 3;
	int Trial, Fail;
	return Kine::VegasIntegration(this, 1, Trial, Fail, Error, Chi2Prob);
}
