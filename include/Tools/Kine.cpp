#include "Tools/Kine.h"

//Kinematic factors, from Shrock (might be wrong)
double Kine::Unhelicity(double M_Meson, double M_Lepton, double M_Sterile, int H) //Correct scaling of flux
{
	if (H == 0)
		return Unhelicity(M_Meson, M_Lepton, M_Sterile,  1) +
		       Unhelicity(M_Meson, M_Lepton, M_Sterile, -1);
	else if (M_Meson >= M_Lepton + M_Sterile)
	{
		double X = pow(M_Lepton/M_Meson, 2.0);
		double Y = pow(M_Sterile/M_Meson, 2.0);	
		return (X+Y - (X-Y)*(X-Y) + H*(X-Y)*sqrt(Kallen(1, X, Y))) * 
			sqrt(Kallen(1, X, Y)) / (2*X * (1-X)*(1-X));
	}
	else return 0;
}

//kallen triangular function, it should never be negative
double Kine::Kallen(double X, double Y, double Z)
{
	double Lambda = X*X + Y*Y + Z*Z - 2*(X*Y + X*Z + Y*Z);
	if (Lambda < 0) return -Lambda;
	else return Lambda;
}


//Bethe blocke radiatio nformula
double Kine::Bethe(double Beta, double Mass, double Density, double I, int Z, int A)
{
	double K = 0.307075;	//From PDG MeV mol-1 cm2
	double Beta2 = Beta*Beta;
	double Gamma = 1.0/sqrt(1-Beta2);
	double Gamma2 = Gamma*Gamma;
	double e2M = Const::fMElectron / Mass;
	double Wmax = (2000 * Const::fMElectron * Beta2 * Gamma2) / (1 + 2*Gamma*e2M + e2M*e2M);
	double LogArg = 2000 * Const::fMElectron * Beta2 * Gamma2 * Wmax / (1e-12*I*I); 	//Everything in MeV
	return 0.1 * Density * (K * Z) / (A * Beta2) * (0.5 * log (LogArg) - Beta2);	//stopping power in GeV/m (0.1*)
}

//raiantion length
double Kine::RadiationLength(double Density, int Z, int A)
{
	double K = 716.408;
	double L0 = log(184.15 * pow(Z, -1.0/3.0));
	double L1 = log(1194 * pow(Z, -2.0/3.0));

	return 0.01 * Density * K*A / (Z*Z*(L0 - Rad(Z)) + Z*L1);		//in m (0.01*)
}

double Kine::Rad(int Z)
{
	double a = Const::fAem * Z;
	return a*a * ( 1.0/(1+a*a) + 0.20206 - 0.0369*a*a + 0.0083*pow(a,4) - 0.002*pow(a,6));
}
