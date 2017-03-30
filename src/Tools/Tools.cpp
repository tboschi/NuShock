#include "Tools.h"

//Kinematic factors
double Kine::ShrockFactor(double M_Meson, double M_Lepton, double M_Sterile) //Correct scaling of flux
{
	double dM_a = pow(M_Lepton/M_Meson, 2.0);
	double dM_i = pow(M_Sterile/M_Meson, 2.0);	

	if (M_Meson >= M_Lepton + M_Sterile)
		return ShrockRho(dM_a, dM_i)/(dM_a * pow(1-dM_a, 2.0));
	else return 0;
}

double Kine::ShrockRho(double X, double Y)
{
	return ShrockFM(X, Y)*sqrt(ShrockLambda(1, X, Y));
}

double Kine::ShrockFM(double X, double Y)
{
	return X+Y - (X-Y)*(X-Y);
}

double Kine::ShrockLambda(double X, double Y, double Z)
{
	double Lambda = X*X + Y*Y + Z*Z - 2*(X*Y + X*Z + Y*Z);
	if (Lambda < 0) return -Lambda;
	else return Lambda;
}

double Kine::I1_f(double t, double X, double Y, double Z)	//To be integrated
{
	return (t-X*X-Y*Y)*(1+Z*Z-t)*sqrt(ShrockLambda(t, X*X, Y*Y))*sqrt(ShrockLambda(1, t, Z*Z)) / t;
}

double Kine::I1_xyz(double X, double Y, double Z)
{
	double A = (1-Z)*(1-Z);
	double B = (X+Y)*(X+Y);
	if (A > B)
	{
		double tmp = B;
		B = A;
		A = tmp;
	}
	double a, b;
	double h = (B-A)/Sample;
	double Integral = 0;	//3/8 Simpson's method for integration
	for (a = A; b < B; a = b)
	{
		b = a + h;
		Integral += 0.125 * (I1_f(a, X, Y, Z) + 3*I1_f((a+2*b)/3.0, X, Y, Z) + 3*I1_f((2*a+b)/3.0, X, Y, Z) + I1_f(b, X, Y, Z));
	}	
	return 12.0 * Integral;
}

double Kine::I1_xy(double X, double Y)
{
	return ((1+X-Y)*(1+X) - 4*X) * sqrt(ShrockLambda(1.0, X, Y));
}

double Kine::I2_f(double t, double X, double Y, double Z)	//To be integrated
{
	return (1+X*X-t)*sqrt(ShrockLambda(t, Y*Y, Z*Z))*sqrt(ShrockLambda(1.0, t, X*X)) / t;
}

double Kine::I2_xyz(double X, double Y, double Z)
{
	double A = (1-X)*(1-X);
	double B = (Y+Z)*(Y+Z);
	if (A > B)
	{
		double tmp = B;
		B = A;
		A = tmp;
	}
	double a, b;
	double h = (B-A)/Sample;
	double Integral = 0.0;	//3/8 Simpson's method for integration
	for (a = A; b < B; a = b)
	{
		b = a + h;
		Integral += 0.125 * (I2_f(a, X, Y, Z) + 3*I2_f((a+2*b)/3.0, X, Y, Z) + 3*I2_f((2*a+b)/3.0, X, Y, Z) + I2_f(b, X, Y, Z));
	}	
	return 24.0*Y*Z * Integral;
}

double Kine::I2_xy(double X, double Y)
{
	return ((1+X-Y)*(1+X+2*Y) - 4*X) * sqrt(ShrockLambda(1.0, X, Y));
}

double Kine::I3_xy(double X, double Y)
{
	return (1+2*Y)*(1-Y) * sqrt(ShrockLambda(1.0, X, Y));
}
