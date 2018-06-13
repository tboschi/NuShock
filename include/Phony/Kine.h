/*
 * Tools/Kine
 * namespace Kine for kinematic functions
 * 
 * Author: Tommaso Boschi
 */

#ifndef TOOLS_H
#define TOOLS_H

#include <cmath>

//Kinematic functions
namespace Kine
{
	//unhelicity
	double Unhelicity(double M_Meson, double M_Lepton, double M_Sterile, int H = 0); 
	double Kallen(double X, double Y, double Z);

	//radiation stuff
	double Bethe(double Beta, double Mass, double Density, double I, int Z, int A);		//GeV/m
	double RadiationLength(double Density, int Z, int A);
	double Rad(int Z);

	//must be moved to decay, here it is just useless
	double I1_f(double t, double X, double Y, double Z);	//To be integrated
	double I1_xyz(double X, double Y, double Z);
	double I1_xy(double X, double Y);
	double I2_f(double t, double X, double Y, double Z);	//To be integrated
	double I2_xyz(double X, double Y, double Z);
	double I2_xy(double X, double Y);
	double I3_xy(double X, double Y);

}

#endif
