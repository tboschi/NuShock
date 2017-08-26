/*
 * Tools for MC
 * Here namespaces are defined, continaing
 * Kine: kinematic formualae
 * Const: standard model constants
 *
 * Author: Tommaso Boschi
 */

#ifndef TOOLS_H
#define TOOLS_H

#include <iostream>
#include <fstream>
#include <cstring>
#include <map>
#include <vector>

//ROOT include
#include "TMath.h"

//GENIE include
#include "Conventions/Constants.h"
#include "PDG/PDGLibrary.h"

//Channel modes
enum ChannelName
{
	_undefined,
	_ALL,		//
	_nnn,		//3 body
	_nGAMMA,	//2 body
	_nEE,		//3 body
	_nEMU,		//3 body
	_nMUE,		//3 body
	_nPI0,		//3 body
	_EPI,		//2 body
	_nMUMU,		//3 body
	_MUPI,		//2 body
	_EKA,		//2 body
	_nKA0,		//3 body
	_Muon,		//3 body
	_Kaon,		//3 body
	_Kaon0		//3 body
};

//Kinematic functions
namespace Kine
{
	double ShrockFactor(double M_Meson, double M_Lepton, double M_Sterile); 
	double ShrockRho(double X, double Y);
	double ShrockFM(double X, double Y);
	double Lambda(double X, double Y, double Z);
	double I1_f(double t, double X, double Y, double Z);	//To be integrated
	double I1_xyz(double X, double Y, double Z);
	double I1_xy(double X, double Y);
	double I2_f(double t, double X, double Y, double Z);	//To be integrated
	double I2_xyz(double X, double Y, double Z);
	double I2_xy(double X, double Y);
	double I3_xy(double X, double Y);

	double Bethe(double Beta, double Mass, double Density, double I, int Z, int A);		//GeV/m
	double RadiationLength(double Density, int Z, int A);
	double Rad(int Z);

	static const int Sample = 1000;	//Sample for integration
	static const int Loop = 1000;	//Sample for integration
}

//Constants
namespace Const
{
	static const double fC = 299792458;		//m/s
	static const double fhBar = 6.5821189916e-25;	//GeV s, from PDG
	static const double fAem = 1.0/137.035999074;	// from PDG

	//Conversion
	static const double fM2GeV = 5.06e15;	//1GeV in 1/m
	static const double fS2GeV = 1.52e24;	//1GeV in 1/s
	static const double fPi = 3.1415926536;	//1GeV in 1/s
	static const double fPi2 = fPi*fPi;	//1GeV in 1/s
	static const double fPi3 = fPi2*fPi;	//1GeV in 1/s
	static const double fDeg = 180.0/fPi;	//Rad to Deg

	//CKM entries
	static const double fV_ud = 0.97417;
	static const double fV_us = 0.2248;
	static const double fV_ub = 0.0409;
	static const double fV_cd = 0.220;
	static const double fV_cs = 0.995;
	static const double fV_cb = 0.0405;
	static const double fV_td = 0.0082;
	static const double fV_ts = 0.0400;
	static const double fV_tb = 1.009;

	//PMNS entries
	static const double fU_e1 = 0.81;
	static const double fU_e2 = 0.54;
	static const double fU_e3 = -0.15;
	static const double fU_m1 = -0.35;
	static const double fU_m2 = 0.70;
	static const double fU_m3 = 0.62;
	static const double fU_t1 = 0.44;
	static const double fU_t2 = -0.45;
	static const double fU_t3 = 0.77;

	//Masses
	static const double fMElectron = genie::PDGLibrary::Instance()->Find(11)->Mass();
	static const double fMMuon = genie::PDGLibrary::Instance()->Find(13)->Mass();
	static const double fMPion = genie::PDGLibrary::Instance()->Find(211)->Mass();
	static const double fMPion0 = genie::PDGLibrary::Instance()->Find(111)->Mass();
	static const double fMKaon = genie::PDGLibrary::Instance()->Find(321)->Mass();
	static const double fMKaon0 = genie::PDGLibrary::Instance()->Find(311)->Mass();
	static const double fMProton = genie::PDGLibrary::Instance()->Find(2212)->Mass();
	static const double fMNeutron = genie::PDGLibrary::Instance()->Find(2211)->Mass();

	//SM constant - PDG 2016
	static const double fGF = 1.16637876e-5;	//GeV-2, from PDG
	static const double fGF2 = fGF*fGF;		//From PDG
	static const double fSin2W = 0.23129;		//Sin weinberg squared - MSbar scheme
	static const double fFPion2 = pow(0.1302, 2.0);	//Decay constant squared of pion
	static const double fFKaon2 = pow(0.1556, 2.0);	//Decay constant squared of kaon
	static const double fLambdaPlus = 0.0297;		//Linear dependence of f+ in Ke3 (PDG)
	static const double fLambdaZero = 0.0196;		//Linear dependence of f0 in Km3 (PDG)
	static const double fVusFKaon = 0.2165;		//From 1607.00299
}

#endif
