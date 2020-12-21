/*
 * Tools
 * Const: standard model constants
 *
 * Author: Tommaso Boschi
 */

#ifndef CONST_H
#define CONST_H

#include <cmath>

//Constants
namespace Const
{
	static const double C = 299792458;		//m/s
	static const double hBar = 6.5821189916e-25;	//GeV s, from PDG
	static const double Aem = 1.0/137.035999074;	// from PDG
	static const double Na = 6.02214085774e23;	//mol-1

	//Conversion
	static const double M2GeV = 5.06e15;		//1GeV in 1/m
	static const double S2GeV = 1.52e24;		//1GeV in 1/s
	static const double GeV2cm = 389.4e-30;	//1GeV-2 in cm2
	static const double GeV2ub = 0.3894e3;		//1GeV-2 in ub
	static const double pi = 3.1415926536;		//pi
	static const double pi2 = pi*pi;		//pi
	static const double pi3 = pi2*pi;		//pi
	static const double pi4 = pi3*pi;		//pi
	static const double pi5 = pi4*pi;		//pi
	static const double Deg = 180.0/pi;		//Rad to Deg

	//CKM entries
	static const double U_ud = 0.97417;
	static const double U_us = 0.2248;
	static const double U_ub = 0.0409;
	static const double U_cd = 0.220;
	static const double U_cs = 0.995;
	static const double U_cb = 0.0405;
	static const double U_td = 0.0082;
	static const double U_ts = 0.0400;
	static const double U_tb = 1.009;

	//PMNS entries
	static const double U_e1 = 0.81;
	static const double U_e2 = 0.54;
	static const double U_e3 = -0.15;
	static const double U_m1 = -0.35;
	static const double U_m2 = 0.70;
	static const double U_m3 = 0.62;
	static const double U_t1 = 0.44;
	static const double U_t2 = -0.45;
	static const double U_t3 = 0.77;

	//Masses in GeV - PDG 2017
	static const double MQuarkU   = 2.2e-3;
	static const double MQuarkD   = 4.7e-3;
	static const double MQuarkS   = 96e-3;
	static const double MQuarkC   = 1.27;
	static const double MQuarkB   = 4.18;
	static const double MQuarkT   = 173.1;
	static const double MPhoton   = 0.0;
	static const double MNeutrino = 0.0;
	static const double MElectron = 0.510999e-3;
	static const double MMuon     = 105.6583e-3;
	static const double MTau      = 1776.86e-3;
	static const double MPion     = 139.57061e-3;
	static const double MPion0    = 134.9770e-3;
	static const double MKaon     = 493.677e-3;
	static const double MKaon0    = 497.611e-3;
	static const double MEta      = 547.862e-3;
	static const double MEtai     = 957.78e-3;
	static const double MRho      = 775.11e-3;
	static const double MRho0     = 775.26e-3;
	static const double MOmega    = 782.65e-3;
	static const double MKaonx    = 891.76e-3;
	static const double MKaon0x   = 682e-3;
	static const double MPhi      = 1019.460e-3;
	static const double MD	      = 1869.56e-3;
	static const double MDs       = 1968.28e-3;
	static const double MProton   = 938.272081e-3;
	static const double MNeutron  = 939.565143e-3;
	static const double MW        = 80.385;
	static const double MZ        = 91.1876;

	//SM constant - PDG 2016
	static const double GF = 1.16637876e-5;	//GeV-2, from PDG
	static const double GF2 = GF*GF;		//From PDG
	static const double sin2W = 0.23129;		//sin weinberg squared - MSbar scheme
	
	//decay constant - 0901.3789
	/*
	static const double DPion   = 0.1327;
	static const double DPion0  = 0.1300;
	static const double DKaon   = 0.1598;
	static const double DRho    = 0.2200;
	static const double DRho0   = 0.2200;
	static const double DEta    = 0.1647;
	static const double DEtai   = 0.1529;
	static const double DOmega  = 0.1950;
	static const double DKaonx  = 0.2170;
	static const double DKaon0x = 0.2170;	
	static const double DPhi    = 0.2290;
	static const double DCharm  = 0.2226;
	*/

	//decay constant - 1805.08567
	static const double DPion   = 0.1302;		//GeV
	//static const double DPion0  = 0.1300;	//GeV
	static const double DKaon   = 0.1556;		//GeV
	static const double DRho    = 0.1620;		//GeV²
	static const double VRho    = 1 - 2 * sin2W;
	//static const double DRho0   = 0.2200;	//GeV²
	static const double DEta    = 0.0817;		//GeV
	static const double DEtai   = 0.0947;		//GeV
	static const double DOmega  = 0.1530;		//GeV²
	static const double VOmega  = 4 * sin2W / 3.0;
	static const double DKaonx  = 0.1933;	//maybe not possible
	//static const double DKaon0x = 0.2170;	
	static const double DPhi    = 0.2340;		//GeV²
	static const double VPhi    = 4 * sin2W / 3.0 - 1;
	static const double DCharm  = 0.2120;		//GeV

	static const double KaPi = 0.9700;	//f(0)
	static const double K0L_ = 0.0267;	//Linear dependence of f+ in K0m3 (PDG)
	static const double K0L0 = 0.0117;	//Linear dependence of f0 in K0m3 (PDG)
	static const double KCL_ = 0.0277;	//Linear dependence of f+ in K+m3 (PDG)
	static const double KCL0 = 0.0183;	//Linear dependence of f0 in K+m3 (PDG)

	static const double MagMuN = -1.9130427345;	//neutron magnetic moment (in nuclear magneton units)
	static const double MagMuP = 2.79284735128;	//proton magnetic moment (in nuclear magneton units);
	static const double MA = 0.990;		//GeV, axial mass, from GENIE
	//static const double MA = 1.032;		//GeV, axial mass, from PRD35, 785 (1987) [Arhens]
	//static const double MA = 1.270;		//GeV, axial mass, from PRD92, 113011 (2015) [lattice]
	//static const double MA = 1.026;		//GeV, axial mass, from Giunti-Kim
	static const double MV = 0.840;		//GeV, vectorial mass, from GENIE
	static const double GA0 = 1.2671;		//Axial form factor at Q2 = 0, from GENIE

	static const double Wwidth = 2.085;		//W decay width, in GeV
	static const double Zwidth = 2.4952;		//W decay width, in GeV
}

#endif
