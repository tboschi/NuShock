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
#include <sstream>
#include <iomanip>
#include <getopt.h>
#include <algorithm>

//ROOT include
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TRandom3.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TDatabasePDG.h"
#include "TParticlePDG.h"

//CUBA
#include "cuba.h"

//Channel modes
enum ChannelName
{
	_undefined,
	_ALL,		//
	_nnn,		//3 body	2
	_nGAMMA,	//2 body	3
	_nEE,		//3 body	4
	_nEMU,		//3 body	5
	_nMUE,		//3 body	6
	_nMUMU,		//3 body	7
	_nET,		//3 body	8
	_nTE,		//3 body	9
	_nMUT,		//3 body	10
	_nTMU,		//3 body	11
	_nPI0,		//2 body	12
	_EPI,		//2 body	13
	_MUPI,		//2 body	14
	_TPI,		//2 body	15
	_EKA,		//2 body	16
	_MUKA,		//2 body	17
	_nRHO0,		//2 body	18
	_ERHO,		//2 body	19
	_MURHO,		//2 body	20
	_EKAx,		//2 body	21
	_nKA0x,		//2 body	22
	_MUKAx,		//2 body	23
	_nOMEGA,	//2 body	24
	_nETA,		//2 body	25
	_nETAi,		//2 body	26
	_nPHI,		//2 body	27
	_ECHARM,	//2 body	28
	//special
	_ExpALL,	//
////////////////////////////////
	_Muon,		//3 body
	_TauE,		//3 body
	_TauM,		//3 body
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


	template<class TempClass>
	int Integrand(const int *nDim, const double x[], const int *nComp, double f[], void *UserData)
	{
		TempClass *TempObject = static_cast<TempClass*>(UserData);

		TempObject->Integrand(nDim, x, nComp, f);
	}

	template<class TempClass>
	double VegasIntegration(TempClass *TempObject, int nDim, int &Trial, int &Fail, double &Error, double &Chi2Prob)
	{
		//input
		double EpsRel = 1.0e-8;		//relative error for each component
		double EpsAbs = 1.0e-12;	//absolute error
		int MinEval = 1e5;		//minimum number of evaluation
		int MaxEval = 1e6;		//maximum number of evaluation
		int nStart = 1000;
		int nIncrease = 500;
		int nBatch = 1000;
		void *UserData = TempObject;
		char *state = NULL;
		void *spin = NULL;

		integrand_t IntCast = reinterpret_cast<integrand_t>(&Integrand<TempClass>); 
	
		//output
		double Integral;

		Vegas(nDim, 1, IntCast, UserData, 1, 	//ndim, ncomp, integrand_t, userdata, nvec
		      EpsRel, EpsAbs, 0, 0, 		//epsrel, epsabs, verbosity, seed
		      MinEval, MaxEval, nStart, nIncrease, nBatch,
		      0, state, spin,			//gridno, statefile, spin
		      &Trial, &Fail, &Integral, &Error, &Chi2Prob);
	
		return Integral;
	}

	template<class TempClass>
	double BooleIntegration(TempClass *TempObject)
	{
		double a = 0, b = 0;
		double h = 1.0/100.0;
		double Integral = 0;
		for (a = 0; b + 1e-12 < 1.0; a = b)
		{
			b = a + h;

			Integral += 7  * TempObject->Variable( a );
			Integral += 32 * TempObject->Variable( (3*a+b)/4.0 );
			Integral += 12 * TempObject->Variable( (a+b)/2.0 );
			Integral += 32 * TempObject->Variable( (a+3*b)/4.0 );
			Integral += 7  * TempObject->Variable( b );
		}	
	
		return Integral * h/90.0;
	}	
}

//Constants
namespace Const
{
	static const double fC = 299792458;		//m/s
	static const double fhBar = 6.5821189916e-25;	//GeV s, from PDG
	static const double fAem = 1.0/137.035999074;	// from PDG
	static const double fNa = 6.02214085774e23;	//mol-1

	//Conversion
	static const double fM2GeV = 5.06e15;		//1GeV in 1/m
	static const double fS2GeV = 1.52e24;		//1GeV in 1/s
	static const double fGeV2cm = 389.4e-30;	//1GeV-2 in cm2
	static const double fGeV2ub = 0.3894e3;		//1GeV-2 in ub
	static const double fPi = 3.1415926536;		//pi
	static const double fPi2 = fPi*fPi;		//pi
	static const double fPi3 = fPi2*fPi;		//pi
	static const double fDeg = 180.0/fPi;		//Rad to Deg

	//CKM entries
	static const double fU_ud = 0.97417;
	static const double fU_us = 0.2248;
	static const double fU_ub = 0.0409;
	static const double fU_cd = 0.220;
	static const double fU_cs = 0.995;
	static const double fU_cb = 0.0405;
	static const double fU_td = 0.0082;
	static const double fU_ts = 0.0400;
	static const double fU_tb = 1.009;

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

	//Masses in GeV - PDG 2017 -- remove genie dep on Tools
	static const double fMQuarkU = 2.2e-3;
	static const double fMQuarkD = 4.7e-3;
	static const double fMQuarkS = 96e-3;
	static const double fMQuarkC = 1.28;
	static const double fMQuarkB = 4.18;
	static const double fMQuarkT = 173.1;
	static const double fMElectron = 0.510999e-3;
	static const double fMMuon = 105.6583e-3;
	static const double fMTau = 1776.86e-3;
	static const double fMPion = 139.57061e-3;
	static const double fMPion0 = 134.9770e-3;
	static const double fMKaon = 493.677e-3;
	static const double fMKaon0 = 497.611e-3;
	static const double fMEta = 547.862e-3;
	static const double fMEtai = 957.78e-3;
	static const double fMRho = 775.11e-3;
	static const double fMRho0 = 775.26e-3;
	static const double fMOmega = 782.65e-3;
	static const double fMKaonx = 891.76e-3;
	static const double fMKaon0x = 682e-3;
	static const double fMPhi = 1019.460e-3;
	static const double fMCharm = 1869.56e-3;
	static const double fMDs = 1968.28e-3;
	static const double fMProton = 938.272081e-3;
	static const double fMNeutron = 939.565143e-3;
	static const double fMW = 80.385;
	static const double fMZ = 91.1876;

	//SM constant - PDG 2016
	static const double fGF = 1.16637876e-5;	//GeV-2, from PDG
	static const double fGF2 = fGF*fGF;		//From PDG
	static const double fSin2W = 0.23129;		//Sin weinberg squared - MSbar scheme
	
	//decay constant
	static const double fDPion2   = pow(0.1307, 2.0);
	static const double fDPion02  = pow(0.1300, 2.0);
	static const double fDKaon2   = pow(0.1598, 2.0);
	static const double fDRho2    = pow(0.2200, 2.0);
	static const double fDRho02   = pow(0.2200, 2.0);
	static const double fDEta2    = pow(0.1647, 2.0);
	static const double fDEtai2   = pow(0.1529, 2.0);
	static const double fDOmega2  = pow(0.1950, 2.0);
	static const double fDKaonx2  = pow(0.2170, 2.0);
	static const double fDKaon0x2 = pow(0.2170, 2.0);	
	static const double fDPhi2    = pow(0.2290, 2.0);
	static const double fDCharm2  = pow(0.2226, 2.0);

	static const double fVLight   = fSin2W / 3.0;		//Vector constant for Rho 0 and Omega
	static const double fVStrange = -0.25+fSin2W/3.0;	//Vector constant for Kaon star and Phi
	static const double fVCharm   = 0.25-2.0*fSin2W / 3.0;	//Vector constant for Charm meson D and JPsi
	static const double fVusFKaon = 0.2165;			//From 1607.00299

	static const double fLambdaPlus = 0.0297;	//Linear dependence of f+ in Ke3 (PDG)
	static const double fLambdaZero = 0.0196;	//Linear dependence of f0 in Km3 (PDG)
	static const double fMagMuN = -1.9130427345;	//neutron magnetic moment (in nuclear magneton units)
	static const double fMagMuP = 2.79284735128;	//proton magnetic moment (in nuclear magneton units);
	static const double fMA = 0.990;		//GeV, axial mass, from GENIE
	//static const double fMA = 1.032;		//GeV, axial mass, from PRD35, 785 (1987) [Arhens]
	//static const double fMA = 1.270;		//GeV, axial mass, from PRD92, 113011 (2015) [lattice]
	//static const double fMA = 1.026;		//GeV, axial mass, from Giunti-Kim
	static const double fMV = 0.840;		//GeV, vectorial mass, from GENIE
	static const double fGA0 = 1.2671;		//Axial form factor at Q2 = 0, from GENIE

	static const double fWwidth = 2.085;		//W decay width, in GeV
	static const double fZwidth = 2.4952;		//W decay width, in GeV
}

#endif
