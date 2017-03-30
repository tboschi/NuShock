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
#include "TH1.h"
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"

//GENIE include
#include "Conventions/Constants.h"
#include "PDG/PDGLibrary.h"

//Kinematic functions
namespace Kine
{
	double ShrockFactor(double M_Meson, double M_Lepton, double M_Sterile); 
	double ShrockRho(double X, double Y);
	double ShrockFM(double X, double Y);
	double ShrockLambda(double X, double Y, double Z);
	double I1_f(double t, double X, double Y, double Z);	//To be integrated
	double I1_xyz(double X, double Y, double Z);
	double I1_xy(double X, double Y);
	double I2_f(double t, double X, double Y, double Z);	//To be integrated
	double I2_xyz(double X, double Y, double Z);
	double I2_xy(double X, double Y);
	double I3_xy(double X, double Y);
	static double Sample = 100.0;
}

//Constants
namespace Const
{
	//Conversion
	static double fM2GeV = 5.06e15;

	//CKM entries
	static double fV_ud = 0.97417;
	static double fV_us = 0.2248;
	static double fV_ub = 0.0409;
	static double fV_cd = 0.220;
	static double fV_cs = 0.995;
	static double fV_cb = 0.0405;
	static double fV_td = 0.0082;
	static double fV_ts = 0.0400;
	static double fV_tb = 1.009;

	//PMNS entries
	static double fU_e1 = 0.81;
	static double fU_e2 = 0.54;
	static double fU_e3 = -0.15;
	static double fU_m1 = -0.35;
	static double fU_m2 = 0.70;
	static double fU_m3 = 0.62;
	static double fU_t1 = 0.44;
	static double fU_t2 = -0.45;
	static double fU_t3 = 0.77;

	//Masses
	static double fMElectron = genie::PDGLibrary::Instance()->Find(11)->Mass();
	static double fMMuon = genie::PDGLibrary::Instance()->Find(13)->Mass();
	static double fMPion = genie::PDGLibrary::Instance()->Find(211)->Mass();
	static double fMPion0 = genie::PDGLibrary::Instance()->Find(111)->Mass();
	static double fMKaon = genie::PDGLibrary::Instance()->Find(321)->Mass();
	static double fMKaon0 = genie::PDGLibrary::Instance()->Find(311)->Mass();

	//SM constant - PDG 2016
	static double fSin2W = 0.23129;		//Sin weinberg squared - MSbar scheme
	static double fFPion2 = pow(0.1302, 2.0);	//Decay constant squared of pion
	static double fFKaon2 = pow(0.1556, 2.0);	//Decay constant squared of kaon
}

#endif
