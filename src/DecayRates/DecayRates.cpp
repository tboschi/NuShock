#include "DecayRates.h"

//M_Sterile > 2 M_Electron (always)
double DecayRates::nEE(double M_Sterile, double U_e, double U_m, double U_t)
{
	double dMe = M_Electron / M_Sterile;
	double dMn = M_Neutrino / M_Sterile;
	double gL = -0.5 + Tools::Const::fSin2w;
	double gR = Tools::Const::fSin2w;
	double Int1 = I1_xyz(dMn, dMe, dMe);
	double Int2 = I2_xyz(dMn, dMe, dMe);
	double KF_e = (gL*gR + gR) * Int2 + (gL*gL + gR*gR + (1+2*gL))*Int1;
	double KF_m = (gL*gR) * Int2 + (gL*gL + gR*gR)*Int1;
//	double KF_t = (gL*gR) * I2_xyz(dMn, dMe, dMe) + (gL*gL + gR*gR)*I1_xyz(dMn, dMe, dMe);

	return genie::constants::kGF2 * pow(M_Sterile, 5) * (U_e*U_e * KF_e + (U_m*U_m + U_t*U_t) * KF_m) / (96.0 * genie::constants::kPi3);
}

double DecayRates::nGAMMA(double M_Sterile, double U_e, double U_m, double U_t)
{
	double AemPi = genie::constants::kAem / genie::constants::kPi;

	return genie::constants::kGF2 * pow(M_Sterile, 5) / (192.0 * genie::constants::kPi3) * (U_e*U_e + U_m*U_m + U_t*U_t) * (27.0/32.0 * AemPi);
}

//M_Sterile > M_Muon
double DecayRates::EMU(double M_Sterile, double U_e, double U_m, double U_t)
{
	double dMe = M_Electron / M_Sterile;
	double dMm = M_Muon / M_Sterile;
	double dMn = M_Neutrino / M_Sterile;

	return genie::constants::kGF2 * pow(M_Sterile, 5) * (U_e*U_e * I1_xyz(dMm, dMn, dMe) + U_m*U_mi * I1_xyz(dMe, dMn, dMm)) / (192.0 * genie::constants::kPi3);
}

//M_Sterile > M_Pion
double DecayRates::EPI(double M_Sterile, double U_e, double U_m, double U_t)
{
	double dMe2 = M_Electron*M_Electron/M_Sterile/M_Sterile;
	double dMp2 = M_Pion*M_Pion/M_Sterile/M_Sterile;

	return genie::constants::kGF2 * pow(M_Sterile, 3) / (16.0 * genie::constants::kPi) * U_e*U_e * pow(Tools::Const::fV_ud, 2.0) * Tools::Const:fFPion2 * I1_xy(dMe2, dMp2); 
}

double DecayRates::nPI0(double M_Sterile, double U_e, double U_m, double U_t)
{
	double dMp2 = M_Pion0*M_Pion0/M_Sterile/M_Sterile;

	return genie::constants::kGF2 * pow(M_Sterile, 3) / (64.0 * genie::constants::kPi) * (U_e*U_e + U_m*U_m + U_t*U_t) * pow((1.0-dMp2), 2.0);
}

//M_Sterile > 2 M_Muon
double DecayRates::nMUMU(double M_Sterile, double U_e, double U_m, double U_t)
{
	double dMm = M_Muon / M_Sterile;
	double dMn = M_Neutrino / M_Sterile;
	double gL = -0.5 + Tools::Const::fSin2w;
	double gR = Tools::Const::fSin2w;
	double Int1 = I1_xyz(dMn, dMm, dMm);
	double Int2 = I2_xyz(dMn, dMm, dMm);
	double KF_m = (gL*gR + gR) * Int2 + (gL*gL + gR*gR + (1+2*gL))*Int1;
	double KF_e = (gL*gR) * Int2 + (gL*gL + gR*gR)*Int1;
//	double KF_t = (gL*gR) * I2_xyz(dMn, dMe, dMe) + (gL*gL + gR*gR)*I1_xyz(dMn, dMe, dMe);

	return genie::constants::kGF2 * pow(M_Sterile, 5) * (U_m*U_m * KF_m + (U_e*U_e + U_t*U_t) * KF_e) / (96.0 * genie::constants::kPi3);
}

//M_Sterile > M_Pion + M_Muon
double DecayRates::MUPI(double M_Sterile, double U_e, double U_m, double U_t)
{
	double dMm2 = M_Muon*M_Muon/M_Sterile/M_Sterile;
	double dMp2 = M_Pion*M_Pion/M_Sterile/M_Sterile;

	return genie::constants::kGF2 * pow(M_Sterile, 3) / (16.0 * genie::constants::kPi) * U_m*U_m * pow(Tools::Const::fV_ud, 2.0) * Tools::Const::fFPion2 * I1_xy(dMm2, dMp2); 
}
