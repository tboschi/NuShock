#include "DecayRates.h"

Decay::Decay(double MSterile, double Ue, double Um, double Ut)	//Decay rates calculator
{
	SetMSterile(MSterile);
	SetUe(Ue);
	SetUm(Um);
	Set(Ut);

	MapInit();
}

//Initialisation of maps
void Decay::MapInit()
{
	mapChannel["ALL"] = ev_ALL;
	mapChannel["nnn"] = ev_nnn;
	mapChannel["nGAMMA"] = ev_nGAMMA;
	mapChannel["nEE"] = ev_nEE;
	mapChannel["nEMU"] = ev_nEMU;
	mapChannel["nPI0"] = ev_nPI0;
	mapChannel["EPI"] = ev_EPI;
	mapChannel["MUPI"] = ev_MUPI;
	mapChannel["nMUMU"] = ev_nMUMU;
	mapChannel["EKA"] = ev_EKA;
	mapChannel["nKA0"] = ev_nKA0;

	mapEnhance["nnn"] = 1.0;
	mapEnhance["nGAMMA"] = 1.0;
	mapEnhance["nEE"] = 1.0;
	mapEnhance["nEMU"] = 1.0;
	mapEnhance["nPI0"] = 1.0;
	mapEnhance["EPI"] = 1.0;
	mapEnhance["MUPI"] = 1.0;
	mapEnhance["nMUMU"] = 1.0;
	mapEnhance["EKA"] = 1.0;
	mapEnhance["nKA0"] = 1.0;
}

double Decay::Gamma(std::string Channel, double B = 1.0)
{
	SetEnhancement(Channel, B);

	double Result;
	switch(mapChannel[Channel])
	{
		case ev_ALL:
			Result = Total();
			break;
		case ev_nnn:
			Result = nnn();
			break;
		case ev_nGAMMA:
			Result = nGAMMA();
			break;
		case ev_nEE:
			Result = nEE();
			break;
		case ev_nEMU:
			Result = nEMU();
			break;
		case ev_nPI0:
			Result = nPI0();
			break;
		case ev_EPI:
			Result = EPI();
			break;
		case ev_MUPI:
			Result = MUPI();
			break;
		case ev_nMUMu:
			Result = nMUMU();
			break;
		case ev_EKA:
			Result = EKA();
			break;
		case ev_nKA0:
			Result = nKA0();
			break;
		default:
			Result = -1.0;
			break;
	}

	SetEnhancement();
	return Result;
}

double Decay::Other(std::string Channel, double A)
{
	SetEnhancement("ALL", A);

	double Result;
	switch(mapChannel[Channel])
	{
		case ev_ALL:
			Result = 0.0;
			break;
		case ev_nnn:
			Result = Total()-nnn();
			break;
		case ev_nGAMMA:
			Result = Total()-nGAMMA();
			break;
		case ev_nEE:
			Result = Total()-nEE();
			break;
		case ev_nEMU:
			Result = Total()-nEMU();
			break;
		case ev_nPI0:
			Result = Total()-nPI0();
			break;
		case ev_EPI:
			Result = Total()-EPI();
			break;
		case ev_MUPI:
			Result = Total()-MUPI();
			break;
		case ev_nMUMu:
			Result = Total()-nMUMU();
			break;
		case ev_EKA:
			Result = Total()-EKA();
			break;
		case ev_nKA0:
			Result = Total()-nKA0();
			break;
		default:
			Result = -1.0;
			break;
	}

	SetEnhancement();
	return Result;
}

double Decay::Branch(std::string Channel, double A, double B)
{
	SetEnhancement("ALL", A);
	SetEnhancement(Channel, B);
	double TotalGamma = Total();

	double Result;
	switch(mapChannel[Channel])
	{
		case ev_ALL:
			Result = 1.0;
			break;
		case ev_nnn:
			Result = nnn()/Total();
			break;
		case ev_nGAMMA:
			Result = nGAMMA()/Total();
			break;
		case ev_nEE:
			Result = nEE()/Total();
			break;
		case ev_nEMU:
			Result = nEMU()/Total();
			break;
		case ev_nPI0:
			Result = nPI0()/Total();
			break;
		case ev_EPI:
			Result = EPI()/Total();
			break;
		case ev_MUPI:
			Result = MUPI()/Total();
			break;
		case ev_nMUMu:
			Result = nMUMU()/Total();
			break;
		case ev_EKA:
			Result = EKA()/Total();
			break;
		case ev_nKA0:
			Result = nKA0()/Total();
			break;
		default:
			Result = -1.0;
			break;
	}

	SetEnhancement();
	return Result;
}

//Controller of decay enhancement
void SetEnhancement(std::string Channel, double K)
{
	std::map<std::string, double>::iterator it;
	if (Channel == "ALL")
		for (it = mapEnhance.begin(); it != mapEnhance.end(); ++it)
			it->second = K;
	else
	{
		it = mapEnhance.find(Channel);
		if (it != mapEnhance.end())
			mapEnhance[Channel] = K;
	}
}

//Very boring stuff
//total decay width
double Decay::Total()
{
	return nnn() + nGAMMA() + nEE() + 2.0*nEMU() + nPI0() +
	       2.0*EPI() + nMUMU() + 2.0*MUPI() + 2.0*EKA() + nKA0();
}


//individual decay channels
//
double Decay::nnn()
{
	if (M_Sterile >= 3.0 * M_Neutrino)
		return mapEnhance["nnn"] * genie::constants::kGF2 * pow(M_Sterile, 5) * 
			(U_e*U_e + U_m*U_m + U_t*U_t) / (96.0 * genie::constants::kPi3);
	else return 0.0;
}

double Decay::nGAMMA()
{
	double AemPi = genie::constants::kAem / genie::constants::kPi;

	if (M_Sterile >= 3.0 * M_Neutrino)
		return mapEnhance["nGAMMA"] * genie::constants::kGF2 * pow(M_Sterile, 5) *
			(U_e*U_e + U_m*U_m + U_t*U_t) * (27.0/32.0 * AemPi) /
			(192.0 * genie::constants::kPi3);
	else return 0.0;
}

//M_Sterile > 2 M_Electron (always)
double Decay::nEE()
{
	if (M_Sterile >= 2 * M_Electron)
	{
		double dMe = M_Electron / M_Sterile;
		double dMn = M_Neutrino / M_Sterile;
		double gL = -0.5 + Tools::Const::fSin2w;
		double gR = Tools::Const::fSin2w;
		double Int1 = I1_xyz(dMn, dMe, dMe);
		double Int2 = I2_xyz(dMn, dMe, dMe);
		double KF_e = (gL*gR + gR) * Int2 + (gL*gL + gR*gR + (1+2*gL))*Int1;
		double KF_m = (gL*gR) * Int2 + (gL*gL + gR*gR)*Int1;
//		double KF_t = (gL*gR) * I2_xyz(dMn, dMe, dMe) + (gL*gL + gR*gR)*I1_xyz(dMn, dMe, dMe);

		return mapEnhance["nEE"] * genie::constants::kGF2 * pow(M_Sterile, 5) * 
			(U_e*U_e * KF_e + (U_m*U_m + U_t*U_t) * KF_m) / 
			(96.0 * genie::constants::kPi3);
	}
	else return 0.0;
}

//M_Sterile > M_Muon
double Decay::nEMU()	//Valid for electron+antimuon and positron+muon
{
	if (M_Sterile >= M_Electron + M_Muon)
	{
		double dMe = M_Electron / M_Sterile;
		double dMm = M_Muon / M_Sterile;
		double dMn = M_Neutrino / M_Sterile;

		return mapEnhance["nEMU"] * genie::constants::kGF2 * pow(M_Sterile, 5) * 
			(U_e*U_e * I1_xyz(dMm, dMn, dMe) + U_m*U_mi * I1_xyz(dMe, dMn, dMm)) / 
			(192.0 * genie::constants::kPi3);
	}
	else return 0.0;
}

//M_Sterile > M_Pion0
double Decay::nPI0()
{
	if (M_Sterile >= M_Pion0)
	{
		double dMp2 = M_Pion0*M_Pion0/M_Sterile/M_Sterile;

		return mapEnhance["nPI0"] * genie::constants::kGF2 * pow(M_Sterile, 3) *
			(U_e*U_e + U_m*U_m + U_t*U_t) * 
			pow((1.0-dMp2), 2.0) * Tools::Const::fFPion2 / 
			(64.0 * genie::constants::kPi);
	}
	else return 0.0;
}

//M_Sterile > M_Pion
double Decay::EPI()
{
	if (M_Sterile >= M_Electron + M_Pion)
	{
		double dMe2 = M_Electron*M_Electron/M_Sterile/M_Sterile;
		double dMp2 = M_Pion*M_Pion/M_Sterile/M_Sterile;

		return mapEnhance["EPI"] * genie::constants::kGF2 * pow(M_Sterile, 3) *
		       U_e*U_e * 
		       pow(Tools::Const::fV_ud, 2.0) * Tools::Const:fFPion2 * I1_xy(dMe2, dMp2) / 
		       (16.0 * genie::constants::kPi);
	}
	else return 0.0;
}

//M_Sterile > 2 M_Muon
double Decay::nMUMU()
{
	if (M_Sterile >= 2 * M_Muon)
	{
		double dMm = M_Muon / M_Sterile;
		double dMn = M_Neutrino / M_Sterile;
		double gL = -0.5 + Tools::Const::fSin2w;
		double gR = Tools::Const::fSin2w;
		double Int1 = I1_xyz(dMn, dMm, dMm);
		double Int2 = I2_xyz(dMn, dMm, dMm);
		double KF_m = (gL*gR + gR) * Int2 + (gL*gL + gR*gR + (1+2*gL))*Int1;
		double KF_e = (gL*gR) * Int2 + (gL*gL + gR*gR)*Int1;
//		double KF_t = (gL*gR) * I2_xyz(dMn, dMe, dMe) + (gL*gL + gR*gR)*I1_xyz(dMn, dMe, dMe);

		return mapEnhance["nMUMU"] * genie::constants::kGF2 * pow(M_Sterile, 5) *
		       (U_m*U_m * KF_m + (U_e*U_e + U_t*U_t) * KF_e) /
		       (96.0 * genie::constants::kPi3);
	}
	else return 0.0;
}

//M_Sterile > M_Pion + M_Muon
double Decay::MUPI()
{
	if (M_Sterile >= M_Muon + M_Pion)
	{
		double dMm2 = M_Muon*M_Muon/M_Sterile/M_Sterile;
		double dMp2 = M_Pion*M_Pion/M_Sterile/M_Sterile;

		return mapEnhance["MUPI"] * genie::constants::kGF2 * pow(M_Sterile, 3) *
		       U_m*U_m * 
		       pow(Tools::Const::fV_ud, 2.0)*Tools::Const::fFPion2 * I1_xy(dMm2, dMp2) /
		       / (16.0 * genie::constants::kPi);
	}
	else return 0.0;
}

//M_Sterile > M_Kaon + M_Electron
double Decay::EKA()
{
	if (M_Sterile >= M_Electron + M_Kaon)
	{
		double dMe2 = M_Electron*M_Electron/M_Sterile/M_Sterile;
		double dMk2 = M_Kaon*M_Kaon/M_Sterile/M_Sterile;

		return mapEnhance["EKA"] * genie::constants::kGF2 * pow(M_Sterile, 3) *
		       U_e*U_e * 
		       pow(Tools::Const::fV_us, 2.0) * Tools::Const::fFKaon2 * I1_xy(dMm2, dMp2) /
		       (16.0 * genie::constants::kPi);
	}
	else return 0.0;
}

//M_Sterile > M_Kaon0
double Decay::nKA0()
{
	if (M_Sterile >= M_Kaon0)
	{
		double dMp2 = M_Kaon0*M_Kaon0/M_Sterile/M_Sterile;

		return mapEnhance["nKA0"] * genie::constants::kGF2 * pow(M_Sterile, 3) *
			(U_e*U_e + U_m*U_m + U_t*U_t)
			Tools::Const::fFKaon2 * pow((1.0-dMp2), 2.0) / 
			(64.0 * genie::constants::kPi);
	}
	else return 0.0;
}


//Get functions
double Decay::GetMSterile()
{
	return M_Sterile;
}

double Decay::GetUe()
{
	return U_e;
}

double Decay::GetUm()
{
	return U_m;
}

double Decay::GetUt()
{
	return U_t;
}

//Set functions
void Decay::SetMSterile(double X)
{
	M_Sterile = X;
}

void Decay::SetUe(double X)
{
	U_e = X;
}

void Decay::SetUm(double X)
{
	U_m = X;
}

void Decay::SetUt(double X)
{
	U_t = X;
}
