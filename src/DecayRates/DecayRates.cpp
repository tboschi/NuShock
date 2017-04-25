#include "DecayRates.h"

Decay::Decay(double MSterile, double Ue, double Um, double Ut)	: //Decay rates calculator
	M_Neutrino(0.0),
	M_Photon(0.0),
	M_Electron(Const::fMElectron),
	M_Muon(Const::fMMuon),
	M_Pion(Const::fMPion),
	M_Pion0(Const::fMPion0),
	M_Kaon(Const::fMKaon),
	M_Kaon0(Const::fMKaon0)
{
	SetMSterile(MSterile);
	SetUe(Ue);
	SetUm(Um);
	SetUt(Ut);

	MapInit();
	SetEnhancement();
	Event = new TGenPhaseSpace;
	N_vec = new TLorentzVector;
}

//Initialisation of map
void Decay::MapInit()
{
	mapChannel["ALL"] = _ALL;
	mapChannel["nnn"] = _nnn;
	mapChannel["nGAMMA"] = _nGAMMA;
	mapChannel["nEE"] = _nEE;
	mapChannel["nEMU"] = _nEMU;
	mapChannel["nPI0"] = _nPI0;
	mapChannel["EPI"] = _EPI;
	mapChannel["nMUMU"] = _nMUMU;
	mapChannel["MUPI"] = _MUPI;
	mapChannel["EKA"] = _EKA;
	mapChannel["nKA0"] = _nKA0;
}

double Decay::Gamma(std::string Channel, double B)
{
	SetEnhancement(Channel, B);

	double Result;
	switch(mapChannel[Channel])
	{
		case _ALL:
			Result = Total();
			break;
		case _nnn:
			Result = nnn();
			break;
		case _nGAMMA:
			Result = nGAMMA();
			break;
		case _nEE:
			Result = nEE();
			break;
		case _nEMU:
			Result = nEMU();
			break;
		case _nPI0:
			Result = nPI0();
			break;
		case _EPI:
			Result = EPI();
			break;
		case _nMUMU:
			Result = nMUMU();
			break;
		case _MUPI:
			Result = MUPI();
			break;
		case _EKA:
			Result = EKA();
			break;
		case _nKA0:
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
		case _ALL:
			Result = 0.0;
			break;
		case _nnn:
			Result = Total()-nnn();
			break;
		case _nGAMMA:
			Result = Total()-nGAMMA();
			break;
		case _nEE:
			Result = Total()-nEE();
			break;
		case _nEMU:
			Result = Total()-nEMU();
			break;
		case _nPI0:
			Result = Total()-nPI0();
			break;
		case _EPI:
			Result = Total()-EPI();
			break;
		case _nMUMU:
			Result = Total()-nMUMU();
			break;
		case _MUPI:
			Result = Total()-MUPI();
			break;
		case _EKA:
			Result = Total()-EKA();
			break;
		case _nKA0:
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

	double Result;
	switch(mapChannel[Channel])
	{
		case _ALL:
			Result = 1.0;
			break;
		case _nnn:
			Result = nnn()/Total();
			break;
		case _nGAMMA:
			Result = nGAMMA()/Total();
			break;
		case _nEE:
			Result = nEE()/Total();
			break;
		case _nEMU:
			Result = nEMU()/Total();
			break;
		case _nPI0:
			Result = nPI0()/Total();
			break;
		case _EPI:
			Result = EPI()/Total();
			break;
		case _nMUMU:
			Result = nMUMU()/Total();
			break;
		case _MUPI:
			Result = MUPI()/Total();
			break;
		case _EKA:
			Result = EKA()/Total();
			break;
		case _nKA0:
			Result = nKA0()/Total();
			break;
		default:
			Result = -1.0;
			break;
	}

	SetEnhancement();
	return Result;
}

int Decay::PhaseSpace(std::string Channel, double &Weight)	//Return number of products 
{								//0 if decay not valid
	SetEnhancement();

	double Mass[3];
	int Products;

	switch(mapChannel[Channel])
	{
		case _nnn:
			Mass[0] = M_Neutrino;
			Mass[1] = M_Neutrino;
			Mass[2] = M_Neutrino;
			Products = 3 * Event->SetDecay(*N_vec, 3, Mass);
			Weight = Event->Generate();
			break;
		case _nGAMMA:
			Mass[0] = M_Neutrino;
			Mass[1] = M_Photon;
			Products = 2 * Event->SetDecay(*N_vec, 3, Mass);
			Weight = Event->Generate();
			break;
		case _nEE:
			Mass[0] = M_Neutrino;
			Mass[1] = M_Electron;
			Mass[2] = M_Electron;
			Products = 3 * Event->SetDecay(*N_vec, 3, Mass);
			Weight = Event->Generate();
			break;
		case _nEMU:
			Mass[0] = M_Neutrino;
			Mass[1] = M_Electron;
			Mass[2] = M_Muon;
			Products = 3 * Event->SetDecay(*N_vec, 3, Mass);
			Weight = Event->Generate();
			break;
		case _nPI0:
			Mass[0] = M_Neutrino;
			Mass[1] = M_Pion0;
			Products = 2 * Event->SetDecay(*N_vec, 2, Mass);
			Weight = Event->Generate();
			break;
		case _EPI:
			Mass[0] = M_Electron;
			Mass[1] = M_Pion;
			Products = 2 * Event->SetDecay(*N_vec, 2, Mass);
			Weight = Event->Generate();
			break;
		case _nMUMU:
			Mass[0] = M_Neutrino;
			Mass[1] = M_Muon;
			Mass[2] = M_Muon;
			Products = 3 * Event->SetDecay(*N_vec, 3, Mass);
			Weight = Event->Generate();
			break;
		case _MUPI:
			Mass[0] = M_Muon;
			Mass[1] = M_Pion;
			Products = 2 * Event->SetDecay(*N_vec, 2, Mass);
			Weight = Event->Generate();
			break;
		case _EKA:
			Mass[0] = M_Electron;
			Mass[1] = M_Kaon;
			Products = 2 * Event->SetDecay(*N_vec, 2, Mass);
			Weight = Event->Generate();
			break;
		case _nKA0:
			Mass[0] = M_Neutrino;
			Mass[1] = M_Kaon0;
			Products = 2 * Event->SetDecay(*N_vec, 2, Mass);
			Weight = Event->Generate();
			break;
		default:
			Products = 0;
			break;
	}

	return Products;
}

//Controller of decay enhancement
void Decay::SetEnhancement(std::string Channel, double K)
{
	std::map<std::string, ChannelName>::iterator it;
	if (Channel == "ALL")
		for (it = mapChannel.begin(); it != mapChannel.end(); ++it)
			mapEnhance[it->first] = K;
	else
	{
		it = mapChannel.find(Channel);
		if (it != mapChannel.end())
			mapEnhance[it->first] = K;
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
	{
		return mapEnhance["nnn"] * genie::constants::kGF2 * pow(M_Sterile, 5) * 
			(U_e*U_e + U_m*U_m + U_t*U_t) / (96.0 * genie::constants::kPi3);
	}
	else return 0.0;
}

double Decay::nGAMMA()
{
	double AemPi = genie::constants::kAem / genie::constants::kPi;

	if (M_Sterile >= M_Neutrino + M_Photon)
	{
		return mapEnhance["nGAMMA"] * genie::constants::kGF2 * pow(M_Sterile, 5) *
			(U_e*U_e + U_m*U_m + U_t*U_t) * (27.0/32.0 * AemPi) /
			(192.0 * genie::constants::kPi3);
	}
	else return 0.0;
}

//M_Sterile > 2 M_Electron (always)
double Decay::nEE()
{
	if (M_Sterile >= 2 * M_Electron)
	{
		double dMe = M_Electron / M_Sterile;
		double dMn = M_Neutrino / M_Sterile;
		double gL = -0.5 + Const::fSin2W;
		double gR = Const::fSin2W;
		double Int1 = Kine::I1_xyz(dMn, dMe, dMe);
		double Int2 = Kine::I2_xyz(dMn, dMe, dMe);
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
			(U_e*U_e * Kine::I1_xyz(dMm, dMn, dMe) +
			U_m*U_m * Kine::I1_xyz(dMe, dMn, dMm)) / 
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
			pow((1.0-dMp2), 2.0) * Const::fFPion2 / 
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
		       pow(Const::fV_ud, 2.0) * Const::fFPion2 * Kine::I1_xy(dMe2, dMp2) / 
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
		double gL = -0.5 + Const::fSin2W;
		double gR = Const::fSin2W;
		double Int1 = Kine::I1_xyz(dMn, dMm, dMm);
		double Int2 = Kine::I2_xyz(dMn, dMm, dMm);
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
		       pow(Const::fV_ud, 2.0)*Const::fFPion2 * Kine::I1_xy(dMm2, dMp2) /
		       (16.0 * genie::constants::kPi);
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
		       pow(Const::fV_us, 2.0) * Const::fFKaon2 * Kine::I1_xy(dMe2, dMk2) /
		       (16.0 * genie::constants::kPi);
	}
	else return 0.0;
}

//M_Sterile > M_Kaon0
double Decay::nKA0()
{
	if (M_Sterile >= M_Kaon0)
	{
		double dMk2 = M_Kaon0*M_Kaon0/M_Sterile/M_Sterile;

		return mapEnhance["nKA0"] * genie::constants::kGF2 * pow(M_Sterile, 3) *
			(U_e*U_e + U_m*U_m + U_t*U_t) * 
			Const::fFKaon2 * pow((1.0-dMk2), 2.0) / 
			(64.0 * genie::constants::kPi);
	}
	else return 0.0;
}


std::vector<std::string> Decay::ListChannels()
{
	std::map<std::string, ChannelName>::iterator it;
	std::vector<std::string> vList;
	for (it = mapChannel.begin(); it != mapChannel.end(); ++it)
		vList.push_back(it->first);
	return vList;
}

//Get functions
TLorentzVector *Decay::GetNvec()
{
	return N_vec;
}

TLorentzVector *Decay::GetDecayProduct(int i)
{
	return Event->GetDecay(i);
}

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
void Decay::SetNvec(TLorentzVector &X)
{
	*N_vec = X;
}

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
