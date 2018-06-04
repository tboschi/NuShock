#include "PhaseSpace.h"

PhaseSpace::PhaseSpace(Neutrino *N) :
	Event(new TGenPhaseSpace),
	GenMT(new TRandom3(0))
{
	SetNeutrino(N);

	TheSpace = new TGenPhaseSpace;
	N_labf = new TLorentzVector;
	N_rest = new TLorentzVector;
}

bool PhaseSpace::SetDecay(Channel Name, TLorentzVector &Rest)
{
	if (Channel_prev != Name)
	{
		LoadMass(Name);
		Channel_prev = Name;
	}

	double MassArray[10];
	nDaughter = vMass.size();
	for (unsigned int i = 0; i < nDaughter; ++i)
		MassArray[i] = Masses.at(i);

	return TheSpace->SetDecay(Rest, nDaughter, MassArray);
}

bool PhaseSpace::Generate(Channel Name)
{
	if (SetDecay(Name))
	{
		double Max = 1.0;
		do
			Event->Generate();
		while (GenMT->Ratio(Name));
	}
	else
		return false;

	return true;
}

double PhaseSpace::Ratio(Channel Name)
{
	double Ret = -1.0;
	switch (Name)
	{
		case _ALL:
		case _nnn:
		case _nGAMMA:
			Ret = nGAMMA()/Max_nGAMMA();
			break;
		case _nEE:
			Ret = nEE()/Max_nEE();
			break;
		case _nEMU:
			Ret = nEMU()/Max_nEMU();
			break;
		case _nMUE:
			Ret = nMUE()/Max_nMUE();
			break;
		case _nMUMU:
			Ret = nMUMU()/Max_nMUMU();
			break;
		case _nET:
			Ret = nET()/Max_nET();
			break;
		case _nTE:
			Ret = nTE()/Max_nTE();
			break;
		case _nMUT:
			Ret = nMUT()/Max_nMUT();
			break;
		case _nTMU:
			Ret = TMU()/Max_TMU();
			break;
		case _nPI0:
			break;
			Ret = nPI0()/Max_nPI0();
		case _EPI:
			break;
			Ret = EPI()/Max_EPI();
		case _MUPI:
			break;
			Ret = MUPI()/Max_MUPI();
		case _TPI:
			break;
			Ret = TPI()/Max_TPI();
		case _EKA:
			break;
			Ret = EKA()/Max_EKA();
		case _MUKA:
			break;
			Ret = MUKA()/Max_MUKA();
		case _nRHO0:
			break;
			Ret = nRHO0()/Max_nRHO0();
		case _ERHO:
			break;
			Ret = ERHO()/Max_ERHO();
		case _MURHO:
			break;
			Ret = MURHO()/Max_MURHO();
		case _EKAx:
			break;
			Ret = EKAx()/Max_EKAx();
		case _MUKAx:
			break;
			Ret = MUKAx()/Max_MUKAx();
		case _nOMEGA:
			break;
			Ret = nOMEGA()/Max_nOMEGA();
		case _nETA:
			break;
			Ret = nETA()/Max_nETA();
		case _nETAi:
			break;
			Ret = nETAi()/Max_nETAi();
		case _nPHI:
			break;
			Ret = nPHI()/Max_nPHI();
		case _ECHARM:
			break;
			Ret = ECHARM()/Max_ECHARM();
		case _MuonE:
			break;
			Ret = MuonE()/Max_MuonE();
		case _MuonM:
			break;
			Ret = MuonM()/Max_MuonM();
		case _TauEE:
			break;
			Ret = TauEE()/Max_TauEE();
		case _TauET:
			break;
			Ret = TauET()/Max_TauET();
		case _TauME:
			break;
			Ret = TauME()/Max_TauME();
		case _TauMT:
			break;
			Ret = TauMT()/Max_TauMT();
		case _PionE:
			break;
			Ret = PionE()/Max_PionE();
		case _PionM:
			break;
			Ret = PionM()/Max_PionM();
		case _KaonE:
			break;
			Ret = KaonE()/Max_KaonE();
		case _KaonM:
			break;
			Ret = KaonM()/Max_KaonM();
		case _CharmE:
			break;
			Ret = CharmE()/Max_CharmE();
		case _CharmM:
			break;
			Ret = CharmN()/Max_CharmN();
		case _CharmT:
			break;
			Ret = CharmT()/Max_CharmT();
		case _Kaon0E:
			break;
			Ret = Kaon0E()/Max_Kaon0E();
		case _Kaon0M:
			break;
			Ret = Kaon0M()/Max_Kaon0M();
		case _KaonCE:
			break;
			Ret = KaonCE()/Max_KaonCE();
		case _KaonCM:
			break;
			Ret = KaonCM()/Max_KaonCM();
		default:
	}
}

double PhaseSpace::nEE()
{
	TLorentzVector vec_n0 = *(Event->GetDecayProduct(0));
	TLorentzVector vec_E1 = *(Event->GetDecayProduct(1));
	TLorentzVector vec_E2 = *(Event->GetDecayProduct(2));

	TLorentzVector vec_ss = N_rest - vec_E1;
	double cos0 = vec_E1->CosTheta();
	double s = vec_ss->M2();

	double fnEE_e, fmEE_mt;
	M2_NeutrinoLeptonAA(fnEE_e, fmEE_mt, M_Neutrino, M_Electron, s, cos0);		//defined in DecayRates
	return fnEE_e * Ue*Ue + fnEE_mt * (Um*Um + Ut*Ut);
}

double PhaseSpace::nEMU()
{
	TLorentzVector vec_n0  = *(Event->GetDecayProduct(0));
	TLorentzVector vec_E1  = *(Event->GetDecayProduct(1));
	TLorentzVector vec_MU2 = *(Event->GetDecayProduct(2));

	TLorentzVector vec_ss = N_rest - vec_E1;
	double cos0 = vec_E1->CosTheta();
	double s = vec_ss->M2();

	return NeutrinoLeptonAB(fAAA, fABB, M_Neutrino, M_Electron, s, cos0);		//defined in DecayRates
}

double PhaseSpace::nEMU()
{
	TLorentzVector vec_n0  = *(Event->GetDecayProduct(0));
	TLorentzVector vec_E1  = *(Event->GetDecayProduct(1));
	TLorentzVector vec_MU2 = *(Event->GetDecayProduct(2));

	TLorentzVector vec_ss = N_rest - vec_E1;
	double cos0 = vec_E1->CosTheta();
	double s = vec_ss->M2();

	return NeutrinoLeptonAB(fAAA, fABB, M_Neutrino, M_Electron, s, cos0);		//defined in DecayRates
}

double PhaseSpace::nMUE()
{
	TLorentzVector vec_n0  = *(Event->GetDecayProduct(0));
	TLorentzVector vec_E1  = *(Event->GetDecayProduct(1));
	TLorentzVector vec_MU2 = *(Event->GetDecayProduct(2));

	TLorentzVector vec_ss = N_rest - vec_E1;
	double cos0 = vec_E1->CosTheta();
	double s = vec_ss->M2();

	double fAAA = NeutrinoLeptonAB(fAAA, fABB, M_Neutrino, M_Electron, s, cos0);		//defined in DecayRates
}

double PhaseSpace::nMUMU()
{
	TLorentzVector vec_n0  = *(Event->GetDecayProduct(0));
	TLorentzVector vec_E1  = *(Event->GetDecayProduct(1));
	TLorentzVector vec_MU2 = *(Event->GetDecayProduct(2));

	TLorentzVector vec_ss = N_rest - vec_E1;
	double cos0 = vec_E1->CosTheta();
	double s = vec_ss->M2();

	double fAAA = NeutrinoLeptonAB(fAAA, fABB, M_Neutrino, M_Electron, s, cos0);		//defined in DecayRates
}

double PhaseSpace::nET()
{
	TLorentzVector vec_n0  = *(Event->GetDecayProduct(0));
	TLorentzVector vec_E1  = *(Event->GetDecayProduct(1));
	TLorentzVector vec_MU2 = *(Event->GetDecayProduct(2));

	TLorentzVector vec_ss = N_rest - vec_E1;
	double cos0 = vec_E1->CosTheta();
	double s = vec_ss->M2();

	double fAAA = NeutrinoLeptonAB(fAAA, fABB, M_Neutrino, M_Electron, s, cos0);		//defined in DecayRates
}

double PhaseSpace::nTE()
{
	TLorentzVector vec_n0  = *(Event->GetDecayProduct(0));
	TLorentzVector vec_E1  = *(Event->GetDecayProduct(1));
	TLorentzVector vec_MU2 = *(Event->GetDecayProduct(2));

	TLorentzVector vec_ss = N_rest - vec_E1;
	double cos0 = vec_E1->CosTheta();
	double s = vec_ss->M2();

	double fAAA = NeutrinoLeptonAB(fAAA, fABB, M_Neutrino, M_Electron, s, cos0);		//defined in DecayRates
}

double PhaseSpace::Max_nEE()
{
	if (fmax_nEE < 0 || IsChanged())
		fmax_nEE = Max_NeutrinoLeptonAA(M_Electron, M_Electron);

	return fmax_nEE;
}

double PhaseSpace::Max_nMUE()
{
	if (fmax_nEE < 0 || IsChanged())
		fmax_nEE = Max_NeutrinoLeptonAB(M_Electron, M_Muon);

	return fmax_nEE;
}

//return the max
double PhaseSpace::NeutrinoLeptonAA(double &fAAA, double &fABB, double M_Neutrino, double M_Lepton, double s, double theta)
{
	double dMN2 = M_Neutrino*M_Neutrino/GetMass()/GetMass();
	double dML2 = M_Lepton*M_Lepton/GetMass()/GetMass();

	Generate;
	GetDecayProduct;
	M2_NeutrinoLeptonAA();
	while (GenMT->Rndm() < M2_NeutrinoLeptonAA / Max);
}

double PhaseSpace::Max_NeutrinoLeptonAB(double M_LeptonA, double M_LeptonB, double M_Neutrino, double s, double theta)
{
	double dMN2 = M_Neutrino*M_Neutrino/GetMass()/GetMass();
	double dML2 = M_Lepton*M_Lepton/GetMass()/GetMass();

	double dMA2 = M_LeptonA*M_LeptonB/GetMass()/GetMass();
	double dMB2 = M_LeptonB*M_LeptonB/GetMass()/GetMass();

	return Const::fGF2 / (16.0 * Const::fPi3) * 
	       pow(GetMass(), 5) * Max_WW(dMN2, dMA2, dMB2, theta);
}

double DecayRates::Max_WW(double x, double y, double z)
{
	I_var.clear();

	I_var.push_back(x);	//0
	I_var.push_back(y);	//1
	I_var.push_back(z);	//2

	SetFunction(&I_WW_s);

	return Inte::Max(this);
}

double PhaseSpace::Max_WW_s(double s)
{
	I_var.push_back(s);	//3

	SetFunction(&I_WW_cos0);
	double Ret = Inte::Max(this);		//max in cos0 for given s

	I_var.pop_back();
	SetFunction(&I_WW_s);
	return Ret;
}

double PhaseSpace::Max_WW_cos0(double cos0)
{
	const double &x = I_var.at(0);
	const double &y = I_var.at(1);
	const double &z = I_var.at(2);
	double s = I_var.at(3);
	double t = 0;

	//define vars
	double Diff = stLimit(double &s, double &t);
	cos0 = -1 + 2*cos0;

	return M2_WW(x, y, z, s, cos0);		//defined in DecayRates
}

double PhaseSpace::I_WW(double x, double y, double z, double theta)
{
	I_var.clear();

	I_var.push_back(x);	//0
	I_var.push_back(y);	//1
	I_var.push_back(z);	//2
	I_var.push_back(cos0);	//3

	SetFunction(&I_WW_s);				//Function will fix the integration volume
	return Inte::BooleIntegration(this); 
}

double PhaseSpace::I_WW_s(double s)
{
	const double &x = I_var.at(0);
	const double &y = I_var.at(1);
	const double &z = I_var.at(2);
	const double &cos0 = I_var.at(3);

	return I_WW_s(s, cos0, x, y, z);
}

double PhaseSpace::LoadPhaseSpace(double Max)
{
	double / dMax;
	return TheSpace->Generate();
}

void PhaseSpace::Boost()
{
	for (unsigned int i = 0; i < nDaughter; ++i)
		TheSpace->GetDecayProduct
	N_vec->Boost();
}

double PhaseSpace::MaxThreeGamma()
{
	SetFunction(&I_WZ_s);
	Inte::MaxMin(this, Max, Min);
}


double PhaseSpace::NeutrinoLeptonAB(double M_LeptonA, double M_LeptonB, double M_Neutrino)
{
	SetFunction(dNeutrinoLeptonAB);
}

double PhaseSpace::dNeutrinoLeptonAB(double M_LeptonA, double M_LeptonB, double M_Neutrino, double s, double cos0)
{
	double dMN2 = M_Neutrino*M_Neutrino/GetMass()/GetMass();
	double dMA2 = M_LeptonA*M_LeptonB/GetMass()/GetMass();
	double dMB2 = M_LeptonB*M_LeptonB/GetMass()/GetMass();

	return GetMult() * Const::fGF2 / (16.0 * Const::fPi3) * 
		pow(GetMass(), 5) * I_WW_s(dMN2, dMA2, dMB2);
}

double PhaseSpace::MaxNeutrinoLeptonAB(double M_LeptonA, double M_LeptonB, double M_Neutrino)
{
	double dMN2 = M_Neutrino*M_Neutrino/GetMass()/GetMass();
	double dMA2 = M_LeptonA*M_LeptonB/GetMass()/GetMass();
	double dMB2 = M_LeptonB*M_LeptonB/GetMass()/GetMass();

	return GetMult() * Const::fGF2 / (16.0 * Const::fPi3) * 
		pow(GetMass(), 5) * Max_WW(dMN2, dMA2, dMB2);
}

TLorentzVector *PhaseSpace::GetNvec()
{
	return N_vec;
}

TLorentzVector PhaseSpace::GetDecayProduct(int i, int &ID)
{
	TLorentzVector Daughter = *(Event->GetDecay(i));
	Daughter.Boost(N_labf->BoostVector());

	return Daughter;
}

void PhaseSpace::SetNLabf(TLorentzVector &Vec)
{
	*N_labf = Vec;
	SetNRest(N_labf->M());
}

void PhaseSpace::SetNRest(double Mass)
{
	N_rest->SetPxPyPzE(0, 0, 0, Mass);
}

