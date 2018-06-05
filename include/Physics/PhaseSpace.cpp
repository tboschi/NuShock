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
			Ret = nPI0()/Max_nPI0();
			break;
		case _EPI:
			Ret = EPI()/Max_EPI();
			break;
		case _MUPI:
			Ret = MUPI()/Max_MUPI();
			break;
		case _TPI:
			Ret = TPI()/Max_TPI();
			break;
		case _EKA:
			Ret = EKA()/Max_EKA();
			break;
		case _MUKA:
			Ret = MUKA()/Max_MUKA();
			break;
		case _nRHO0:
			Ret = nRHO0()/Max_nRHO0();
			break;
		case _ERHO:
			Ret = ERHO()/Max_ERHO();
			break;
		case _MURHO:
			Ret = MURHO()/Max_MURHO();
			break;
		case _EKAx:
			Ret = EKAx()/Max_EKAx();
			break;
		case _MUKAx:
			Ret = MUKAx()/Max_MUKAx();
			break;
		case _nOMEGA:
			Ret = nOMEGA()/Max_nOMEGA();
			break;
		case _nETA:
			Ret = nETA()/Max_nETA();
			break;
		case _nETAi:
			Ret = nETAi()/Max_nETAi();
			break;
		case _nPHI:
			Ret = nPHI()/Max_nPHI();
			break;
		case _ECHARM:
			Ret = ECHARM()/Max_ECHARM();
			break;
		case _MuonE:
			Ret = MuonE()/Max_MuonE();
			break;
		case _MuonM:
			Ret = MuonM()/Max_MuonM();
			break;
		case _TauEE:
			Ret = TauEE()/Max_TauEE();
			break;
		case _TauET:
			Ret = TauET()/Max_TauET();
			break;
		case _TauME:
			Ret = TauME()/Max_TauME();
			break;
		case _TauMT:
			Ret = TauMT()/Max_TauMT();
			break;
		case _PionE:
			Ret = PionE()/Max_PionE();
			break;
		case _PionM:
			Ret = PionM()/Max_PionM();
			break;
		case _KaonE:
			Ret = KaonE()/Max_KaonE();
			break;
		case _KaonM:
			Ret = KaonM()/Max_KaonM();
			break;
		case _CharmE:
			Ret = CharmE()/Max_CharmE();
			break;
		case _CharmM:
			Ret = CharmN()/Max_CharmN();
			break;
		case _CharmT:
			Ret = CharmT()/Max_CharmT();
			break;
		case _Kaon0E:
			Ret = Kaon0E()/Max_Kaon0E();
			break;
		case _Kaon0M:
			Ret = Kaon0M()/Max_Kaon0M();
			break;
		case _KaonCE:
			Ret = KaonCE()/Max_KaonCE();
			break;
		case _KaonCM:
			Ret = KaonCM()/Max_KaonCM();
			break;
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

	double gL = -0.5 + Const::fSin2W;
	double gR = Const::fSin2W;

	double MaxWW = Max_WW(dMN2, dML2, dML2);	//W or Z mediation
	double MaxWZ = Max_WZ(dMN2, dML2, dML2);	//cross term for W + Z

	double Amp_AAA = (gL*gL + gR*gR + 1 + 2*gL) * MaxWW + (gL*gR + gR) * MaxWW;	//both W and Z
	double Amp_ABB = (gL*gL + gR*gR           ) * MaxWW + (gL*gR     ) * MaxWZ;	//only Z

	fAAA = Const::fGF2 * pow(GetMass(), 5) * Amp_AAA / (16.0 * Const::fPi3);	//nu flavour is the same of leptons
	fABB = Const::fGF2 * pow(GetMass(), 5) * Amp_AAA / (16.0 * Const::fPi3);	//nu flavour is different from leptons

	return 0.0;
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

double PhaseSpace::Max_WW(double x, double y, double z)	//depends on s, t, cos0, cos1
{
	I_var.clear();

	I_var.push_back(x);	//0
	I_var.push_back(y);	//1
	I_var.push_back(z);	//2

	SetFunction_D(&Max_WW_D);
	std::vector<double? xPos;
	return Inte::NeldMedSolver(this, xPos, 4, Inte::Maximum);
}

double PhaseSpace::Max_WW_D(double *x)
{
	const double &x = I_var.at(0);
	const double &y = I_var.at(1);
	const double &z = I_var.at(2);
	double &s = x[0];
	double &t = x[1];
	double &cos0 = x[2];
	double &cos1 = x[3];

	//define vars
	//double Diff = stLimit(double &s, double &t);
	//double cos0 = -1 + 2*cos0;

	return M2_WW(x, y, z, s, t, cos0, cos1);		//defined in DecayRates
}

double DecayRates::Max_WW(double x, double y, double z)	//depends on s, t, cos0, cos1
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

	SetFunction(&I_WW_t);
	double Ret = Inte::Max(this);		//max in cos0 for given s

	I_var.pop_back();
	SetFunction(&I_WW_s);
	return Ret;
}

double PhaseSpace::Max_WW_t(double t)
{
	I_var.push_back(t);	//3

	SetFunction(&I_WW_cos0);
	double Ret = Inte::Max(this);		//max in cos0 for given s

	I_var.pop_back();
	SetFunction(&I_WW_t);
	return Ret;
}

double PhaseSpace::Max_WW_t(double t)
{
	I_var.push_back(t);	//3

	SetFunction(&I_WW_cos0);
	double Ret = Inte::Max(this);		//max in cos0 for given s

	I_var.pop_back();
	SetFunction(&I_WW_t);
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

double DecayRates::Max_WZ(double x, double y, double z)	//depends on s, cos0
{
	I_var.clear();

	I_var.push_back(x);	//0
	I_var.push_back(y);	//1
	I_var.push_back(z);	//2

	SetFunction(&I_WZ_s);

	return Inte::Max(this);
}

double PhaseSpace::Max_WZ_s(double s)
{
	I_var.push_back(s);	//3

	SetFunction(&I_WZ_cos0);
	double Ret = Inte::Max(this);		//max in cos0 for given s

	I_var.pop_back();
	SetFunction(&I_WZ_s);
	return Ret;
}

double PhaseSpace::Max_WZ_cos0(double cos0)
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


//////////
//	//
//	//
//////////

void PhaseSpace::Boost()
{
	for (unsigned int i = 0; i < nDaughter; ++i)
		TheSpace->GetDecayProduct
	N_vec->Boost();
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
