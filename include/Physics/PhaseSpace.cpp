#include "PhaseSpace.h"

PhaseSpace::PhaseSpace() :
	Event(new TGenPhaseSpace),
	GenMT(new TRandom3(0)),
	N_labf(new TLorentzVector()),
	N_rest(new TLorentzVector())
{
}

bool PhaseSpace::SetDecay(Channel Name)
{
	if (Channel_prev != Name)
	{
		LoadMass(Name);
		Channel_prev = Name;
	}

	double MassArray[10];
	for (unsigned int i = 0; i < Daughter(); ++i)
		MassArray[i] = vMass.at(i);

	return Event->SetDecay(*Rest(), Daughter(), MassArray);
}

bool PhaseSpace::Generate(Channel Name)
{
	if (SetDecay(Name))
	{
		do
			Event->Generate();
		while (GenMT->Rndm() > Ratio(Name));

		return true;
	}
	else
		return false;
}

double PhaseSpace::Ratio(Channel Name)
{
	double Ret = 1.0;
	switch (Name)
	{
		case _ALL:
		case _nnn:
		case _nGAMMA:
			break;
		case _nEE:
			Ret = nEE_ratio();
			break;
		case _nEM:
			Ret = nEM_ratio();
			break;
		case _nME:
			Ret = nME_ratio();
			break;
		case _nMM:
			Ret = nMM_ratio();
			break;
		case _nET:
			Ret = nET_ratio();
			break;
		case _nTE:
			Ret = nTE_ratio();
			break;
		case _nMT:
			Ret = nMT_ratio();
			break;
		case _nTM:
			Ret = nTM_ratio();
			break;
		case _nPI0:
			Ret = nPI0_ratio();
			break;
		case _EPI:
			Ret = EPI_ratio();
			break;
		case _MPI:
			Ret = MPI_ratio();
			break;
		case _TPI:
			Ret = TPI_ratio();
			break;
		case _EKA:
			Ret = EKA_ratio();
			break;
		case _MKA:
			Ret = MKA_ratio();
			break;
		case _nRHO0:
			Ret = nRHO0_ratio();
			break;
		case _ERHO:
			Ret = ERHO_ratio();
			break;
		case _MRHO:
			Ret = MRHO_ratio();
			break;
		case _EKAx:
			Ret = EKAx_ratio();
			break;
		case _MKAx:
			Ret = MKAx_ratio();
			break;
		case _nOMEGA:
			Ret = nOMEGA_ratio();
			break;
		case _nETA:
			Ret = nETA_ratio();
			break;
		case _nETAi:
			Ret = nETAi_ratio();
			break;
		case _nPHI:
			Ret = nPHI_ratio();
			break;
		case _ECHARM:
			Ret = ECHARM_ratio();
			break;
		case _MuonE:
		case _MuonM:
		case _TauEE:
		case _TauET:
		case _TauMM:
		case _TauMT:
		case _PionE:
		case _PionM:
		case _KaonE:
		case _KaonM:
		case _CharmE:
		case _CharmM:
		case _CharmT:
		case _Kaon0E:
		case _Kaon0M:
		case _KaonCE:
		case _KaonCM:
			break;
		default:
			break;
	}

	return Ret;
}

double PhaseSpace::nEE_ratio()
{
	if (maxnEE_e < 0 || maxnEE_a || IsChanged())
		Max_NeutrinoLeptonAA(maxnEE_e, maxnEE_a, M_Neutrino, M_Electron);

	if (fnEE_e < 0 || fnEE_a < 0 || IsChanged())
		NeutrinoLeptonAA(fnEE_e, fnEE_a, M_Neutrino, M_Electron);

	return ( (fnEE_e + fnEE_a) * Ue(2) + fnEE_a * (Um(2) + Ut(2)) ) / 
	       ( (maxnEE_e + maxnEE_a) * Ue(2) + maxnEE_a * (Um(2) + Ut(2)) );
}

double PhaseSpace::nMM_ratio()
{
	if (maxnMM_m < 0 || maxnMM_a || IsChanged())
		Max_NeutrinoLeptonAA(maxnMM_m, maxnMM_a, M_Neutrino, M_Electron);

	if (fnMM_m < 0 || fnMM_a < 0 || IsChanged())
		NeutrinoLeptonAA(fnMM_m, fnMM_a, M_Neutrino, M_Electron);

	return ( (fnMM_m + fnMM_a) * Um(2) + fnMM_a * (Ue(2) + Ut(2)) ) / 
	       ( (maxnMM_m + maxnMM_a) * Um(2) + maxnMM_a * (Ue(2) + Ut(2)) );
}

double PhaseSpace::nEM_ratio()
{
	if (maxnEM || IsChanged())
		maxnEM = Max_NeutrinoLeptonAB(M_Neutrino, M_Electron, M_Muon);

	if (fnEM < 0 || IsChanged())
		fnEM = NeutrinoLeptonAB(M_Neutrino, M_Electron, M_Muon);

	return fnEM / maxnEM;
}

double PhaseSpace::nME_ratio()
{
	if (maxnME || IsChanged())
		maxnME = Max_NeutrinoLeptonAB(M_Neutrino, M_Muon, M_Electron);

	if (fnME < 0 || IsChanged())
		fnME = NeutrinoLeptonAB(M_Neutrino, M_Muon, M_Electron);

	return fnME / maxnME;
}

double PhaseSpace::nET_ratio()
{
	if (maxnET || IsChanged())
		maxnET = Max_NeutrinoLeptonAB(M_Neutrino, M_Electron, M_Tau);

	if (fnET < 0 || IsChanged())
		fnET = NeutrinoLeptonAB(M_Neutrino, M_Electron, M_Tau);

	return fnET / maxnET;
}

double PhaseSpace::nTE_ratio()
{
	if (maxnTE || IsChanged())
		maxnTE = Max_NeutrinoLeptonAB(M_Neutrino, M_Tau, M_Electron);

	if (fnTE < 0 || IsChanged())
		fnTE = NeutrinoLeptonAB(M_Neutrino, M_Tau, M_Electron);

	return fnET / maxnET;
}

double PhaseSpace::nMT_ratio()
{
	if (maxnMT || IsChanged())
		maxnMT = Max_NeutrinoLeptonAB(M_Neutrino, M_Muon, M_Tau);

	if (fnMT < 0 || IsChanged())
		fnMT = NeutrinoLeptonAB(M_Neutrino, M_Muon, M_Tau);

	return fnMT / maxnMT;
}

double PhaseSpace::nTM_ratio()
{
	if (maxnTM || IsChanged())
		maxnTM = Max_NeutrinoLeptonAB(M_Neutrino, M_Tau, M_Muon);

	if (fnTM < 0 || IsChanged())
		fnTM = NeutrinoLeptonAB(M_Neutrino, M_Tau, M_Muon);

	return fnET / maxnET;
}

double PhaseSpace::nPI0_ratio()
{
	if (maxnPI0 || IsChanged())
		maxnPI0 = Max_NeutrinoPseudoMeson(M_Neutrino, M_Pion0);

	if (fnPI0 < 0 || IsChanged())
		fnPI0 = NeutrinoPseudoMeson(M_Neutrino, M_Pion0);

	return fnPI0 / maxnPI0;
}

double PhaseSpace::nETA_ratio()
{
	if (maxnETA || IsChanged())
		maxnETA = Max_NeutrinoPseudoMeson(M_Neutrino, M_Eta);

	if (fnETA < 0 || IsChanged())
		fnETA = NeutrinoPseudoMeson(M_Neutrino, M_Eta);

	return fnETA / maxnETA;
}

double PhaseSpace::nETAi_ratio()
{
	if (maxnETAi || IsChanged())
		maxnETAi = Max_NeutrinoPseudoMeson(M_Neutrino, M_Etai);

	if (fnETAi < 0 || IsChanged())
		fnETAi = NeutrinoPseudoMeson(M_Neutrino, M_Etai);

	return fnETAi / maxnETAi;
}

double PhaseSpace::EPI_ratio()
{
	if (maxEPI || IsChanged())
		maxEPI = Max_LeptonPseudoMeson(M_Electron, M_Pion);

	if (fEPI < 0 || IsChanged())
		fEPI = LeptonPseudoMeson(M_Electron, M_Pion);

	return fEPI / maxEPI;
}

double PhaseSpace::MPI_ratio()
{
	if (maxMPI || IsChanged())
		maxMPI = Max_LeptonPseudoMeson(M_Muon, M_Pion);

	if (fMPI < 0 || IsChanged())
		fMPI = LeptonPseudoMeson(M_Muon, M_Pion);

	return fMPI / maxMPI;
}

double PhaseSpace::TPI_ratio()
{
	if (maxTPI || IsChanged())
		maxTPI = Max_LeptonPseudoMeson(M_Tau, M_Pion);

	if (fTPI < 0 || IsChanged())
		fTPI = LeptonPseudoMeson(M_Tau, M_Pion);

	return fTPI / maxTPI;
}

double PhaseSpace::EKA_ratio()
{
	if (maxEKA || IsChanged())
		maxEKA = Max_LeptonPseudoMeson(M_Electron, M_Kaon);

	if (fEKA < 0 || IsChanged())
		fEKA = LeptonPseudoMeson(M_Electron, M_Kaon);

	return fEKA / maxEKA;
}

double PhaseSpace::MKA_ratio()
{
	if (maxMKA || IsChanged())
		maxMKA = Max_LeptonPseudoMeson(M_Muon, M_Kaon);

	if (fMKA < 0 || IsChanged())
		fMKA = LeptonPseudoMeson(M_Muon, M_Kaon);

	return fMKA / maxMKA;
}

double PhaseSpace::ECHARM_ratio()
{
	if (maxECHARM || IsChanged())
		maxECHARM = Max_LeptonPseudoMeson(M_Electron, M_Charm);

	if (fECHARM < 0 || IsChanged())
		fECHARM = LeptonPseudoMeson(M_Electron, M_Charm);

	return fECHARM / maxECHARM;
}

double PhaseSpace::nRHO0_ratio()
{
	if (maxnRHO0 || IsChanged())
		maxnRHO0 = Max_NeutrinoVectorMeson(M_Neutrino, M_Rho0);

	if (fnRHO0 < 0 || IsChanged())
		fnRHO0 = NeutrinoVectorMeson(M_Neutrino, M_Rho0);

	return fnRHO0 / maxnRHO0;
}

double PhaseSpace::nOMEGA_ratio()
{
	if (maxnOMEGA || IsChanged())
		maxnOMEGA = Max_NeutrinoVectorMeson(M_Neutrino, M_Omega);

	if (fnOMEGA < 0 || IsChanged())
		fnOMEGA = NeutrinoVectorMeson(M_Neutrino, M_Omega);

	return fnOMEGA / maxnOMEGA;
}

double PhaseSpace::nPHI_ratio()
{
	if (maxnPHI || IsChanged())
		maxnPHI = Max_NeutrinoVectorMeson(M_Neutrino, M_Phi);

	if (fnPHI < 0 || IsChanged())
		fnPHI = NeutrinoVectorMeson(M_Neutrino, M_Phi);

	return fnPHI / maxnPHI;
}

double PhaseSpace::ERHO_ratio()
{
	if (maxERHO || IsChanged())
		maxERHO = Max_LeptonVectorMeson(M_Electron, M_Rho);

	if (fERHO < 0 || IsChanged())
		fERHO = LeptonVectorMeson(M_Electron, M_Rho);

	return fERHO / maxERHO;
}

double PhaseSpace::MRHO_ratio()
{
	if (maxMRHO || IsChanged())
		maxMRHO = Max_LeptonVectorMeson(M_Muon, M_Rho);

	if (fMRHO < 0 || IsChanged())
		fMRHO = LeptonVectorMeson(M_Muon, M_Rho);

	return fMRHO / maxMRHO;
}

double PhaseSpace::EKAx_ratio()
{
	if (maxEKAx || IsChanged())
		maxEKAx = Max_LeptonVectorMeson(M_Electron, M_Kaonx);

	if (fEKAx < 0 || IsChanged())
		fEKAx = LeptonVectorMeson(M_Electron, M_Kaonx);

	return fEKAx / maxEKAx;
}

double PhaseSpace::MKAx_ratio()
{
	if (maxMKAx || IsChanged())
		maxMKAx = Max_LeptonVectorMeson(M_Muon, M_Kaonx);

	if (fMKAx < 0 || IsChanged())
		fMKAx = LeptonVectorMeson(M_Muon, M_Kaonx);

	return fMKAx / maxMKAx;
}


//////////////////////
//GENERIC DIFF WIDTH//
//////////////////////
//
double PhaseSpace::NeutrinoLeptonAA(double &d_Ul, double &d_Un, double M_Neut, double M_Lepton)
{
	double dMN2 = M_Neut*M_Neut/Mass(2);
	double dML2 = M_Lepton*M_Lepton/Mass(2);

	double s, t, cos0, cos1;
	Kinematic_3B(s, t, cos0, cos1);

	double gL_CC = 2 - 4*GetParticle();	//times U(lepton flavour)
	double gR_CC = 0;
	double gL_NC = -1 + 2*Const::fSin2W;	//times U(neutrino flavour)
	double gR_NC = 2*Const::fSin2W;

	d_Ul = NeutrinoLeptonLepton(dMN2, dML2, dML2, gL_CC, gR_CC, s, t, cos0, cos1);
	d_Un = NeutrinoLeptonLepton(dMN2, dML2, dML2, gL_NC, gR_NC, s, t, cos0, cos1);

	return 0.0;
}

double PhaseSpace::Max_NeutrinoLeptonAA(double &max_Ul, double &max_Ua, double M_Neut, double M_Lepton)
{
	double dMN2 = M_Neut*M_Neut/Mass(2);
	double dML2 = M_Lepton*M_Lepton/Mass(2);

	double gL_CC = 2 - 4*GetParticle();	//times U(lepton flavour)
	double gR_CC = 0;
	double gL_NC = -1 + 2*Const::fSin2W;	//times U(neutrino flavour)
	double gR_NC = 2*Const::fSin2W;

	max_Ul = max_NeutrinoLeptonLepton(dMN2, dML2, dML2, gL_CC, gR_CC);
	max_Ua = max_NeutrinoLeptonLepton(dMN2, dML2, dML2, gL_NC, gR_NC);

	return 0.0;
}

double PhaseSpace::NeutrinoLeptonAB(double M_Neut, double M_LeptonA, double M_LeptonB)
{
	double dMN2 = M_Neut*M_Neut/Mass(2);
	double dMA2 = M_LeptonA*M_LeptonB/Mass(2);
	double dMB2 = M_LeptonB*M_LeptonB/Mass(2);

	double s, t, cos0, cos1;
	Kinematic_3B(s, t, cos0, cos1);

	double gL = 2;	//times U(lepton flavour)
	double gR = 0;

	return NeutrinoLeptonLepton(dMN2, dMA2, dMB2, gL, gR, s, t, cos0, cos1);

	return 0.0;
}

double PhaseSpace::Max_NeutrinoLeptonAB(double M_Neut, double M_LeptonA, double M_LeptonB)
{
	double dMN2 = M_Neut*M_Neut/Mass(2);
	double dMA2 = M_LeptonA*M_LeptonB/Mass(2);
	double dMB2 = M_LeptonB*M_LeptonB/Mass(2);

	double gL = 2;	//times U(lepton flavour)
	double gR = 0;

	return max_NeutrinoLeptonLepton(dMN2, dMA2, dMB2, gL, gR);
}

double PhaseSpace::max_NeutrinoLeptonLepton(double x, double y, double z, double gL, double gR)
{
	F_var.clear();

	F_var.push_back(x);	//0
	F_var.push_back(y);	//1
	F_var.push_back(z);	//2
	F_var.push_back(gL);	//3
	F_var.push_back(gR);	//4

	std::vector<double> xMin;
	SetFunction_D(&PhaseSpace::max_NeutrinoLeptonLepton_D);
	return Inte::NelMedSolver(this, xMin, 4, Inte::Max);
}

double PhaseSpace::max_NeutrinoLeptonLepton_D(double *p)
{
	double &x  = F_var.at(0);
	double &y  = F_var.at(1);
	double &z  = F_var.at(2);
	double &gL = F_var.at(3);
	double &gR = F_var.at(4);

	double &s    = p[0];
	double &t    = p[1];
	double &cos0 = p[2];
	double &cos1 = p[3];

	return NeutrinoLeptonLepton(x, y, z, s, t, cos0, cos1, gL, gR);
}

double PhaseSpace::NeutrinoLeptonLepton(double x, double y, double z, double s, double t, double cos0, double cos1, double gL, double gR)
{
	double M2 = (gL * gL + gR * gR) * M2_WW(x, y, z, s, cos0) +
		    (2 * gL * gR) * M2_WZ(x, y, z, s, t, cos0, cos1);

	return dGammad5_3B(M2);
}

double PhaseSpace::NeutrinoPseudoMeson(double M_Neut, double M_Meson)
{
	double dMN2 = M_Neut*M_Neut/Mass(2);
	double dMM2 = M_Meson*M_Meson/Mass(2);

	double cos0;
	Kinematic_2B(cos0);

	return 2 * LeptonPseudo(dMN2, dMM2, cos0);
}	//     2 is factor from decay constant which is sqrt(2) wrt to charged meson

double PhaseSpace::Max_NeutrinoPseudoMeson(double M_Neut, double M_Meson)
{
	double dMN2 = M_Neut*M_Neut/Mass(2);
	double dMM2 = M_Meson*M_Meson/Mass(2);

	return 2 * max_LeptonPseudo(dMN2, dMM2);
}

double PhaseSpace::LeptonPseudoMeson(double M_Neut, double M_Meson)
{
	double dMN2 = M_Neut*M_Neut/Mass(2);
	double dMM2 = M_Meson*M_Meson/Mass(2);

	double cos0;
	Kinematic_2B(cos0);

	return LeptonPseudo(dMN2, dMM2, cos0);
}

double PhaseSpace::Max_LeptonPseudoMeson(double M_Neut, double M_Meson)
{
	double dMN2 = M_Neut*M_Neut/Mass(2);
	double dMM2 = M_Meson*M_Meson/Mass(2);

	return max_LeptonPseudo(dMN2, dMM2);
}

double PhaseSpace::max_LeptonPseudo(double x, double y)
{
	F_var.clear();

	F_var.push_back(x);
	F_var.push_back(y);

	SetFunction(&PhaseSpace::max_LeptonPseudo_cos0);
	return Inte::GoldRatioSolver(this, Inte::Max);
}

double PhaseSpace::max_LeptonPseudo_cos0(double cos0)
{
	double &x  = F_var.at(0);
	double &y  = F_var.at(1);

	return LeptonPseudo(x, y, cos0);
}

double PhaseSpace::LeptonPseudo(double x, double y, double cos0)
{
	double M2 = M2_LeptonPseudoMeson(x, y, cos0);
	return dGammad2_2B(M2, x, y); 
}

double PhaseSpace::NeutrinoVectorMeson(double M_Neut, double M_Meson)
{
	double dMN2 = M_Neut*M_Neut/Mass(2);
	double dMM2 = M_Meson*M_Meson/Mass(2);

	double cos0;
	Kinematic_2B(cos0);

	return 2 * LeptonVector(dMN2, dMM2, cos0);
}	//     2 is factor from decay constant which is sqrt(2) wrt to charged meson

double PhaseSpace::Max_NeutrinoVectorMeson(double M_Neut, double M_Meson)
{
	double dMN2 = M_Neut*M_Neut/Mass(2);
	double dMM2 = M_Meson*M_Meson/Mass(2);

	return 2 * max_LeptonVector(dMN2, dMM2);
}

double PhaseSpace::LeptonVectorMeson(double M_Neut, double M_Meson)
{
	double dMN2 = M_Neut*M_Neut/Mass(2);
	double dMM2 = M_Meson*M_Meson/Mass(2);

	double cos0;
	Kinematic_2B(cos0);

	return LeptonVector(dMN2, dMM2, cos0);
}

double PhaseSpace::Max_LeptonVectorMeson(double M_Neut, double M_Meson)
{
	double dMN2 = M_Neut*M_Neut/Mass(2);
	double dMM2 = M_Meson*M_Meson/Mass(2);

	return max_LeptonVector(dMN2, dMM2);
}

double PhaseSpace::max_LeptonVector(double x, double y)
{
	F_var.clear();

	F_var.push_back(x);
	F_var.push_back(y);

	SetFunction(&PhaseSpace::max_LeptonVector_cos0);
	return Inte::GoldRatioSolver(this, Inte::Max);
}

double PhaseSpace::max_LeptonVector_cos0(double cos0)
{
	double &x  = F_var.at(0);
	double &y  = F_var.at(1);

	return LeptonVector(x, y, cos0);
}

double PhaseSpace::LeptonVector(double x, double y, double cos0)
{
	double M2 = M2_LeptonVectorMeson(x, y, cos0);
	return dGammad2_2B(M2, x, y); 
}


void PhaseSpace::Kinematic_2B(double &cos0)
{
	TLorentzVector vec0 = Daughter(0);
	TLorentzVector vec1 = Daughter(1);

	cos0 = vec1.CosTheta();
}

void PhaseSpace::Kinematic_3B(double &s, double &t, double &cos0, double &cos1)
{
	TLorentzVector vec0 = Daughter(0);
	TLorentzVector vec1 = Daughter(1);
	TLorentzVector vec2 = Daughter(2);

	TLorentzVector vec_s = *Rest() - vec1;
	TLorentzVector vec_t = *Rest() - vec2;

	s = vec_s.M2();
	t = vec_t.M2();

	cos0 = vec1.CosTheta();
	cos1 = vec2.CosTheta();
}

//////////
//	//
//	//
//////////

TLorentzVector PhaseSpace::Daughter(unsigned int i)
{
	TLorentzVector D_vec;

	if (i < Daughter())
	{
		TLorentzVector D_vec = *Event->GetDecay(i);
		D_vec.Boost(LabF()->BoostVector());
	}

	return D_vec;
}

unsigned int PhaseSpace::Daughter()
{
	return vMass.size();
}

TLorentzVector* PhaseSpace::Rest()
{
	return N_rest;
}

TLorentzVector* PhaseSpace::LabF()
{
	return N_labf;
}

void PhaseSpace::SetNLabf(TLorentzVector &Vec)
{
	N_labf = &Vec;
	SetNRest(N_labf->M());
}

void PhaseSpace::SetNRest(double Mass)
{
	N_rest->SetPxPyPzE(0, 0, 0, Mass);
}

void PhaseSpace::Reset()
{
	fnnn    = -1.0;
	fnGAMMA = -1.0;
	fnEE_e  = -1.0;
	fnEE_a  = -1.0;
	fnEM    = -1.0;
	fnME    = -1.0;
	fnMM_m  = -1.0;
	fnMM_a  = -1.0;
	fnET    = -1.0;
	fnTE    = -1.0;
	fnMT    = -1.0;
	fnTM    = -1.0;
	fnPI0   = -1.0;
	fEPI    = -1.0;
	fMPI    = -1.0;
	fTPI    = -1.0;
	fEKA    = -1.0;
	fMKA    = -1.0;
	fnRHO0  = -1.0;
	fERHO   = -1.0;
	fMRHO   = -1.0;
	fEKAx   = -1.0;
	fMKAx   = -1.0;
	fnETA   = -1.0;
	fnETAi  = -1.0;
	fnOMEGA = -1.0;
	fnPHI   = -1.0;
	fECHARM = -1.0;

	maxnnn    = -1.0;
	maxnGAMMA = -1.0;
	maxnEE_e  = -1.0;
	maxnEE_a  = -1.0;
	maxnEM    = -1.0;
	maxnME    = -1.0;
	maxnMM_m  = -1.0;
	maxnMM_a  = -1.0;
	maxnET    = -1.0;
	maxnTE    = -1.0;
	maxnMT    = -1.0;
	maxnTM    = -1.0;
	maxnPI0   = -1.0;
	maxEPI    = -1.0;
	maxMPI    = -1.0;
	maxTPI    = -1.0;
	maxEKA    = -1.0;
	maxMKA    = -1.0;
	maxnRHO0  = -1.0;
	maxERHO   = -1.0;
	maxMRHO   = -1.0;
	maxEKAx   = -1.0;
	maxMKAx   = -1.0;
	maxnETA   = -1.0;
	maxnETAi  = -1.0;
	maxnOMEGA = -1.0;
	maxnPHI   = -1.0;
	maxECHARM = -1.0;
}

void PhaseSpace::SetFunction(double (PhaseSpace::*FF)(double))
{
	double (PhaseSpace::*Function)(double) = 
		static_cast<double (Amplitude::*)(double)>(FF); // ok!
	SetFunction(Function);
}

void PhaseSpace::SetFunction_D(double (PhaseSpace::*FF)(double*))
{
	double (PhaseSpace::*Function)(double*) = 
		static_cast<double (Amplitude::*)(double*)>(FF); // ok!
	SetFunction_D(Function);
}
