#include "PhaseSpace.h"

PhaseSpace::PhaseSpace() :
	Event(new TGenPhaseSpace),
	GenMT(new TRandom3(0)),
	P_labf(new TLorentzVector()),
	P_rest(new TLorentzVector())
{
}

bool PhaseSpace::SetDecay(Channel Name)
{
	bool Type;
	if (Channel_prev != Name)
	{
		Type = LoadMass(Name);
		Channel_prev = Name;
	}

	double MassArray[10];
	if (Type)			//it is decay, first entry is lightest child
		MassArray[0] = vMass.at(0);
	else				//it is production, first entry is neutrino
		MassArray[0] = MassN();

	for (unsigned int i = 1; i < Daughter(); ++i)
		MassArray[i] = vMass.at(i);

	return Event->SetDecay(*Parent(RestFrame), Daughter(), MassArray);
}

bool PhaseSpace::Generate(Channel Name)
{
	unsigned int cc = 0;
	if (SetDecay(Name))
	{
		do
		{
			++cc;
			Event->Generate();
		}
		while (GenMT->Rndm() < Ratio(Name));

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
		case _MuonE:	//these are needed as well
			Ret = MuonE_ratio();
			break;
		case _MuonM:
			Ret = MuonM_ratio();
			break;
		case _TauEE:
			Ret = TauEE_ratio();
			break;
		case _TauET:
			Ret = TauET_ratio();
			break;
		case _TauMM:
			Ret = TauMM_ratio();
			break;
		case _TauMT:
			Ret = TauMT_ratio();
			break;
		case _TauPI:				//1
			//Ret = TauPI();
			break;
		case _Tau2PI:				//1
			//Ret = Tau2PI();
			break;
		case _PionE:				//1
			//Ret = PionE_ratio();
			break;
		case _PionM:				//1
			//Ret = PionM_ratio();
			break;
		case _KaonE:				//1
			//Ret = KaonE_ratio();
			break;
		case _KaonM:				//1
			//Ret = KaonM_ratio();
			break;
		case _CharmE:				//1
			//Ret = CharmE_ratio();
			break;
		case _CharmM:				//1
			//Ret = CharmE_ratio();
			break;
		case _CharmT:				//1
			//Ret = CharmT_ratio();
			break;
		case _Kaon0E:
			Ret = Kaon0E_ratio();
			break;
		case _Kaon0M:
			Ret = Kaon0M_ratio();
			break;
		case _KaonCE:
			Ret = KaonCE_ratio();
			break;
		case _KaonCM:
			Ret = KaonCM_ratio();
			break;
		default:
			break;
	}

	return Ret;
}

//Neutrino LeptonLepton AA
double PhaseSpace::nEE_ratio()
{
	return NeutrinoLeptonAA_ratio(maxnEE_e, maxnEE_a, &PhaseSpace::Ue);
}

double PhaseSpace::nMM_ratio()
{
	return NeutrinoLeptonAA_ratio(maxnMM_m, maxnMM_a, &PhaseSpace::Um);
}

double PhaseSpace::NeutrinoLeptonAA_ratio(double &maxName_u, double &maxName_a, double (PhaseSpace::*Uu)(int))
{
	if (maxName_u < 0 || maxName_a < 0 || IsChanged())
		Max_NeutrinoLeptonAA(maxName_u, maxName_a, vMass.at(0), vMass.at(1));

	double fName_u, fName_a;
	NeutrinoLeptonAA(fName_u, fName_a, vMass.at(0), vMass.at(1));

	return ( fName_u   * (this->*Uu)(2) + fName_a   * (Ue(2) + Um(2) + Ut(2)) ) / 
	       ( maxName_u * (this->*Uu)(2) + maxName_a * (Ue(2) + Um(2) + Ut(2)) );
}

//Neutrino LeptonLepton AB
double PhaseSpace::nEM_ratio()
{
	return NeutrinoLeptonAB_ratio(maxnEM);
}

double PhaseSpace::nME_ratio()
{
	return NeutrinoLeptonAB_ratio(maxnME);
}

double PhaseSpace::nET_ratio()
{
	return NeutrinoLeptonAB_ratio(maxnET);
}

double PhaseSpace::nTE_ratio()
{
	return NeutrinoLeptonAB_ratio(maxnTE);
}

double PhaseSpace::nMT_ratio()
{
	return NeutrinoLeptonAB_ratio(maxnMT);
}

double PhaseSpace::nTM_ratio()
{
	return NeutrinoLeptonAB_ratio(maxnTM);
}

double PhaseSpace::NeutrinoLeptonAB_ratio(double &maxName)
{
	if (maxName < 0 || IsChanged())
		maxName = Max_NeutrinoLeptonAB(vMass.at(0), vMass.at(1), vMass.at(2));

	return NeutrinoLeptonAB(vMass.at(0), vMass.at(1), vMass.at(2)) / maxName;
}

//neutrino psuedomeson

double PhaseSpace::nPI0_ratio()
{
	return NeutrinoPseudoMeson_ratio(maxnPI0);
}

double PhaseSpace::nETA_ratio()
{
	return NeutrinoPseudoMeson_ratio(maxnETA);
}

double PhaseSpace::nETAi_ratio()
{
	return NeutrinoPseudoMeson_ratio(maxnETAi);
}

double PhaseSpace::NeutrinoPseudoMeson_ratio(double &maxName)
{
	if (maxName < 0 || IsChanged())
		maxName = Max_NeutrinoPseudoMeson(vMass.at(0), vMass.at(1));

	return NeutrinoPseudoMeson(vMass.at(0), vMass.at(1)) / maxName;
}

//lepton psuedomeson

double PhaseSpace::EPI_ratio()
{
	return LeptonPseudoMeson_ratio(maxEPI);
}

double PhaseSpace::MPI_ratio()
{
	return LeptonPseudoMeson_ratio(maxMPI);
}

double PhaseSpace::TPI_ratio()
{
	return LeptonPseudoMeson_ratio(maxTPI);
}

double PhaseSpace::EKA_ratio()
{
	return LeptonPseudoMeson_ratio(maxEKA);
}

double PhaseSpace::MKA_ratio()
{
	return LeptonPseudoMeson_ratio(maxMKA);
}

double PhaseSpace::ECHARM_ratio()
{
	return LeptonPseudoMeson_ratio(maxECHARM);
}

double PhaseSpace::LeptonPseudoMeson_ratio(double &maxName)
{
	if (maxName < 0 || IsChanged())
		maxName = Max_LeptonPseudoMeson(vMass.at(0), vMass.at(1));

	return LeptonPseudoMeson(vMass.at(0), vMass.at(1)) / maxName;
}

//neutrino vector meson

double PhaseSpace::nRHO0_ratio()
{
	return NeutrinoVectorMeson_ratio(maxnRHO0);
}

double PhaseSpace::nOMEGA_ratio()
{
	return NeutrinoVectorMeson_ratio(maxnOMEGA);
}

double PhaseSpace::nPHI_ratio()
{
	return NeutrinoVectorMeson_ratio(maxnPHI);
}

double PhaseSpace::NeutrinoVectorMeson_ratio(double &maxName)
{
	if (maxName < 0 || IsChanged())
		maxName = Max_NeutrinoVectorMeson(vMass.at(0), vMass.at(1));

	return NeutrinoVectorMeson(vMass.at(0), vMass.at(1)) / maxName;
}

//lepton vector meson

double PhaseSpace::ERHO_ratio()
{
	return LeptonVectorMeson_ratio(maxERHO);
}

double PhaseSpace::MRHO_ratio()
{
	return LeptonVectorMeson_ratio(maxMRHO);
}

double PhaseSpace::EKAx_ratio()
{
	return LeptonVectorMeson_ratio(maxEKAx);
}

double PhaseSpace::MKAx_ratio()
{
	return LeptonVectorMeson_ratio(maxMKAx);
}

double PhaseSpace::LeptonVectorMeson_ratio(double &maxName)
{
	if (maxName < 0 || IsChanged())
		maxName = Max_LeptonVectorMeson(vMass.at(0), vMass.at(1));

	return LeptonVectorMeson(vMass.at(0), vMass.at(1)) / maxName;
}

////PRODUCTION
//
double PhaseSpace::MuonE_ratio()
{
	if (maxMuonE < 0 || IsChanged())
		maxMuonE = Max_AntiLeptonNeutrino(M_Muon, M_Electron, M_Neutrino);

	return AntiLeptonNeutrino(M_Muon, M_Electron, M_Neutrino) / maxMuonE;
}

double PhaseSpace::MuonM_ratio()
{
	if (maxMuonM < 0 || IsChanged())
		maxMuonM = Max_LeptonNeutrino(M_Muon, M_Electron, M_Neutrino);

	return LeptonNeutrino(M_Muon, M_Electron, M_Neutrino) / maxMuonM;
}


double PhaseSpace::TauEE_ratio()
{
	if (maxTauEE < 0 || IsChanged())
		maxTauEE = Max_AntiLeptonNeutrino(M_Tau, M_Electron, M_Neutrino);

	return AntiLeptonNeutrino(M_Tau, M_Electron, M_Neutrino) / maxTauEE;
}

double PhaseSpace::TauET_ratio()
{
	if (maxTauET < 0 || IsChanged())
		maxTauET = Max_LeptonNeutrino(M_Tau, M_Electron, M_Neutrino);

	return LeptonNeutrino(M_Tau, M_Electron, M_Neutrino) / maxTauET;
}

double PhaseSpace::TauMM_ratio()
{
	if (maxTauMM < 0 || IsChanged())
		maxTauMM = Max_AntiLeptonNeutrino(M_Tau, M_Muon, M_Neutrino);

	return AntiLeptonNeutrino(M_Tau, M_Muon, M_Neutrino) / maxTauMM;
}

double PhaseSpace::TauMT_ratio()
{
	if (maxTauMT < 0 || IsChanged())
		maxTauMT = Max_LeptonNeutrino(M_Tau, M_Muon, M_Neutrino);

	return LeptonNeutrino(M_Tau, M_Muon, M_Neutrino) / maxTauMT;
}

/*
double PhaseSpace::TauPI_ratio()
{
	if (maxTauPI || IsChanged())
		maxTauPI = Max_LeptonMeson(M_Tau, M_Pion);

	return LeptonMeson(M_Tau, M_Pion) / maxTauPI;
}

double PhaseSpace::LeptonMesonDecay(double M_Lepton, double M_Meson)
{
	SetMass(M_Lepton);
	double dMN2 = MassN(2)/Mass(2);
	double dMM2 = M_Meson*M_Meson/Mass(2);

	double M2 = M2_LeptonMeson(dMN2, dMM2);
	return dGammad2_2B(M2, dMN2, dMM2);
}

double PhaseSpace::Max_LeptonMeson(double M_Lepton, double M_Meson)
{
	SetMass(M_Lepton);
	double M2 = M2_LeptonMeson(x, y);
	return dGammad2_2B(M2, x, y);
}

double PhaseSpace::PionE_ratio()
{
	if (maxPionE || IsChanged())
		maxPionE = Max_MesonTwo(M_Pion, M_Electron);

	return MesonTwo(M_Pion, M_Electron) / maxPionE;
}

double PhaseSpace::PionM_ratio()
{
	if (maxPionM || IsChanged())
		maxPionM = Max_MesonTwo(M_Pion, M_Electron);

	return MesonTwo(M_Pion, M_Electron) / maxPionM;
}

double PhaseSpace::KaonE_ratio()
{
	if (maxKaonE || IsChanged())
		maxKaonE = Max_MesonTwo(M_Kaon, M_Electron);

	return MesonTwo(M_Kaon, M_Electron) / maxKaonE;
}

double PhaseSpace::KaonM_ratio()
{
	if (maxKaonM || IsChanged())
		maxKaonM = Max_MesonTwo(M_Kaon, M_Muon);

	return MesonTwo(M_Kaon, M_Muon) / maxKaonM;
}

double PhaseSpace::CharmE_ratio()
{
	if (maxCharmE || IsChanged())
		maxCharmE = Max_MesonTwo(M_CharmS, M_Electron);

	return MesonTwo(M_CharmS, M_Electron) / maxCharmE;
}

double PhaseSpace::CharmM_ratio()
{
	if (maxCharmM || IsChanged())
		maxCharmM = Max_MesonTwo(M_CharmS, M_Muon);

	return MesonTwo(M_CharmS, M_Muon) / maxCharmM;
}

double PhaseSpace::CharmT_ratio()
{
	if (maxCharmT || IsChanged())
		maxCharmT = Max_MesonTwo(M_CharmS, M_Tau);

	return MesonTwo(M_CharmS, M_Tau) / maxCharmT;
}
*/

double PhaseSpace::Kaon0E_ratio()
{
	if (maxKaon0E < 0 || IsChanged())
		maxKaon0E = Max_MesonThree(M_Kaon0, M_Pion, M_Electron, Const::fK0L_, Const::fK0L0);

	return MesonThree(M_Kaon0, M_Pion, M_Electron, Const::fK0L_, Const::fK0L0) / maxKaon0E;
}

double PhaseSpace::Kaon0M_ratio()
{
	if (maxKaon0M < 0 || IsChanged())
		maxKaon0M = Max_MesonThree(M_Kaon0, M_Pion, M_Muon, Const::fK0L_, Const::fK0L0);

	return MesonThree(M_Kaon0, M_Pion, M_Muon, Const::fK0L_, Const::fK0L0) / maxKaon0M;
}

double PhaseSpace::KaonCE_ratio()
{
	if (maxKaonCE < 0 || IsChanged())
		maxKaonCE = Max_MesonThree(M_Kaon, M_Pion0, M_Electron, Const::fKCL_, Const::fKCL0);

	return MesonThree(M_Kaon, M_Pion0, M_Electron, Const::fKCL_, Const::fKCL0) / maxKaonCE;
}

double PhaseSpace::KaonCM_ratio()
{
	if (maxKaonCM < 0 || IsChanged())
		maxKaonCM = Max_MesonThree(M_Kaon, M_Pion0, M_Muon, Const::fKCL_, Const::fKCL0);

	return MesonThree(M_Kaon, M_Pion0, M_Muon, Const::fKCL_, Const::fKCL0) / maxKaonCM;
}


///////////////
//DECAY MODES//
///////////////
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

double PhaseSpace::max_NeutrinoLeptonLepton_D(const double *p)
{
	double x  = F_var.at(0);
	double y  = F_var.at(1);
	double z  = F_var.at(2);
	double gL = F_var.at(3);
	double gR = F_var.at(4);

	double s_    = p[0];
	double t_    = p[1];
	double cos0_ = p[2];
	double cos1_ = p[3];

	return NeutrinoLeptonLepton(s_, t_, cos0_, cos1_, x, y, z, gL, gR);
}

double PhaseSpace::NeutrinoLeptonLepton(double s, double t, double cos0, double cos1, double x, double y, double z, double gL, double gR)
{
	double M2 = (gL * gL + gR * gR) * M2_WW(x, y, z, s, cos0) +
		    (2 * gL * gR) * M2_WZ(x, y, z, s, t, cos0, cos1);

	return dGammad5_3B(M2);
}

/////lepton and pseudomeson
////phase space is in return of ToPsuedonMeson/ToVectorMeson
//
double PhaseSpace::NeutrinoPseudoMeson(double M_Neut, double M_Meson)
{
	SetMass(MassN());
	double dMN2 = M_Neut*M_Neut/Mass(2);
	double dMM2 = M_Meson*M_Meson/Mass(2);

	double cos0;
	Kinematic_2B(cos0);

	M2_Function = &Amplitude::M2_NeutrinoPseudoMeson;
	return ToPseudoMeson(dMN2, dMM2, cos0);
}	//     2 is factor from decay constant which is sqrt(2) wrt to charged meson

double PhaseSpace::Max_NeutrinoPseudoMeson(double M_Neut, double M_Meson)
{
	SetMass(MassN());
	double dMN2 = M_Neut*M_Neut/Mass(2);
	double dMM2 = M_Meson*M_Meson/Mass(2);

	M2_Function = &Amplitude::M2_NeutrinoPseudoMeson;
	return max_ToPseudoMeson(dMN2, dMM2);
}

double PhaseSpace::LeptonPseudoMeson(double M_Lepton, double M_Meson)
{
	SetMass(MassN());
	double dML2 = M_Lepton*M_Lepton/Mass(2);
	double dMM2 = M_Meson*M_Meson/Mass(2);

	double cos0;
	Kinematic_2B(cos0);

	M2_Function = &Amplitude::M2_LeptonPseudoMeson;
	return ToPseudoMeson(dML2, dMM2, cos0);
}

double PhaseSpace::Max_LeptonPseudoMeson(double M_Lepton, double M_Meson)
{
	SetMass(MassN());
	double dML2 = M_Lepton*M_Lepton/Mass(2);
	double dMM2 = M_Meson*M_Meson/Mass(2);

	M2_Function = &Amplitude::M2_LeptonPseudoMeson;
	return max_ToPseudoMeson(dML2, dMM2);
}

double PhaseSpace::max_ToPseudoMeson(double x, double y)
{
	F_var.clear();

	F_var.push_back(x);
	F_var.push_back(y);

	SetFunction(&PhaseSpace::max_ToPseudoMeson_cos0);
	return Inte::GoldRatioSolver(this, Inte::Max);
}

double PhaseSpace::max_ToPseudoMeson_cos0(const double cos0)
{
	double x = F_var.at(0);
	double y = F_var.at(1);
	double cos0_ = -1 + 2 * cos0;

	return ToPseudoMeson(x, y, cos0_);
}

double PhaseSpace::ToPseudoMeson(double x, double y, double cos0)
{
	double M2 = (this->*M2_Function)(x, y, cos0);	//either M2_NeutrinoPseudoMeson 
	return dGammad2_2B(M2, x, y);			//or	 M2_LeptonPseudoMeson
}

/////vector meson
//
double PhaseSpace::NeutrinoVectorMeson(double M_Neut, double M_Meson)
{
	SetMass(MassN());
	double dMN2 = M_Neut*M_Neut/Mass(2);
	double dMM2 = M_Meson*M_Meson/Mass(2);

	double cos0;
	Kinematic_2B(cos0);

	M2_Function = &Amplitude::M2_NeutrinoVectorMeson;
	return ToVectorMeson(dMN2, dMM2, cos0);
}

double PhaseSpace::Max_NeutrinoVectorMeson(double M_Neut, double M_Meson)
{
	SetMass(MassN());
	double dMN2 = M_Neut*M_Neut/Mass(2);
	double dMM2 = M_Meson*M_Meson/Mass(2);

	M2_Function = &Amplitude::M2_NeutrinoVectorMeson;
	return max_ToVectorMeson(dMN2, dMM2);
}

double PhaseSpace::LeptonVectorMeson(double M_Neut, double M_Meson)
{
	SetMass(MassN());
	double dMN2 = M_Neut*M_Neut/Mass(2);
	double dMM2 = M_Meson*M_Meson/Mass(2);

	double cos0;
	Kinematic_2B(cos0);

	M2_Function = &Amplitude::M2_LeptonVectorMeson;
	return ToVectorMeson(dMN2, dMM2, cos0);
}

double PhaseSpace::Max_LeptonVectorMeson(double M_Neut, double M_Meson)
{
	SetMass(MassN());
	double dMN2 = M_Neut*M_Neut/Mass(2);
	double dMM2 = M_Meson*M_Meson/Mass(2);

	M2_Function = &Amplitude::M2_LeptonVectorMeson;
	return max_ToVectorMeson(dMN2, dMM2);
}

double PhaseSpace::max_ToVectorMeson(double x, double y)
{
	F_var.clear();

	F_var.push_back(x);
	F_var.push_back(y);

	SetFunction(&PhaseSpace::max_ToVectorMeson_cos0);
	return Inte::GoldRatioSolver(this, Inte::Max);
}

double PhaseSpace::max_ToVectorMeson_cos0(const double cos0)
{
	double x = F_var.at(0);
	double y = F_var.at(1);
	double cos0_ = -1 + 2 * cos0;

	std::cout << "cos0 " << cos0_ << std::endl;
	return ToVectorMeson(x, y, cos0_);
}

double PhaseSpace::ToVectorMeson(double x, double y, double cos0)
{
	double M2 = (this->*M2_Function)(x, y, cos0);
	return dGammad2_2B(M2, x, y); 
}

//////////////
//PRODUCTION//
//////////////
//
double PhaseSpace::LeptonNeutrino(double M_Lepton0, double M_Lepton, double M_Neut)
{
	SetMass(M_Lepton0);
	double dML2 = M_Lepton*M_Lepton/Mass(2);
	double dMn2 = M_Neut*M_Neut/Mass(2);
	double dMN2 = MassN(2)/Mass(2);

	double s, t, cos0, cos1;
	Kinematic_3B(s, t, cos0, cos1);
	double u = 1 + dML2 + dMn2 + dMN2 - s - t;

	double M2 = M2_LeptonNeutrino(dMn2, dML2, dMN2, u);
	return dGammad5_3B(M2);
}

double PhaseSpace::Max_LeptonNeutrino(double M_Lepton0, double M_Lepton, double M_Neut)
{
	SetMass(M_Lepton0);
	double dML2 = M_Lepton*M_Lepton/Mass(2);
	double dMn2 = M_Neut*M_Neut/Mass(2);
	double dMN2 = MassN(2)/Mass(2);

	double M2 = max_LeptonNeutrino(dMn2, dML2, dMN2);
	return dGammad5_3B(M2);
}

double PhaseSpace::max_LeptonNeutrino(double x, double y, double z)	//vars is s 
{
	F_var.clear();

	F_var.push_back(x);	//0	//x
	F_var.push_back(y);	//1	//y
	F_var.push_back(z);	//2	//z

	SetFunction(&PhaseSpace::max_LeptonNeutrino_u);
	return Inte::GoldRatioSolver(this, Inte::Max);
}

double PhaseSpace::max_LeptonNeutrino_u(double u)	//vars is s 
{
	double x = F_var.at(0);
	double y = F_var.at(1);
	double z = F_var.at(2);

	double u_ = u;
	double fc = Limit(u_, y, z, x);

	return fc * M2_LeptonNeutrino(x, y, z, u_);
}

double PhaseSpace::AntiLeptonNeutrino(double M_Lepton0, double M_Lepton, double M_Neut)
{
	SetMass(M_Lepton0);
	double dML2 = M_Lepton*M_Lepton/Mass(2);
	double dMn2 = M_Neut*M_Neut/Mass(2);
	double dMN2 = MassN(2)/Mass(2);

	double s, t, cos0, cos1;
	Kinematic_3B(s, t, cos0, cos1);
	double u = 1 + dML2 + dMn2 + dMN2 - s - t;

	double M2 = M2_AntiLeptonNeutrino(dMn2, dML2, dMN2, u);
	return dGammad5_3B(M2);
}

double PhaseSpace::Max_AntiLeptonNeutrino(double M_Lepton0, double M_Lepton, double M_Neut)
{
	SetMass(M_Lepton0);
	double dML2 = M_Lepton*M_Lepton/Mass(2);
	double dMn2 = M_Neut*M_Neut/Mass(2);
	double dMN2 = MassN(2)/Mass(2);

	double M2 = max_AntiLeptonNeutrino(dMn2, dML2, dMN2);
	return dGammad5_3B(M2);
}

double PhaseSpace::max_AntiLeptonNeutrino(double x, double y, double z)	//vars is s 
{
	F_var.clear();

	F_var.push_back(x);	//0	//x
	F_var.push_back(y);	//1	//y
	F_var.push_back(z);	//2	//z

	SetFunction(&PhaseSpace::max_AntiLeptonNeutrino_s);
	return Inte::GoldRatioSolver(this, Inte::Max);
}

double PhaseSpace::max_AntiLeptonNeutrino_s(double s)	//vars is s 
{
	double x = F_var.at(0);
	double y = F_var.at(1);
	double z = F_var.at(2);

	double s_ = s;
	double fc = Limit(s_, x, y, z);

	return fc * M2_AntiLeptonNeutrino(x, y, z, s_);
}

double PhaseSpace::MesonThree(double M_Meson0, double M_Meson, double M_Lepton, double L_, double L0)	//decay constant not important
{
	SetMass(M_Meson0);
	double dMM2 = M_Meson*M_Meson/Mass(2);
	double dML2 = M_Lepton*M_Lepton/Mass(2);
	double dMN2 = MassN(2)/Mass(2);	

	double s_, t_, cos0, cos1;
	Kinematic_3B(s_, t_, cos0, cos1);

	double M2 = M2_MesonThree(s_, t_, dMM2, dML2, dMN2, L_, L0);
	return dGammad5_3B(M2);
}

double PhaseSpace::Max_MesonThree(double M_Meson0, double M_Meson, double M_Lepton, double L_, double L0)	//no angle dependence
{
	SetMass(M_Meson0);
	double dMM2 = M_Meson*M_Meson/Mass(2);
	double dML2 = M_Lepton*M_Lepton/Mass(2);
	double dMN2 = MassN(2)/Mass(2);	

	double M2 = max_MesonThree(dMM2, dML2, dMN2, L_, L0);
	return dGammad5_3B(M2);
}

double PhaseSpace::max_MesonThree(double x, double y, double z, double L_, double L0)	//no angle dependence
{
	F_var.clear();

	F_var.push_back(x);	//0	//c2
	F_var.push_back(y);	//1	//b2
	F_var.push_back(z);	//2	//a2
	//F_var.push_back(cos0);	//3	//theta
	F_var.push_back(L_);	//4	//linear dep for decay constant
	F_var.push_back(L0);	//5	//linear dep for decay constant

	std::vector<double> xMin;
	SetFunction_D(&PhaseSpace::max_MesonThree_D);
	return Inte::NelMedSolver(this, xMin, 2, Inte::Max);
}

double PhaseSpace::max_MesonThree_D(const double *p)
{
	double x  = F_var.at(0);
	double y  = F_var.at(1);
	double z  = F_var.at(2);
	double gL = F_var.at(3);
	double gR = F_var.at(4);

	double s_    = p[0];
	double t_    = p[1];

	return M2_MesonThree(s_, t_, x, y, z, gL, gR);
}

///KINEMATICS/////
/////obtaining kinematics from decay products
//
void PhaseSpace::Kinematic_2B(double &cos0)
{
	//TLorentzVector *vec0 = Daughter(0);
	TLorentzVector *vec1 = DaughterVector(1, RestFrame);

	cos0 = cos(vec1->Theta());

	//delete vec1;
}

void PhaseSpace::Kinematic_3B(double &s, double &t, double &cos0, double &cos1)
{
	//TLorentzVector *vec0 = Daughter(0);
	TLorentzVector *vec1 = DaughterVector(1, RestFrame);
	TLorentzVector *vec2 = DaughterVector(2, RestFrame);

	TLorentzVector vec_s = *Parent(RestFrame) - *vec1;
	TLorentzVector vec_t = *Parent(RestFrame) - *vec2;

	s = vec_s.M2();
	t = vec_t.M2();

	cos0 = cos(vec1->Theta());
	cos1 = cos(vec2->Theta());

	//delete vec1;
	//delete vec2;
}

//////////
//	//
//	//
//////////

TLorentzVector* PhaseSpace::DaughterVector(unsigned int i, Reference Frame)
{
	if (i < Daughter())
	{
		TLorentzVector *D_vec = Event->GetDecay(i);
		D_vec->Boost(Parent(Frame)->BoostVector());

		return D_vec;
	}
	else
		return 0;
}

Particle* PhaseSpace::Daughter(unsigned int i, Reference Frame)
{
	if (i < Daughter())
		return new Particle(vPdg.at(i), DaughterVector(i, Frame));
	else
		return 0;
}

unsigned int PhaseSpace::Daughter()
{
	return vMass.size();
}

TLorentzVector* PhaseSpace::Parent(Reference Frame)
{
	switch (Frame)
	{
		case RestFrame:
			return Rest();
		case LabFrame:
			return LabF();
		default:
			return 0;
	}
}

TLorentzVector* PhaseSpace::Rest()
{
	return P_rest;
}

TLorentzVector* PhaseSpace::LabF()
{
	return P_labf;
}

void PhaseSpace::SetLabf(TLorentzVector &Vec)	//parent labframe
{
	P_labf = &Vec;
	SetRest(P_labf->M());
}

void PhaseSpace::SetRest(double Mass)		//parent rest frame
{
	P_rest->SetPxPyPzE(0, 0, 0, Mass);
}

void PhaseSpace::Reset()
{
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

	maxMuonE  = -1.0;
        maxMuonM  = -1.0;
        maxTauEE  = -1.0;
        maxTauET  = -1.0;
        maxTauMM  = -1.0;
        maxTauMT  = -1.0;
        maxTauPI  = -1.0;
        maxTau2PI = -1.0;
        maxPionE  = -1.0;
        maxPionM  = -1.0;
        maxKaonE  = -1.0;
        maxKaonM  = -1.0;
        maxCharmE = -1.0;
        maxCharmM = -1.0;
        maxCharmT = -1.0;
        maxKaon0E = -1.0;
        maxKaon0M = -1.0;
        maxKaonCE = -1.0;
        maxKaonCM = -1.0;
}

void PhaseSpace::SetFunction(double (PhaseSpace::*FF)(const double))
{
	double (Amplitude::*Function)(const double) = 
		static_cast<double (Amplitude::*)(const double)>(FF); // ok!
	Amplitude::SetFunction(Function);
}

void PhaseSpace::SetFunction_D(double (PhaseSpace::*FF)(const double*))
{
	double (Amplitude::*Function)(const double*) = 
		static_cast<double (Amplitude::*)(const double*)>(FF); // ok!
	Amplitude::SetFunction_D(Function);
}
