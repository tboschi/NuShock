#include "PhaseSpace.h"

PhaseSpace::PhaseSpace() :
	Event(new TGenPhaseSpace),
	GenMT(new TRandom3(0))
	//P_labf(new TLorentzVector()),
	//P_rest(new TLorentzVector())
{
	Reset();
}

PhaseSpace::~PhaseSpace()
{
	delete Event;
	delete GenMT;
	//delete P_labf;
	//delete P_rest;
}

bool PhaseSpace::SetDecay(Channel Name)
{
	double MassArray[10];
	switch (LoadMass(Name))
	{
		case DecayRates:
			MassArray[0] = vMass.at(0);
			break;
		case Production:
			MassArray[0] = MassN();
			break;
		case Undefined:
			return false;
	}

	for (int i = 1; i < Daughters(); ++i)
		MassArray[i] = vMass.at(i);

	TLorentzVector vec = Parent(restFrame);
	return Event->SetDecay(vec, Daughters(), MassArray);
}

bool PhaseSpace::Generate(Channel Name, double &val)
{
	int cc = 0;
	if (SetDecay(Name))
	{
		do
		{
			++cc;
			Event->Generate();
			val = Ratio(Name);
		}
		while (val >= 0 && GenMT->Rndm() > val && cc < 1e4);

		if (cc > 1e3)
			return false;
		else
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
		case _nMM:
			Ret = nMM_ratio();
			break;
		case _nET:
			Ret = nET_ratio();
			break;
		case _nMT:
			Ret = nMT_ratio();
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
			Ret = TauPI_ratio();
			break;
		case _Tau2PI:				//1
			//Ret = Tau2PI_ratio();
			break;
		case _PionE:				//1
			Ret = PionE_ratio();
			break;
		case _PionM:				//1
			Ret = PionM_ratio();
			break;
		case _KaonE:				//1
			Ret = KaonE_ratio();
			break;
		case _KaonM:				//1
			Ret = KaonM_ratio();
			break;
		case _CharmE:				//1
			Ret = CharmE_ratio();
			break;
		case _CharmM:				//1
			Ret = CharmE_ratio();
			break;
		case _CharmT:				//1
			Ret = CharmT_ratio();
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
	return NeutrinoLeptonAA_ratio(maxnEE_e, maxnEE_o, &PhaseSpace::Ue);
}

double PhaseSpace::nMM_ratio()
{
	return NeutrinoLeptonAA_ratio(maxnMM_m, maxnMM_o, &PhaseSpace::Um);
}

double PhaseSpace::NeutrinoLeptonAA_ratio(double &maxName_u, double &maxName_o, double (PhaseSpace::*Uu)(int))
{
	if (maxName_u < 0 || maxName_o < 0 || IsChanged())
		Max_NeutrinoLeptonAA(maxName_u, maxName_o, vMass.at(0), vMass.at(1));

	double fName_u, fName_o;
	NeutrinoLeptonAA(fName_u, fName_o, vMass.at(0), vMass.at(1));

	if (fName_u < 0)
		fName_u = 0;
	if (fName_o < 0)
		fName_o = 0;

	return (maxName_u > 1e-25 || maxName_o > 1e-25) ? 
		( fName_u   * (this->*Uu)(2) + fName_o   * (Ue(2) + Um(2) + Ut(2) - (this->*Uu)(2)) ) / 
		( maxName_u * (this->*Uu)(2) + maxName_o * (Ue(2) + Um(2) + Ut(2) - (this->*Uu)(2)) ) : -1.0;
}

//Neutrino LeptonLepton AB
double PhaseSpace::nEM_ratio()
{
	return NeutrinoLeptonAB_ratio(maxnEM_e, maxnEM_m, &PhaseSpace::Ue, &PhaseSpace::Um);
}

double PhaseSpace::nET_ratio()
{
	return NeutrinoLeptonAB_ratio(maxnET_e, maxnET_t, &PhaseSpace::Ue, &PhaseSpace::Ut);
}

double PhaseSpace::nMT_ratio()
{
	return NeutrinoLeptonAB_ratio(maxnMT_m, maxnMT_t, &PhaseSpace::Um, &PhaseSpace::Ut);
}

double PhaseSpace::NeutrinoLeptonAB_ratio(double &maxName_a, double &maxName_b, 
					  double (PhaseSpace::*Ua)(int), double (PhaseSpace::*Ub)(int))
{
	if (maxName_a < 0 || maxName_b < 0 || IsChanged())
		Max_NeutrinoLeptonAB(maxName_a, maxName_b, vMass.at(0), vMass.at(1), vMass.at(2));

	double fName_a, fName_b;
	NeutrinoLeptonAB(fName_a, fName_b, vMass.at(0), vMass.at(1), vMass.at(2));

	if (fName_a < 0)
		fName_a = 0;
	if (fName_b < 0)
		fName_b = 0;

	return (maxName_a > 1e-25 || maxName_b > 1e-25) ? 
		( fName_a   * (this->*Ua)(2) + fName_b   * (this->*Ub)(2) ) / 
		( maxName_a * (this->*Ua)(2) + maxName_b * (this->*Ub)(2) ) : -1.0;
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

	return (maxName > 1e-25) ? 
		NeutrinoPseudoMeson(vMass.at(0), vMass.at(1)) / maxName : -1.0;
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

	return (maxName > 1e-25) ? 
		LeptonPseudoMeson(vMass.at(0), vMass.at(1)) / maxName : -1.0;
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

	return (maxName > 1e-25) ? 
		NeutrinoVectorMeson(vMass.at(0), vMass.at(1)) / maxName : -1.0;
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

	return (maxName > 1e-25) ?
		LeptonVectorMeson(vMass.at(0), vMass.at(1)) / maxName : -1.0;
}

////PRODUCTION
//
//pure leptonic decays
//
double PhaseSpace::MuonE_ratio()
{
	return AntiLeptonNeutrino_ratio(maxMuonE);
}

double PhaseSpace::MuonM_ratio()
{
	return LeptonNeutrino_ratio(maxMuonM);
}

double PhaseSpace::TauEE_ratio()
{
	return AntiLeptonNeutrino_ratio(maxTauEE);
}

double PhaseSpace::TauET_ratio()
{
	return LeptonNeutrino_ratio(maxTauET);
}

double PhaseSpace::TauMM_ratio()
{
	return AntiLeptonNeutrino_ratio(maxTauMM);
}

double PhaseSpace::TauMT_ratio()
{
	return LeptonNeutrino_ratio(maxTauMT);
}

double PhaseSpace::AntiLeptonNeutrino_ratio(double &maxName)
{
	if (maxName < 0 || IsChanged())
		maxName = Max_AntiLeptonNeutrino(vMass.at(0), vMass.at(1), vMass.at(2));

	return (maxName > 1e-25) ?
		 AntiLeptonNeutrino(vMass.at(0), vMass.at(1), vMass.at(2)) / maxName : -1.0;
}

double PhaseSpace::LeptonNeutrino_ratio(double &maxName)
{
	if (maxName < 0 || IsChanged())
		maxName = Max_LeptonNeutrino(vMass.at(0), vMass.at(1), vMass.at(2));

	return (maxName > 1e-25) ?
		 LeptonNeutrino(vMass.at(0), vMass.at(1), vMass.at(2)) / maxName : -1.0;
}

//lepton decay into meson stuff
//
double PhaseSpace::TauPI_ratio()
{
	return LeptonMeson_ratio(maxTauPI);
}

/*
double PhaseSpace::Tau2PI_ratio()
{
	return LeptonMeson_ratio(maxTauPI);
}
*/

double PhaseSpace::LeptonMeson_ratio(double &maxName)
{
	if (maxName < 0 || IsChanged())
		maxName = Max_LeptonMeson(vMass.at(0), vMass.at(1));

	return (maxName > 1e-25) ?
		 LeptonMeson(vMass.at(0), vMass.at(1)) / maxName : -1.0;
}

//meson two body decay
//
double PhaseSpace::PionE_ratio()
{
	return MesonTwo_ratio(maxPionE);
}

double PhaseSpace::PionM_ratio()
{
	return MesonTwo_ratio(maxPionM);
}

double PhaseSpace::KaonE_ratio()
{
	return MesonTwo_ratio(maxKaonE);
}

double PhaseSpace::KaonM_ratio()
{
	return MesonTwo_ratio(maxKaonM);
}

double PhaseSpace::CharmE_ratio()
{
	return MesonTwo_ratio(maxCharmE);
}

double PhaseSpace::CharmM_ratio()
{
	return MesonTwo_ratio(maxCharmM);
}

double PhaseSpace::CharmT_ratio()
{
	return MesonTwo_ratio(maxCharmT);
}

double PhaseSpace::MesonTwo_ratio(double &maxName)
{
	if (maxName < 0 || IsChanged())
		maxName = Max_MesonTwo(vMass.at(0), vMass.at(1));

	return (maxName > 1e-25) ?
		MesonTwo(vMass.at(0), vMass.at(1)) / maxName : -1.0;
}

//three body decays of meson
//
double PhaseSpace::Kaon0E_ratio()
{
	return MesonThree_ratio(maxKaon0E, Const::KCL_, Const::KCL0);
}

double PhaseSpace::Kaon0M_ratio()
{
	return MesonThree_ratio(maxKaon0M, Const::KCL_, Const::KCL0);
}

double PhaseSpace::KaonCE_ratio()
{
	return MesonThree_ratio(maxKaonCE, Const::KCL_, Const::KCL0);
}

double PhaseSpace::KaonCM_ratio()
{
	return MesonThree_ratio(maxKaonCM, Const::KCL_, Const::KCL0);
}

double PhaseSpace::MesonThree_ratio(double &maxName, double L_, double L0)
{
	if (maxName < 0 || IsChanged())
		maxName = Max_MesonThree(vMass.at(0), vMass.at(1), vMass.at(2), L_, L0);

	return (maxName > 1e-25) ?
		MesonThree(vMass.at(0), vMass.at(1), vMass.at(2), L_, L0) / maxName : -1.0;
}

///////////////
//DECAY MODES//
///////////////
//
double PhaseSpace::NeutrinoLeptonAA(double &d_Uu, double &d_Uo, double M_Neut, double M_Lepton)
{
	SetMass(MassN());
	double dMN2 = M_Neut*M_Neut/Mass(2);
	double dML2 = M_Lepton*M_Lepton/Mass(2);

	double s_, t_, u_, coss_, cost_, cosu_;
	Kinematic_3B(s_, t_, u_, coss_, cost_, cosu_);

	double gL_CC = 2 - 4*GetParticle();	//times U(lepton flavour)
	double gR_CC = 0;
	double gL_NC = -1 + 2*Const::sin2W;	//times U(neutrino flavour)
	double gR_NC = 2*Const::sin2W;

	d_Uu = NeutrinoLeptonLepton(s_, u_, coss_, cosu_, dMN2, dML2, dML2, gL_CC, gR_CC);
	d_Uo = NeutrinoLeptonLepton(s_, u_, coss_, cosu_, dMN2, dML2, dML2, gL_NC, gR_NC);

	return 0.0;
}

double PhaseSpace::Max_NeutrinoLeptonAA(double &max_Uu, double &max_Uo, double M_Neut, double M_Lepton)
{
	SetMass(MassN());
	double dMN2 = M_Neut*M_Neut/Mass(2);
	double dML2 = M_Lepton*M_Lepton/Mass(2);

	double gL_CC = -0.5 + Const::sin2W + (2*GetParticle() - 1);	//times U(lepton flavour)
	double gR_CC = Const::sin2W;

	double gL_NC = -0.5 + Const::sin2W;	//times U(neutrino flavour)
	double gR_NC = Const::sin2W;

	max_Uu = max_NeutrinoLeptonLepton(dMN2, dML2, dML2, gL_CC, gR_CC);
	max_Uo = max_NeutrinoLeptonLepton(dMN2, dML2, dML2, gL_NC, gR_NC);

	return 0.0;
}

double PhaseSpace::NeutrinoLeptonAB(double &d_Ua, double &d_Ub, double M_Neut, double M_LeptonA, double M_LeptonB)
{
	SetMass(MassN());
	double dMN2 = M_Neut*M_Neut/Mass(2);
	double dMA2 = M_LeptonA*M_LeptonB/Mass(2);
	double dMB2 = M_LeptonB*M_LeptonB/Mass(2);

	double s_, t_, u_, coss_, cost_, cosu_;
	Kinematic_3B(s_, t_, u_, coss_, cost_, cosu_);

	double gL = 1.0;	//times U(lepton flavour)
	double gR = 0.0;

	d_Ua = NeutrinoLeptonLepton(s_, u_, coss_, cosu_, dMN2, dMA2, dMB2, gL, gR);
	d_Ub = NeutrinoLeptonLepton(u_, s_, cosu_, coss_, dMN2, dMB2, dMA2, gL, gR);

	return 0.0;
}

double PhaseSpace::Max_NeutrinoLeptonAB(double &max_Ua, double &max_Ub, double M_Neut, double M_LeptonA, double M_LeptonB)
{
	SetMass(MassN());
	double dMN2 = M_Neut*M_Neut/Mass(2);
	double dMA2 = M_LeptonA*M_LeptonB/Mass(2);
	double dMB2 = M_LeptonB*M_LeptonB/Mass(2);

	double gL = 1.0;	//times U(lepton flavour)
	double gR = 0.0;

	max_Ua = max_NeutrinoLeptonLepton(dMN2, dMA2, dMB2, gL, gR);
	max_Ub = max_NeutrinoLeptonLepton(dMN2, dMB2, dMA2, gL, gR);
}

double PhaseSpace::max_NeutrinoLeptonLepton(double x, double y, double z, double gL, double gR)
{
	F_var.clear();

	F_var.push_back(x);	//0
	F_var.push_back(y);	//1
	F_var.push_back(z);	//2
	F_var.push_back(gL);	//3
	F_var.push_back(gR);	//4

	SetFunction_D(&PhaseSpace::max_NeutrinoLeptonLepton_D);
	return NelMedSolver(this, 4, -1);	//invert for max
}

double PhaseSpace::max_NeutrinoLeptonLepton_D(double *p)
{
	double x  = F_var.at(0);
	double y  = F_var.at(1);
	double z  = F_var.at(2);
	double gL = F_var.at(3);
	double gR = F_var.at(4);

	double s_    = p[0];
	double u_    = p[1];
	double cos0_ = p[2];
	double cos1_ = p[3];

	if (fabs(cos0_) <= 1 && fabs(cos1_) <= 1)
		return NeutrinoLeptonLepton(s_, u_, cos0_, cos1_, x, y, z, gL, gR);
	else
		return 0.0;
}

double PhaseSpace::NeutrinoLeptonLepton(double s, double u, double cos0, double cos1, double x, double y, double z, double gL, double gR)
{
	double M2 = (gL * gL + gR * gR) * M2_WW(s, cos0, x, y, z) +
		    (2 * gL * gR) * M2_WZ(u, cos1, x, y, z);

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

	double cos0_;
	Kinematic_2B(cos0_);

	M2_Function = &Amplitude::M2_NeutrinoPseudoMeson;
	return ToPseudoMeson(cos0_, dMN2, dMM2);
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

	double cos0_;
	Kinematic_2B(cos0_);

	M2_Function = &Amplitude::M2_LeptonPseudoMeson;
	return ToPseudoMeson(cos0_, dML2, dMM2);
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
	return GoldRatioSolver(this, -1);	//-1 to invert function
}

double PhaseSpace::max_ToPseudoMeson_cos0(double cos0)
{
	double x = F_var.at(0);
	double y = F_var.at(1);
	double cos0_ = -1 + 2 * cos0;

	return ToPseudoMeson(cos0_, x, y);
}

double PhaseSpace::ToPseudoMeson(double cos0, double x, double y)
{
	double M2 = (this->*M2_Function)(cos0, x, y);	//either M2_NeutrinoPseudoMeson 
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
	return ToVectorMeson(cos0, dMN2, dMM2);
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
	return ToVectorMeson(cos0, dMN2, dMM2);
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
	return GoldRatioSolver(this, -1);	//-1 to invert function
}

double PhaseSpace::max_ToVectorMeson_cos0(double cos0)
{
	double x = F_var.at(0);
	double y = F_var.at(1);
	double cos0_ = -1 + 2 * cos0;

	return ToVectorMeson(cos0_, x, y);
}

double PhaseSpace::ToVectorMeson(double cos0, double x, double y)
{
	double M2 = (this->*M2_Function)(cos0, x, y);
	return dGammad2_2B(M2, x, y); 
}

//////////////
//PRODUCTION//
//////////////
//
//	p1,u = N, p2,t = L, p3,s = n
double PhaseSpace::LeptonNeutrino(double M_Lepton0, double M_Lepton, double M_Neut)
{
	SetMass(M_Lepton0);
	double dML2 = M_Lepton*M_Lepton/Mass(2);
	double dMn2 = M_Neut*M_Neut/Mass(2);
	double dMN2 = MassN(2)/Mass(2);

	double s_, u_, t_, coss_, cost_, cosu_;
	Kinematic_3B(s_, t_, u_, coss_, cost_, cosu_);

	double M2 = M2_LeptonNeutrino(s_, dMn2, dML2, dMN2);
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

double PhaseSpace::max_LeptonNeutrino(double x, double y, double z)
{
	F_var.clear();

	F_var.push_back(x);	//0	//x
	F_var.push_back(y);	//1	//y
	F_var.push_back(z);	//2	//z

	SetFunction(&PhaseSpace::max_LeptonNeutrino_u);
	return GoldRatioSolver(this, -1);	//-1 to invert function
}

double PhaseSpace::max_LeptonNeutrino_u(double u)	//vars is s 
{
	double x = F_var.at(0);
	double y = F_var.at(1);
	double z = F_var.at(2);

	double u_ = u;
	double fc = Limit(u_, y, z, x);

	return M2_LeptonNeutrino(u_, x, y, z);
}

double PhaseSpace::AntiLeptonNeutrino(double M_Lepton0, double M_Lepton, double M_Neut)
{
	SetMass(M_Lepton0);
	double dML2 = M_Lepton*M_Lepton/Mass(2);
	double dMn2 = M_Neut*M_Neut/Mass(2);
	double dMN2 = MassN(2)/Mass(2);

	double s_, u_, t_, coss_, cost_, cosu_;
	Kinematic_3B(s_, u_, t_, coss_, cost_, cosu_);

	double M2 = M2_AntiLeptonNeutrino(s_, dMn2, dML2, dMN2);
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
	return GoldRatioSolver(this, -1);	//-1 to invert function
}

double PhaseSpace::max_AntiLeptonNeutrino_s(double s)	//vars is s 
{
	double x = F_var.at(0);
	double y = F_var.at(1);
	double z = F_var.at(2);

	double s_ = s;
	double fc = Limit(s_, x, y, z);

	return M2_AntiLeptonNeutrino(s_, x, y, z);
}

double PhaseSpace::LeptonMeson(double M_Lepton, double M_Meson)
{
	SetMass(M_Lepton);
	double dMN2 = MassN(2)/Mass(2);
	double dMM2 = M_Meson*M_Meson/Mass(2);

	double M2 = M2_LeptonTwo(dMN2, dMM2);
	return dGammad2_2B(M2, dMN2, dMM2);
}

double PhaseSpace::Max_LeptonMeson(double M_Lepton, double M_Meson)
{
	return LeptonMeson(M_Lepton, M_Meson);
}

double PhaseSpace::MesonTwo(double M_Meson, double M_Lepton)
{
	SetMass(M_Meson);
	double dMN2 = MassN(2)/Mass(2);
	double dML2 = M_Lepton*M_Lepton/Mass(2);

	double M2 = M2_MesonTwo(dMN2, dML2);
	return dGammad2_2B(M2, dMN2, dML2);
}

double PhaseSpace::Max_MesonTwo(double M_Meson, double M_Lepton)
{
	return MesonTwo(M_Meson, M_Lepton);
}

double PhaseSpace::MesonThree(double M_Meson0, double M_Meson, double M_Lepton, double L_, double L0)	//decay constant not important
{
	SetMass(M_Meson0);
	double dMM2 = M_Meson*M_Meson/Mass(2);
	double dML2 = M_Lepton*M_Lepton/Mass(2);
	double dMN2 = MassN(2)/Mass(2);	

	double s_, u_, t_, coss_, cost_, cosu_;
	Kinematic_3B(s_, u_, t_, coss_, cost_, cosu_);

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
	F_var.push_back(L_);	//4	//linear dep for decay constant
	F_var.push_back(L0);	//5	//linear dep for decay constant

	SetFunction_D(&PhaseSpace::max_MesonThree_D);
	return NelMedSolver(this, 2, -1);
}

double PhaseSpace::max_MesonThree_D(double *p)
{
	double x  = F_var.at(0);
	double y  = F_var.at(1);
	double z  = F_var.at(2);
	double gL = F_var.at(3);
	double gR = F_var.at(4);

	double s_ = p[0];
	double t_ = p[1];

	return M2_MesonThree(s_, t_, x, y, z, gL, gR);
}

///KINEMATICS/////
/////obtaining kinematics from decay products
//
void PhaseSpace::Kinematic_2B(double &cos0)
{
	if (Daughters() <= 2)
	{
		TLorentzVector vec1 = DaughterVector(1, restFrame);

		cos0 = cos(vec1.Theta());
	}
	else
		std::cout << "Kinematic_2B error" << std::endl;

	//delete vec1;
}

//the mandelstam variables refer to a specific vector
//s -> p3
//t -> p2
//u -> p1
void PhaseSpace::Kinematic_3B(double &s, double &t, double &u, double &cos0s, double &cos0t, double &cos0u)
{
	if (Daughters() <= 3)
	{
		TLorentzVector vec1 = DaughterVector(0, restFrame);
		TLorentzVector vec2 = DaughterVector(1, restFrame);
		TLorentzVector vec3 = DaughterVector(2, restFrame);

		TLorentzVector vec_u = Parent(restFrame) - vec1;
		TLorentzVector vec_t = Parent(restFrame) - vec2;
		TLorentzVector vec_s = Parent(restFrame) - vec3;

		u = vec_u.M2()/Mass(2);
		t = vec_t.M2()/Mass(2);
		s = vec_s.M2()/Mass(2);

		cos0u = cos(vec1.Theta());
		cos0t = cos(vec2.Theta());
		cos0s = cos(vec3.Theta());
	}
	else
		std::cout << "Kinematic_3B Error" << std::endl;
}

//////////
//	//
//	//
//////////

TLorentzVector PhaseSpace::DaughterVector(int i, Reference frame)
{
	TLorentzVector D_vec = *Event->GetDecay(i);
	D_vec.Boost(Parent(frame).BoostVector());

	return D_vec;
}

Particle PhaseSpace::Daughter(int i, Reference frame)
{
	TLorentzVector vec;
	return Particle(vPdg.at(i), DaughterVector(i, frame));
}

int PhaseSpace::Daughters()
{
	return vMass.size();
}

TLorentzVector PhaseSpace::Parent(Reference frame)
{
	switch (frame)
	{
		case restFrame:
			return RestFrame();
		case labFrame:
			return LabFrame();
	}
}

TLorentzVector PhaseSpace::RestFrame()
{
	return pRestFrame;
}

TLorentzVector PhaseSpace::LabFrame()
{
	return pLabFrame;
}

void PhaseSpace::SetLabFrame(TLorentzVector vec)	//parent labframe
{
	pLabFrame = vec;
	SetRestFrame(vec.M());
}

void PhaseSpace::SetRestFrame(double Mass)		//parent rest frame
{
	pRestFrame.SetPxPyPzE(0, 0, 0, Mass);
}

void PhaseSpace::Reset()
{
	maxnnn    = -1.0;
	maxnGAMMA = -1.0;
	maxnEE_e  = -1.0;
	maxnEE_o  = -1.0;
	maxnEM_e  = -1.0;
	maxnEM_m  = -1.0;
	maxnMM_m  = -1.0;
	maxnMM_o  = -1.0;
	maxnET_e  = -1.0;
	maxnET_t  = -1.0;
	maxnMT_m  = -1.0;
	maxnMT_t  = -1.0;
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

void PhaseSpace::SetFunction(double (PhaseSpace::*FF)(double))
{
	double (Amplitude::*Function)(double) = 
		static_cast<double (Amplitude::*)(double)>(FF); // ok!
	Amplitude::SetFunction(Function);
}

void PhaseSpace::SetFunction_D(double (PhaseSpace::*FF)(double*))
{
	double (Amplitude::*Function)(double*) = 
		static_cast<double (Amplitude::*)(double*)>(FF); // ok!
	Amplitude::SetFunction_D(Function);
}
