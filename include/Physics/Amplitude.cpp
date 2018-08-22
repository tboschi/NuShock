/* 
 * unpolarised amplitudes for processes involving N
 * also other util procedures for derived classes
 * and definition of decay width from unpolarised amplitude
 */

#include "Amplitude.h"

Amplitude::Amplitude()	: //Decay rates calculator
	M_Neutrino(0.0),
	M_Photon(0.0),
	M_Electron(Const::fMElectron),
	M_Muon(Const::fMMuon),
	M_Pion(Const::fMPion),
	M_Pion0(Const::fMPion0),
	M_Kaon(Const::fMKaon),
	M_Kaon0(Const::fMKaon0),
	M_Eta(Const::fMEta),
        M_Rho(Const::fMRho),
        M_Rho0(Const::fMRho0),
        M_Omega(Const::fMOmega),
        M_Kaonx(Const::fMKaonx),
        M_Etai(Const::fMEtai),
        M_Phi(Const::fMPhi),
        M_Tau(Const::fMTau),
        M_Charm(Const::fMCharm),
        M_CharmS(Const::fMDs)
{
	M_Sterile_prev = -1.0;
	Channel_prev = _undefined;
}

void Amplitude::LoadMap()
{
	chMap[_undefined] = "undefined";	//0
	chMap[_ALL]    = "ALL";			//1

	//decay channels
	chMap[_nnn]    = "nnn";			//2
	chMap[_nGAMMA] = "nGAMMA";
	chMap[_nEE]    = "nEE";
	chMap[_nEM]    = "nEM";
	chMap[_nMM]    = "nMM";
	chMap[_nET]    = "nET";
	chMap[_nMT]    = "nMT";
	chMap[_nPI0]   = "nPI0";
	chMap[_EPI]    = "EPI";
	chMap[_MPI]    = "MPI";
	chMap[_TPI]    = "TPI";
	chMap[_EKA]    = "EKA";
	chMap[_MKA]    = "MKA";
	chMap[_nRHO0]  = "nRHO0";
	chMap[_ERHO]   = "ERHO";
	chMap[_MRHO]   = "MRHO";
	chMap[_EKAx]   = "EKAx";
	chMap[_MKAx]   = "MKAx";
	chMap[_nOMEGA] = "nOMEGA";
	chMap[_nETA]   = "nETA";
	chMap[_nETAi]  = "nETAi";
	chMap[_nPHI]   = "nPHI";
	chMap[_ECHARM] = "ECHARM";
	chMap[_ExpALL] = "ExpALL";		//25

	//production channels
	chMap[_MuonE]  = "MuonE";		//26
	chMap[_MuonM]  = "MuonM";
	chMap[_TauEE]  = "TauEE";
	chMap[_TauET]  = "TauET";
	chMap[_TauMM]  = "TauMM";
	chMap[_TauMT]  = "TauMT";
	chMap[_TauPI]  = "TauPI";
	chMap[_Tau2PI] = "Tau2PI";
	chMap[_PionE]  = "PionE";
	chMap[_PionM]  = "PionM";
	chMap[_KaonE]  = "KaonE";
	chMap[_KaonM]  = "KaonM";
	chMap[_CharmE] = "CharmE";
	chMap[_CharmM] = "CharmM";
	chMap[_CharmT] = "CharmT";
	chMap[_Kaon0E] = "Kaon0E";
	chMap[_Kaon0M] = "Kaon0M";
	chMap[_KaonCE] = "KaonCE";
	chMap[_KaonCM] = "KaonCM";
}

std::string Amplitude::ShowChannel(Channel Name)
{
	if (chMap.size() == 0)
		LoadMap();

	return chMap[Name];
}

Amplitude::Process Amplitude::LoadMass(Channel Name)	//return true if Decay, false if Production
{
	vMass.clear();
	vPdg.clear();
	switch(Name)
	{
		//masses in channelname order
		case _ALL:
		case _nnn:
			vMass.push_back(M_Neutrino);
			vMass.push_back(M_Neutrino);
			vMass.push_back(M_Neutrino);
			vPdg.push_back(12);
			vPdg.push_back(12);
			vPdg.push_back(12);
			return DecayRates;
		case _nGAMMA:
			vMass.push_back(M_Neutrino);
			vMass.push_back(M_Photon);
			vPdg.push_back(12);
			vPdg.push_back(22);
			return DecayRates;
		case _ExpALL:
			//neutrino lepton lepton AA
		case _nEE:
			vMass.push_back(M_Neutrino);
			vMass.push_back(M_Electron);
			vMass.push_back(M_Electron);
			vPdg.push_back(12);
			vPdg.push_back(11);
			vPdg.push_back(11);
			return DecayRates;
		case _nMM:
			vMass.push_back(M_Neutrino);
			vMass.push_back(M_Muon);
			vMass.push_back(M_Muon);
			vPdg.push_back(12);
			vPdg.push_back(13);
			vPdg.push_back(13);
			return DecayRates;
			//neutrino lepton lepton AB
		case _nEM:
			vMass.push_back(M_Neutrino);
			vMass.push_back(M_Electron);
			vMass.push_back(M_Muon);
			vPdg.push_back(12);
			vPdg.push_back(11);
			vPdg.push_back(13);
			return DecayRates;
		case _nET:
			vMass.push_back(M_Neutrino);
			vMass.push_back(M_Electron);
			vMass.push_back(M_Tau);
			vPdg.push_back(12);
			vPdg.push_back(11);
			vPdg.push_back(15);
			return DecayRates;
		case _nMT:
			vMass.push_back(M_Neutrino);
			vMass.push_back(M_Muon);
			vMass.push_back(M_Tau);
			vPdg.push_back(12);
			vPdg.push_back(13);
			vPdg.push_back(15);
			return DecayRates;
			//neutrino psuedomeson
		case _nPI0:
			vMass.push_back(M_Neutrino);
			vMass.push_back(M_Pion0);
			vPdg.push_back(12);
			vPdg.push_back(111);
			return DecayRates;
		case _nETA:
			vMass.push_back(M_Neutrino);
			vMass.push_back(M_Eta);
			vPdg.push_back(12);
			vPdg.push_back(221);
			return DecayRates;
		case _nETAi:
			vMass.push_back(M_Neutrino);
			vMass.push_back(M_Etai);
			vPdg.push_back(12);
			vPdg.push_back(331);
			return DecayRates;
			//lepton psuedomeson
		case _EPI:
			vMass.push_back(M_Electron);
			vMass.push_back(M_Pion);
			vPdg.push_back(11);
			vPdg.push_back(211);
			return DecayRates;
		case _MPI:
			vMass.push_back(M_Muon);
			vMass.push_back(M_Pion);
			vPdg.push_back(13);
			vPdg.push_back(211);
			return DecayRates;
		case _TPI:
			vMass.push_back(M_Tau);
			vMass.push_back(M_Pion);
			vPdg.push_back(211);
			vPdg.push_back(15);
			return DecayRates;
		case _EKA:
			vMass.push_back(M_Electron);
			vMass.push_back(M_Kaon);
			vPdg.push_back(11);
			vPdg.push_back(321);
			return DecayRates;
		case _MKA:
			vMass.push_back(M_Muon);
			vMass.push_back(M_Kaon);
			vPdg.push_back(13);
			vPdg.push_back(321);
			return DecayRates;
		case _ECHARM:
			vMass.push_back(M_Electron);
			vMass.push_back(M_Charm);
			vPdg.push_back(12);
			vPdg.push_back(411);
			return DecayRates;
			//neutrino vectormeson
		case _nRHO0:
			vMass.push_back(M_Neutrino);
			vMass.push_back(M_Rho0);
			vPdg.push_back(12);
			vPdg.push_back(113);
			return DecayRates;
		case _nOMEGA:
			vMass.push_back(M_Neutrino);
			vMass.push_back(M_Omega);
			vPdg.push_back(12);
			vPdg.push_back(223);
			return DecayRates;
		case _nPHI:
			vMass.push_back(M_Neutrino);
			vMass.push_back(M_Phi);
			vPdg.push_back(12);
			vPdg.push_back(333);
			return DecayRates;
			//lepton vectormeson
		case _ERHO:
			vMass.push_back(M_Electron);
			vMass.push_back(M_Rho);
			vPdg.push_back(11);
			vPdg.push_back(213);
			return DecayRates;
		case _MRHO:
			vMass.push_back(M_Muon);
			vMass.push_back(M_Rho0);
			vPdg.push_back(13);
			vPdg.push_back(213);
			return DecayRates;
		case _EKAx:
			vMass.push_back(M_Electron);
			vMass.push_back(M_Kaonx);
			vPdg.push_back(11);
			vPdg.push_back(9000321);
			return DecayRates;
		case _MKAx:
			vMass.push_back(M_Muon);
			vMass.push_back(M_Kaonx);
			vPdg.push_back(13);
			vPdg.push_back(9000321);
			return DecayRates;
		//PRODUCTION
		//Parent first
		case _MuonE:
		case _MuonM:
			vMass.push_back(M_Muon);
			vMass.push_back(M_Electron);
			vMass.push_back(M_Neutrino);
			vPdg.push_back(13);
			vPdg.push_back(11);
			vPdg.push_back(12);
			return Production;
		case _TauEE:
		case _TauET:
			vMass.push_back(M_Tau);
			vMass.push_back(M_Electron);
			vMass.push_back(M_Neutrino);
			vPdg.push_back(15);
			vPdg.push_back(11);
			vPdg.push_back(12);
			return Production;
		case _TauMM:
		case _TauMT:
			vMass.push_back(M_Tau);
			vMass.push_back(M_Muon);
			vMass.push_back(M_Neutrino);
			vPdg.push_back(15);
			vPdg.push_back(13);
			vPdg.push_back(12);
			return Production;
		case _TauPI:
			vMass.push_back(M_Tau);
			vMass.push_back(M_Pion);
			vPdg.push_back(15);
			vPdg.push_back(211);
			return Production;
		case _Tau2PI:
			vMass.push_back(M_Tau);
			vMass.push_back(M_Pion);
			vMass.push_back(M_Pion0);
			vPdg.push_back(15);
			vPdg.push_back(211);
			vPdg.push_back(111);
			return Production;
		case _PionE:
			vMass.push_back(M_Pion);
			vMass.push_back(M_Electron);
			vPdg.push_back(211);
			vPdg.push_back(11);
			return Production;
		case _PionM:
			vMass.push_back(M_Pion);
			vMass.push_back(M_Muon);
			vPdg.push_back(211);
			vPdg.push_back(13);
			return Production;
		case _KaonE:
			vMass.push_back(M_Kaon);
			vMass.push_back(M_Electron);
			vPdg.push_back(321);
			vPdg.push_back(11);
			return Production;
		case _KaonM:
			vMass.push_back(M_Kaon);
			vMass.push_back(M_Muon);
			vPdg.push_back(321);
			vPdg.push_back(13);
			return Production;
		case _CharmE:
			vMass.push_back(M_CharmS);
			vMass.push_back(M_Electron);
			vPdg.push_back(431);
			vPdg.push_back(11);
			return Production;
		case _CharmM:
			vMass.push_back(M_CharmS);
			vMass.push_back(M_Muon);
			vPdg.push_back(431);
			vPdg.push_back(13);
			return Production;
		case _CharmT:
			vMass.push_back(M_CharmS);
			vMass.push_back(M_Tau);
			vPdg.push_back(431);
			vPdg.push_back(15);
			return Production;
		case _Kaon0E:
			vMass.push_back(M_Kaon0);
			vMass.push_back(M_Pion);
			vMass.push_back(M_Electron);
			vPdg.push_back(130);
			vPdg.push_back(211);
			vPdg.push_back(11);
			return Production;
		case _Kaon0M:
			vMass.push_back(M_Kaon0);
			vMass.push_back(M_Pion);
			vMass.push_back(M_Muon);
			vPdg.push_back(130);
			vPdg.push_back(211);
			vPdg.push_back(13);
			return Production;
		case _KaonCE:
			vMass.push_back(M_Kaon);
			vMass.push_back(M_Pion0);
			vMass.push_back(M_Electron);
			vPdg.push_back(321);
			vPdg.push_back(111);
			vPdg.push_back(11);
			return Production;
		case _KaonCM:
			vMass.push_back(M_Kaon);
			vMass.push_back(M_Pion0);
			vMass.push_back(M_Muon);
			vPdg.push_back(321);
			vPdg.push_back(111);
			vPdg.push_back(13);
			return Production;
		default:
			std::cerr << ShowChannel(Name) << " : Channel not recognised" << std::endl;
			return Undefined;
	}
}

//Kinematic function
double Amplitude::Kallen(double X, double Y, double Z)
{
	return X*X + Y*Y + Z*Z - 2*(X*Y + X*Z + Y*Z);
}

double Amplitude::SqrtKallen(double X, double Y, double Z)
{
	double Kal = Kallen(X, Y, Z);
	if (Kal < 0)
		return sqrt(-Kal);
	else
		return sqrt(Kal);
}

/////////////////////
//Diff decay widths//
////////////////////

//	to be integrated over ds/m0², dt/m0², da, dcosb, dc
//	ds = -2M dE
double Amplitude::dGammad5_3B(double M2)
{
	return M2 * Mass() / (2048 * Const::fPi5);	//integration over dsigma dtau
}

//	after integrating angular depend., to be integrated over ds, dt
double Amplitude::dGammad2_3B(double M2)
{
	return 8 * Const::fPi2 * dGammad5_3B(M2);
}

//	to be integrated over da, dcosb		(m1/m0)²  (m2/m0)²
double Amplitude::dGammad2_2B(double M2, double x, double y)
{
	return M2 * SqrtKallen(1, x, y) / (64 * Const::fPi2 * Mass());
}

//	constant, after integrating angular depend.
double Amplitude::dGammad0_2B(double M2, double x, double y)
{
	return 4 * Const::fPi * dGammad2_2B(M2, x, y);
}

//integration limits for three body when is constant wrt to t
double Amplitude::Limit(double &s, double x, double y, double z)
{
	double t = 0;
	return Limit(s, t, x, y, z);
}

//integration limits for three body, where z refers to s, y referst to t, and x refers to u
double Amplitude::Limit(double &s, double &t, double x, double y, double z)
{
	double s0 = s;
	double sInf = x + y + 2*sqrt(x*y);
	double sSup = 1 + z - 2*sqrt(z);
	s = sInf + (sSup - sInf) * s;	//this is s

	double Kal = SqrtKallen(s, y, x) * SqrtKallen(1, s, z);
	double tInf = z + x + ((1 - s - z) * (s - y + x) - Kal) / (2 * s);
	double tSup = z + x + ((1 - s - z) * (s - y + x) + Kal) / (2 * s);
	t = tInf + (tSup - tInf) * t;	//this is t

	return (sSup - sInf) * (tSup - tInf);
}

//////////////////////
//Generic amplitudes//
//////////////////////

//					      angle lep    lepton    meson	
double Amplitude::M2_LeptonPseudoMeson(double cos0, double x, double y)
{
	return (GetFermion() ? 0.0 : M2_LeptonPseudoMeson(-Helicity(), cos0, x, y)) +	//added only if majorana
				     M2_LeptonPseudoMeson( Helicity(), cos0, x, y);
}

double Amplitude::M2_LeptonPseudoMeson(int Hel, double cos0, double x, double y)
{
	return Const::fGF2 * Mass(4) * 
		(pow(1 - x, 2) - y * (1 + x) - (1 - x) * Hel * SqrtKallen(1, x, y) * cos0);
}

//						angle lep    lepton    meson	
double Amplitude::M2_NeutrinoPseudoMeson(double cos0, double x, double y)
{
	return (GetFermion() ? 0.0 : M2_NeutrinoPseudoMeson(-Helicity(), cos0, x, y)) +	//added only if majorana
				     M2_NeutrinoPseudoMeson( Helicity(), cos0, x, y);
}

double Amplitude::M2_NeutrinoPseudoMeson(int Hel, double cos0, double x, double y)
{
	return M2_LeptonPseudoMeson(Hel, x, y, cos0) / 2.0;
}

double Amplitude::M2_LeptonVectorMeson(double cos0, double x, double y)
{
	return (GetFermion() ? 0.0 : M2_LeptonVectorMeson(-Helicity(), cos0, x, y)) +	//added only if majorana
				     M2_LeptonVectorMeson( Helicity(), cos0, x, y);
}

//					lepton		meson	
double Amplitude::M2_LeptonVectorMeson(int Hel, double cos0, double x, double y)	//must be divided by vector meson mass
{
	return Const::fGF2 * Mass(4) * 
		(pow(1 - x, 2) + y * (1 + x) - 2 * y*y - (1 - x - 2*y) *  Hel * SqrtKallen(1, x, y) * cos0);
}

double Amplitude::M2_NeutrinoVectorMeson(double cos0, double x, double y)
{
	return (GetFermion() ? 0.0 : M2_NeutrinoVectorMeson(-Helicity(), cos0, x, y)) +	//added only if majorana
				     M2_NeutrinoVectorMeson( Helicity(), cos0, x, y);
}

//					lepton		meson	angle
double Amplitude::M2_NeutrinoVectorMeson(int Hel, double cos0, double x, double y)
{
	return M2_LeptonVectorMeson(Hel, cos0, x, y) / 2.0;
}

double Amplitude::M2_WW(double s, double cos0, double x, double y, double z)	//gL^2 + gR^2
{
	return (GetFermion() ? 0.0 : M2_WW(-Helicity(), s, cos0, x, y, z)) +	//added only if majorana
				     M2_WW( Helicity(), s, cos0, x, y, z);
}
//			lepton energy is s = (p0-p2)² and cos0s the angle wrt z-axis
//					       neutrino, letpon,   lepton
double Amplitude::M2_WW(int Hel, double s, double cos0s, double x, double y, double z)	//gL^2 + gR^2
{
	return 16 * Const::fGF2 * Mass(4) *
	       (s - x - y) * (1 + z - s - Hel * SqrtKallen(1, z, s) * cos0s);
}

double Amplitude::M2_WZ(double u, double cos0u, double x, double y, double z)	//gL^2 + gR^2
{
	return (GetFermion() ? 0.0 : M2_WZ(-Helicity(), u, cos0u, x, y, z)) +	//added only if majorana
				     M2_WZ( Helicity(), u, cos0u, x, y, z);
}
//			lepton energiy is u = (p0-p1)² and cos0u the angle wrt z-axis
//			       neutrino, letpon,   lepton
double Amplitude::M2_WZ(int Hel, double u, double cos0u, double x, double y, double z)	//2gL*gR
{
	return 16 * Const::fGF2 * Mass(4) *
		sqrt(y * z) * (1 + x - u - Hel * SqrtKallen(1, x, u) * cos0u);
}

double Amplitude::M2_WZ(double s, double t, double cos0s, double cos0t, double x, double y, double z)	//2gL*gR
{
	return (GetFermion() ? 0.0 : M2_WZ(-Helicity(), s, t, cos0s, cos0t, x, y, z)) +	//added only if majorana
				     M2_WZ( Helicity(), s, t, cos0s, cos0t, x, y, z);
}
//			lepton energies are s = (p0-p2)², t = (p0-p3)² and cos0s,t the angles wrt z-axis
//			       neutrino, letpon,   lepton
double Amplitude::M2_WZ(int Hel, double s, double t, double cos0s, double cos0t, double x, double y, double z)	//2gL*gR
{
	double u = 1 + x + y + z - s - t;
	double cos0u = (SqrtKallen(1, y, s) * cos0s + SqrtKallen(1, z, t) * cos0t) / 
			SqrtKallen(1, x, u) ;

	return M2_WZ(Hel, u, cos0u, x, y, z);
}

//
//////////////
//production//
//////////////
//
//		This amplitude is to be used if the mixing comes from the decaying flavour
//		   production is from lepton, described in Jackson Frame!
//				           neutrino  lepton    neutrino
double Amplitude::M2_LeptonNeutrino(double u, double x, double y, double z)
{
	return (GetFermion() ? 0.0 : M2_LeptonNeutrino(-Helicity(), u, x, y, z)) +	//added only if majorana
				     M2_LeptonNeutrino( Helicity(), u, x, y, z);
}

double Amplitude::M2_LeptonNeutrino(int Hel, double u, double x, double y, double z)
{
	return 16 * Const::fGF2 * Mass(4) *
		(1 + x - u) * (u - y - z - Hel * SqrtKallen(u, y, z));
}

//	This amplitude is to be used if the mixing comes from the flavour in final state
//				     neutrino  lepton    neutrino  neutrino  angle betw. lepton and neutr
//	production is from antilepeton
double Amplitude::M2_AntiLeptonNeutrino(double s, double x, double y, double z)
{
	return (GetFermion() ? 0.0 : M2_AntiLeptonNeutrino(-Helicity(), s, x, y, z)) +	//added only if majorana
				     M2_AntiLeptonNeutrino( Helicity(), s, x, y, z);
}
double Amplitude::M2_AntiLeptonNeutrino(int Hel, double s, double x, double y, double z)
{
	return 16 * Const::fGF2 * Mass(4) * 
		(s - x - y) * (1 + z - s - Hel * SqrtKallen(1, s, z));
}

double Amplitude::M2_LeptonTwo(double x, double y)
{
	return (GetFermion() ? 0.0 : M2_LeptonTwo(-Helicity(), x, y)) +	//added only if majorana
				     M2_LeptonTwo( Helicity(), x, y);
}
//					      neutrino	meson
double Amplitude::M2_LeptonTwo(int Hel, double x, double y)	//y is the meson
{
	return Const::fGF2 * Mass(4) * 
		(pow(1 - x, 2) - y * (1 + x) - (1 - x) * Hel * SqrtKallen(1, x, y));
}

double Amplitude::M2_LeptonThree(double x, double y, double z)
{
	return (GetFermion() ? 0.0 : M2_LeptonThree(-Helicity(), x, y, z)) +	//added only if majorana
				     M2_LeptonThree( Helicity(), x, y, z);
}
//	not implemented
double Amplitude::M2_LeptonThree(int Hel, double x, double y, double z)	
{
	return Const::fGF2 * Mass(4);
}

double Amplitude::M2_MesonTwo(double x, double y)
{
	return (GetFermion() ? 0.0 : M2_MesonTwo(-Helicity(), x, y)) +	//added only if majorana
				     M2_MesonTwo( Helicity(), x, y);
}
//					   neutrino  lepton
double Amplitude::M2_MesonTwo(int Hel, double x, double y)
{
	return Const::fGF2 * Mass(4) * 
		(x + y - pow(x - y, 2) - Hel * (y - x) * SqrtKallen(1, x, y));
}

double Amplitude::M2_MesonThree(double s, double t, double x, double y, double z, double L_, double L0)
{
	return (GetFermion() ? 0.0 : M2_MesonThree(-Helicity(), s, t, x, y, z, L_, L0)) +	//added only if majorana
				     M2_MesonThree( Helicity(), s, t, x, y, z, L_, L0);
}
//Jackson frame??
double Amplitude::M2_MesonThree(int Hel, double s, double t, double x, double y, double z, double L_, double L0)
{
	double u = 1 + x + y + z - s - t;

	double F = 2 * (1 + L_ * u / x);
	double G = (1 + L_ * u / x) - (L_ - L0) * (1 + 1 / x);

	double A = (1 + y - t)*(1 + z - s - Hel * SqrtKallen(1, z, s)) -
	           (u - y - z - Hel * SqrtKallen(u, y, z));
	double B = (y + z) * (u - y - z) + 4 * y * z -
	           (y - z) * Hel * SqrtKallen(u, y, z);
	double C = (1 + y - t) * 2*z + (1 + z - s) * (2*y + Hel * SqrtKallen(u, y, z)) -
		   Hel * (u - z + y) * SqrtKallen(1, z, s);

	return Const::fGF2 * ( (F*F) * A + (G*G) * B - (F*G) * C ) / 2.0;
}

//Generic function set up for template analysis
//
void Amplitude::SetFunction(double (Amplitude::*FF)(double))
{
	fFunction = FF;
}

void Amplitude::SetFunction_D(double (Amplitude::*FF)(double*))
{
	fFunction_D = FF;
}

double Amplitude::Function(double x)
{
	return (this->*fFunction)(x);
}

double Amplitude::Function_D(double *x)
{
	return (this->*fFunction_D)(x);
}

bool Amplitude::IsChanged()
{
	bool Ret = (fabs(MassN() - M_Sterile_prev) > 1e-9) ||
		   (Helicity() != iHel_prev);

	M_Sterile_prev = MassN();
	iHel_prev = Helicity();

	//Reset decay widths if changed
	if (Ret)
		Reset();

	return Ret;
}

//Get functions
double Amplitude::Mass(int E)		//decaying particle mass
{
	return pow(M_Parent, E);
}

double Amplitude::MassN(int E)		//neutrino mass
{
	return pow(M_Sterile, E);
}

double Amplitude::Ue(int E)
{
	return pow(fUe, E);
}

double Amplitude::Um(int E)
{
	return pow(fUm, E);
}

double Amplitude::Ut(int E)
{
	return pow(fUt, E);
}

int Amplitude::Helicity()
{
	return GetParticle() ? iHel : -iHel;
}

bool Amplitude::GetFermion()
{
	return bFermion;
}

bool Amplitude::GetParticle()
{
	return bParticle;
}

//Set functions
void Amplitude::SetMass(double Mass)		//mass of decaying particle
{
	M_Parent = Mass;
}

void Amplitude::SetMassN(double Mass)		//mass of neutrino
{
	M_Sterile = Mass;
}

void Amplitude::SetUe(double Ue)
{
	fUe = Ue;
}

void Amplitude::SetUm(double Um)
{
	fUm = Um;
}

void Amplitude::SetUt(double Ut)
{
	fUt = Ut;
}

void Amplitude::SetFermion(bool Fermion)
{
	bFermion = Fermion;	//true for Dirac, false for Majorana
}

void Amplitude::SetParticle(bool Particle)
{
	bParticle = Particle;	//true for Dirac, false for Majorana
}

void Amplitude::SetHelicity(int Helicity)
{
	iHel = Helicity;	//-1 for Left, +1 for Right, 0 for unpolarised
}

void Amplitude::SetNeutrino(double Mass, double* Mixings, bool Fermion, bool Particle, int Helix)
{
	SetMassN(Mass);
	SetUe(Mixings[0]);
	SetUm(Mixings[1]);
	SetUt(Mixings[2]);
	SetFermion(Fermion);
	SetParticle(Particle);
	SetHelicity(Helix);
}
