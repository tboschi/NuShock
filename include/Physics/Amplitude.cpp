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
        M_Charm(Const::fMCharm)
{
	M_Sterile_prev = -1.0;
	Channel_prev = _undefined;
}

void Amplitude::LoadMass(Channel Name)
{
	vMass.clear();
	switch(Name)
	{
		case _ALL:
		case _nnn:
			vMass.push_back(M_Neutrino);
			vMass.push_back(M_Neutrino);
			vMass.push_back(M_Neutrino);
			break;
		case _nGAMMA:
			vMass.push_back(M_Neutrino);
			vMass.push_back(M_Photon);
			break;
		case _nEE:
		case _ExpALL:
			vMass.push_back(M_Neutrino);
			vMass.push_back(M_Electron);
			vMass.push_back(M_Electron);
			break;
		case _nEMU:
			vMass.push_back(M_Neutrino);
			vMass.push_back(M_Electron);
			vMass.push_back(M_Muon);
			break;
		case _nMUE:
			vMass.push_back(M_Neutrino);
			vMass.push_back(M_Electron);
			vMass.push_back(M_Muon);
			break;
		case _nMUMU:
			vMass.push_back(M_Neutrino);
			vMass.push_back(M_Muon);
			vMass.push_back(M_Muon);
			break;
		case _nET:
			vMass.push_back(M_Neutrino);
			vMass.push_back(M_Electron);
			vMass.push_back(M_Tau);
			break;
		case _nTE:
			vMass.push_back(M_Neutrino);
			vMass.push_back(M_Electron);
			vMass.push_back(M_Tau);
			break;
		case _nMUT:
			vMass.push_back(M_Neutrino);
			vMass.push_back(M_Muon);
			vMass.push_back(M_Tau);
			break;
		case _nTMU:
			vMass.push_back(M_Neutrino);
			vMass.push_back(M_Muon);
			vMass.push_back(M_Tau);
			break;
		case _nPI0:
			vMass.push_back(M_Neutrino);
			vMass.push_back(M_Pion0);
			break;
		case _EPI:
			vMass.push_back(M_Electron);
			vMass.push_back(M_Pion);
			break;
		case _MUPI:
			vMass.push_back(M_Muon);
			vMass.push_back(M_Pion);
			break;
		case _TPI:
			vMass.push_back(M_Pion);
			vMass.push_back(M_Tau);
			break;
		case _EKA:
			vMass.push_back(M_Electron);
			vMass.push_back(M_Kaon);
			break;
		case _MUKA:
			vMass.push_back(M_Muon);
			vMass.push_back(M_Kaon);
			break;
		case _EKAx:
			vMass.push_back(M_Electron);
			vMass.push_back(M_Kaonx);
			break;
		case _MUKAx:
			vMass.push_back(M_Muon);
			vMass.push_back(M_Kaonx);
			break;
		case _nRHO0:
			vMass.push_back(M_Neutrino);
			vMass.push_back(M_Rho0);
			break;
		case _ERHO:
			vMass.push_back(M_Electron);
			vMass.push_back(M_Rho);
			break;
		case _MURHO:
			vMass.push_back(M_Muon);
			vMass.push_back(M_Rho0);
			break;
		case _nETA:
			vMass.push_back(M_Neutrino);
			vMass.push_back(M_Eta);
			break;
		case _nETAi:
			vMass.push_back(M_Neutrino);
			vMass.push_back(M_Etai);
			break;
		case _nOMEGA:
			vMass.push_back(M_Neutrino);
			vMass.push_back(M_Omega);
			break;
		case _nPHI:
			vMass.push_back(M_Neutrino);
			vMass.push_back(M_Phi);
			break;
		case _ECHARM:
			vMass.push_back(M_Electron);
			vMass.push_back(M_Charm);
			break;
		case _MuonE:
		case _MuonM:
			vMass.push_back(M_Muon);
			vMass.push_back(M_Electron);
			vMass.push_back(M_Neutrino);
			break;
		case _TauEE:
		case _TauET:
			vMass.push_back(M_Tau);
			vMass.push_back(M_Electron);
			vMass.push_back(M_Neutrino);
			break;
		case _TauMM:
		case _TauMT:
			vMass.push_back(M_Tau);
			vMass.push_back(M_Muon);
			vMass.push_back(M_Neutrino);
			break;
		case _PionE:
			vMass.push_back(M_Pion);
			vMass.push_back(M_Electron);
			break;
		case _PionM:
			vMass.push_back(M_Pion);
			vMass.push_back(M_Muon);
			break;
		case _KaonE:
			vMass.push_back(M_Kaon);
			vMass.push_back(M_Electron);
			break;
		case _KaonM:
			vMass.push_back(M_Kaon);
			vMass.push_back(M_Muon);
			break;
		case _CharmE:
			vMass.push_back(M_CharmS);
			vMass.push_back(M_Electron);
			break;
		case _CharmM:
			vMass.push_back(M_CharmS);
			vMass.push_back(M_Muon);
			break;
		case _CharmT:
			vMass.push_back(M_CharmS);
			vMass.push_back(M_Tau);
			break;
		case _Kaon0E:
			vMass.push_back(M_Kaon0);
			vMass.push_back(M_Pion);
			vMass.push_back(M_Electron);
			break;
		case _Kaon0M:
			vMass.push_back(M_Kaon0);
			vMass.push_back(M_Pion);
			vMass.push_back(M_Muon);
			break;
		case _KaonCE:
			vMass.push_back(M_Kaon);
			vMass.push_back(M_Pion0);
			vMass.push_back(M_Electron);
			break;
		case _KaonCM:
			vMass.push_back(M_Kaon);
			vMass.push_back(M_Pion0);
			vMass.push_back(M_Muon);
			break;
		default:
			std::cerr << "Channel not recognised" << std::endl;
			break;
	}
}

//Kinematic function
double Amplitude::Kallen(double X, double Y, double Z)
{
	return X*X + Y*Y + Z*Z - 2*(X*Y + X*Z + Y*Z);
}

double Amplitude::tKallen(double X, double Y, double Z)
{
	return sqrt(Kallen(X, Y, Z));
}

/////////////////////
//Diff decay widths//
////////////////////

//	to be integrated over ds/m0², dt/m0², da, dcosb, dc
//	ds = -2M dE
double Amplitude::dGammad5_3B(double M2)
{
	//return M2 / (512 * Const::fPi3 * Const::fPi2 * GetMass());
	return M2 * GetMass() / (512 * Const::fPi3 * Const::fPi2);	//integration over dsigma dtau
}

//	after integrating angular depend., to be integrated over ds, dt
double Amplitude::dGammad2_3B(double M2)
{
	return 8 * Const::fPi2 * dGammad5_3B(M2);
}

//	to be integrated over da, dcosb		(m1/m0)²  (m2/m0)²
double Amplitude::dGammad2_2B(double M2, double x, double y)
{
	return M2 * SqrtKallen(1, x, y) / (64 * Const::fPi2 * GetMass());
}

//	constant, after integrating angular depend.
double Amplitude::dGammad0_2B(double M2, double x, double y)
{
	return 4 * Const::fPi * dGammad0_2B(M2, x, y);
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
	double sInf = x + y + 2*sqrt(x*y);
	double sSup = 1 + z - sqrt(z);
	s = xInf + (xSup - xInf) * s;	//this is s

	double Kal = SqrtKallen(s, y, x) * SqrtKallen(1, s, z);
	double tInf = z + x + ((1 - s - z) * (s - y + x) - Kal) / (2 * s);
	double tSup = z + x + ((1 - s - z) * (s - y + x) + Kal) / (2 * s);
	//double tSup = ( (2 - s)*(s + y - x) + sqrt(Kine::Kallen(s, y, x)*Kine::Kallen(1, s, z)) ) / (2*s);
	t = tInf + (tSup - tInf) * t;	//this is t

	return (sSup - sInf) * (tSup - tInf);
}

//////////////////////
//Generic amplitudes//	//check prefactors!
//////////////////////

//					lepton		meson	angle
double Amplitude::M2_LeptonPseudoMeson(double x, double y, double cos0)
{
	return Const::fGF2 * pow(GetMass(), 4) * 
		(pow(1 - x, 2) - y * (1 + x) - (1 - x) * GetHelicity() * SqrtKallen(1, x, y) * cos0);
}

//					neutrino	meson	angle
double Amplitude::M2_NeutrinoPseudoMeson(double x, double y, double cos0)
{
	return M2_LeptonPseudoMeson(x, y, cos0) / 2.0;
}

//					lepton		meson	angle
double Amplitude::M2_LeptonVectorMeson(double x, double y, double cos0)	//must be divided by vector meson mass
{
	return Const::fGF2 * pow(GetMass(), 4) * 
		(pow(1 - x, 2) + y * (1 + x) - 2*y - (1 - x - 2*y) *  GetHelicity() * SqrtKallen(1, x, y) * cos0);
}

//					lepton		meson	angle
double Amplitude::M2_NeutrinoVectorMeson(double x, double y, double cos0)
{
	M2_NeutrinoVectorMeson(x, y, cos0) / 2.0;
}

//			lepton energy is s = (p0-p2)² and cos0s the angle wrt z-axis
//			       neutrino, letpon,   lepton
double Amplitude::M2_WW(double x, double y, double z, double s, double cos0)	//gL^2 + gR^2
{
	return 4 * Const::fGF2 * pow(GetMass(), 4) *
	       (s - x - y) * (1 + z - s - GetHelicity() * SqrtKallen(1, z, s) * cos0);
}

//			lepton energies are s = (p0-p2)², t = (p0-p3)² and cos0s,t the angles wrt z-axis
//			       neutrino, letpon,   lepton
double Amplitude::M2_WZ(double x, double y, double z, double s, double t, double cos0s, double cos0t)	//2gL*gR
{
	double u = 1 + x + y + z - s - t;
	double cos0u = (SqrtKallen(1, y, s) * cos0s + SqrtKallen(1, z, t) * cos0t) / 
			SqrtKallen(1, x, u) ;

	return M2_WZ(x, y, z, u, cos0u);
}

//			lepton energiy is u = (p0-p1)² and cos0u the angle wrt z-axis
//			       neutrino, letpon,   lepton
double Amplitude::M2_WZ(double x, double y, double z, double u, double cos0u)	//2gL*gR
{
	return 4 * Const::fGF2 * pow(GetMass(), 4) *
		sqrt(y * z) * (1 + x - u - GetHelicity() * SqrtKallen(1, x, u) * cos0u);
}

//
//////////////
//production//
//////////////
//
//		This amplitude is to be used if the mixing comes from the decaying flavour
//					   neutrino lepton  neutrino  neutrino  lepton     angle betw. lepton and neutr
double Amplitude::M2_LeptonNeutrino(double x, double y, double z, double s, double t, double cos0)
{
	double u = 1 + x + y + z - s - t;

	return 16 * Const::fGF2 * pow(GetMass(), 4) *
		(1 + x - u) * (u + y + z -
		GetHelicity() * (u + y + z - (1 + y - t + SqrtKallen(1, y, t) * cos0) * (1 + z - s - SqrtKallen(1, z, s)) / 2) );
}

//		This amplitude is to be used if the mixing comes from the flavour in final state
//						     neutrino  lepton    neutrino  neutrino  angle betw. lepton and neutr
double ProductionRates::M2_LeptonAntineutrino(double x, double y, double z, double s, double cos0)
{
	return 16 * Const::fGF2 * pow(GetMass(), 4) * 
		(s - x - y) * (1 + z - s - GetHelicity() * SqrtKallen(1, s, z));
}

//					      neutrino	meson
double ProductionRates::M2_LeptonMeson(double x, double y)	//y is the meson
{
	return Const::fGF2 * pow(GetMass(), 4) * 
		(pow(1 - x, 2) - y * (1 + x) - (1 - x) * GetHelicity() SqrtKallen(1, x, y));
}

//					   neutrino  lepton
double ProductionRates::M2_MesonTwo(double x, double y)
{
	return Const::fGF2 * pow(GetMass(), 4) * 
		(x + y - pow(x - y, 2) - GetHelicity() * (y - x) * SqrtKallen(1, x, y));
}

double Production::M2_MesonThree(double s, double t, double x, double y, double cos0, double L_, double L0)
{
	double F = 2 + 2 * (1 + x + y + z - s - t) * L_ / x;
	double G = 1 + (2 + y + z - s - t) * L_ / x - (1 - x) * L0 / X;

	double A = (1 + z - s) * (1 + y - t) - (1 + x - s - t) -
		GetHelicity() * ( SqrtKallen(1, z, s) * (1 + y - t) + SqrtKallen(1, y, t) * (1 + z - s) * cos0 ) / 2.0;
	double B = (y + z) * (1 + x - s - t) / + 4 * y * z -
		GetHelicity() * (y - z) * ( SqrtKallen(1, z, s) * (1 + y - t) - SqrtKallen(1, y, t) * (1 + z - s) * cos0 ) / 2.0;
	double C = y * (1 + z - s) + z * (1 + y - t) -
	       	GetHelicity() * (y * SqrtKallen(1, z, s) + z * SqrtKallen(1, y, t) * cos0);

	return Const::fGF2 * ( (F*F) * A + (G*G) * B - 2*(F*G) * C ) / 2.0;
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
	return (*fFunction)(x);
}

double Amplitude::Function_D(double *x)
{
	return (*fFunction_D)(x);
}

/*
std::vector<std::string> Amplitude::ListChannels()
{
	std::map<std::string, ChannelName>::iterator it;
	std::vector<std::string> vList;
	for (it = mapChannel.begin(); it != mapChannel.end(); ++it)
		vList.push_back(it->first);
	return vList;
}
*/
bool Amplitude::IsChanged()
{
	bool Ret = (fabs(GetMass() - M_Sterile_prev) > 1e-9);

	M_Sterile_prev = GetMass();

	//Reset decay widths if changed
	if (Ret)
		Reset();

	return Ret;
}

//Get functions
double Amplitude::GetMass()
{
	return M_Sterile;
}

double Amplitude::GetUe()
{
	return Ue;
}

double Amplitude::GetUm()
{
	return Um;
}

double Amplitude::GetUt()
{
	return Ut;
}

bool Amplitude::GetFermion()
{
	return bFermion;
}

bool Amplitude::GetParticle()
{
	return bParticle;
}

int Amplitude::GetHelicity()
{
	return iHel;
}

//Set functions
void Amplitude::SetMass(double Mass)
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

void Amplitude::SetNeutrino(double Mass, double* Mixings, bool Fermion, bool Particle, bool Helicity)
{
	SetMass(Mass);
	SetUe(Mixings[0]);
	SetUm(Mixings[1]);
	SetUt(Mixings[2]);
	SetFermion(Fermion);
	SetParticle(Particle);
	SetHelicity(Helicity);
}
