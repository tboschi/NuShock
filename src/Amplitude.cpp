/* 
 * unpolarised amplitudes for processes involving N
 * also other util procedures for derived classes
 * and definition of decay width from unpolarised amplitude
 */

#include "Amplitude.h"

Amplitude::Amplitude(const Neutrino &N) :
	_channel(Channel::_undefined),
	_N(N),
	_m_parent(0)
{
}

void Amplitude::LoadMass(const Channel::Name &chan)
{
	switch(chan)
	{
		//masses in channelname order
		case _ALL:
		case _nnn:
			_masspdg = {{Const::MNeutrino, 12},
				   {Const::MNeutrino, 12},
				   {Const::MNeutrino, 12}};
			break;
		case _nGAMMA:
			_masspdg = {{Const::MNeutrino, 12},
				   {Const::MPhoton, 22}};
			break;
		case _ExpALL:
			//neutrino lepton lepton AA
		case _nEE:
			_masspdg = {{Const::MNeutrino, 12}
				   {Const::MElectron, -11}
				   {Const::MElectron, 11}};
			break;
		case _nMM:
			_masspdg = {{Const::MNeutrino, 12},
				   {Const::MMuon, -13},
				   {Const::MMuon, 13}};
			break;
			//neutrino lepton lepton AB
		case _nEM:
			_masspdg = {{Const::MNeutrino, 12},
				   {Const::MElectron, 11},
				   {Const::MMuon, 13}};
			break;
		case _nET:
			_masspdg = {{Const::MNeutrino, 12},
				   {Const::MElectron, 11},
				   {Const::MTau, 15}};
			break;
		case _nMT:
			_masspdg = {{Const::MNeutrino, 12},
				   {Const::MMuon, 13},
				   {Const::MTau, 15}};
			break;
			//neutrino psuedomeson
		case _nPI0:
			_masspdg = {{Const::MNeutrino, 12},
				   {Const::MPion0, 111}};
			break;
		case _nETA:
			_masspdg = {{Const::MNeutrino, 12},
				   {Const::MEta, 221}};
			break;
		case _nETAi:
			_masspdg = {{Const::MNeutrino, 12},
				   {Const::MEtai, 331}};
			break;
			//lepton psuedomeson
		case _EPI:
			_masspdg = {{Const::MElectron, 11}
				   {Const::MPion, 211}};
			break;
		case _MPI:
			_masspdg = {{Const::MMuon, 13}
				   {Const::MPion, 211}};
			break;
		case _TPI:
			_masspdg = {{Const::MTau, 15},
				   {Const::MPion, 211}};
			break;
		case _EKA:
			_masspdg = {{Const::MElectron, 11},
				   {Const::MKaon, 321}};
			break;
		case _MKA:
			_masspdg = {{Const::MMuon, 13},
				   {Const::MKaon, 321}};
			break;
		case _ECHARM:
			_masspdg = {{Const::MElectron, 12},
				   {Const::MD, 411}};
			_pdgs = {12, 411};
			break;
			//neutrino vectormeson
		case _nRHO0:
			_masspdg = {{Const::MNeutrino, 12},
				   {Const::MRho0, 113}};
			break;
		case _nOMEGA:
			_masspdg = {{Const::MNeutrino, 12},
				   {Const::MOmega, 223}};
			break;
		case _nPHI:
			_masspdg = {{Const::MNeutrino, 12},
				   {Const::MPhi, 333}};
			break;
			//lepton vectormeson
		case _ERHO:
			_masspdg = {{Const::MElectron, 11},
				   {Const::MRho, 213}};
			break;
		case _MRHO:
			_masspdg = {{Const::MMuon, 13},
				   {Const::MRho0, 213}};
			break;
		case _EKAx:
			_masspdg = {{Const::MElectron, 11},
				   {Const::MKaonx, 9000321}};
			break;
		case _MKAx:
			_masspdg = {{Const::MMuon, 13},
				   {Const::MKaonx, 9000321}}
			break;
		//PRODUCTION
		//Parent first
		case _MuonE:
		case _MuonM:
			_masspdg = {{Const::MMuon, 13},
				   {Const::MElectron, 11},
				   {Const::MNeutrino, 12}};
			break;
		case _TauEE:
		case _TauET:
			_masspdg = {{Const::MTau, 15},
				   {Const::MElectron, 11},
				   {Const::MNeutrino, 12}};
			break;
		case _TauMM:
		case _TauMT:
			_masspdg = {{Const::MTau, 15},
				   {Const::MMuon, 13},
				   {Const::MNeutrino, 12}};
			break;
		case _TauPI:
			_masspdg = {{Const::MTau, 15},
				   {Const::MPion, 211}};
			break;
		case _Tau2PI:
			_masspdg = {{Const::MTau, 15},
				   {Const::MPion, 211},
				   {Const::MPion0, 111}};
			break;
		case _PionE:
			_masspdg = {{Const::MPion, 211},
				   {Const::MElectron, 11}};
			break;
		case _PionM:
			_masspdg = {{Const::MPion, 211},
				   {Const::MMuon, 13}};
			break;
		case _KaonE:
			_masspdg = {{Const::MKaon, 321},
				   {Const::MElectron, 11}};
			break;
		case _KaonM:
			_masspdg = {{Const::MKaon, 321},
				   {Const::MMuon, 13}};
			break;
		case _CharmE:
			_masspdg = {{Const::MDs, 431},
				   {Const::MElectron, 11}};
			break;
		case _CharmM:
			_masspdg = {{Const::MDs, 431},
				   {Const::MMuon, 13}};
			break;
		case _CharmT:
			_masspdg = {{Const::MDs, 431},
				   {Const::MTau, 15}};
			break;
		case _Kaon0E:
			_masspdg = {{Const::MKaon0, 130},
				   {Const::MPion, 211},
				   {Const::MElectron, 11}};
			break;
		case _Kaon0M:
			_masspdg = {{Const::MKaon0, 130},
				   {Const::MPion, 211},
				   {Const::MMuon, 13}};
			break;
		case _KaonCE:
			_masspdg = {{Const::MKaon, 321},
				   {Const::MPion0, 111},
				   {Const::MElectron, 11}};
			break;
		case _KaonCM:
			_masspdg = {{Const::MKaon, 321},
				   {Const::MPion0, 111},
				   {Const::MMuon, 13}};
			break;
		default:
			throw std::invalid_argument("Channel " + toString(chan) + " unknown");
	}
}

//Kinematic function
double Amplitude::Kallen(double x, double y, double z)
{
	return x*x + y*y + z*z - 2*(x*y + x*z + y*z);
}

double Amplitude::SqrtKallen(double x, double y, double z)
{
	double kal = Kallen(x, y, z);
	if (kal < 0)
		return std::sqrt(-kal);
	return std::sqrt(kal);
}

/////////////////////
//Diff decay widths//
////////////////////

//	to be integrated over ds/m0², dt/m0², da, dcosb, dc
//	ds = -2M dE
double Amplitude::dGammad5_3B(double M2)
{
	return M2 * _m_parent / (2048. * Const::pi5);	//integration over dsigma dtau
}

//	after integrating angular depend., to be integrated over ds, dt
double Amplitude::dGammad2_3B(double M2)
{
	return 8. * Const::pi2 * dGammad5_3B(M2);
}

//	to be integrated over da, dcosb		(m1/m0)²  (m2/m0)²
double Amplitude::dGammad2_2B(double M2, double x, double y)
{
	return M2 * std::sqrt(Kallen(1, x, y)) / (64. * Const::pi2 * _m_parent);
}

//	constant, after integrating angular depend.
double Amplitude::dGammad0_2B(double M2, double x, double y)
{
	return 4. * Const::pi * dGammad2_2B(M2, x, y);
}

//integration limits for three body when is constant wrt to t
double Amplitude::Limit(double &s, double x, double y, double z)
{
	double t = 0.;
	return Limit(s, t, x, y, z);
}

//integration limits for three body, where z refers to s, y referst to t, and x refers to u
double Amplitude::Limit(double &s, double &t, double x, double y, double z)
{
	double s0 = s;
	double slow = x + y + 2 * std::sqrt(x * y);
	double ssup = 1 + z - 2 * std::sqrt(z);
	s = slow + (ssup - slow) * s;	//this is s

	double kal = std::sqrt(Kallen(s, y, x) * Kallen(1, s, z));
	double tlow = z + x + ((1 - s - z) * (s - y + x) - kal) / (2 * s);
	double tsup = z + x + ((1 - s - z) * (s - y + x) + kal) / (2 * s);
	if (s == 0)
		tlow = tsup = z + x;
	t = tlow + (tsup - tlow) * t;

	return (ssup - slow) * (tsup - tlow);
}

//////////////////////
//Generic amplitudes//
//////////////////////

//					      angle lep    lepton    meson	
double Amplitude::M2_LeptonPseudoMeson(double cos0, double x, double y)
{
		// added only if Majorana
	return (IsMajorana() ? 0. : M2_LeptonPseudoMeson(-Neutrino::Helicity(_opts), cos0, x, y))
				  + M2_LeptonPseudoMeson( Neutrino::Helicity(_opts), cos0, x, y);
}

double Amplitude::M2_LeptonPseudoMeson(int hel, double cos0, double x, double y)
{
	return Const::GF2 * std::pow(_m_parent, 4) * 
		(pow(1 - x, 2) - y * (1 + x) - (1 - x) * hel * std::sqrt(Kallen(1, x, y)) * cos0);
}

//						angle lep    lepton    meson	
double Amplitude::M2_NeutrinoPseudoMeson(double cos0, double x, double y)
{
		// added only if Majorana
	return (IsMajorana() ? 0. : M2_NeutrinoPseudoMeson(-Neutrino::Helicity(_opts),
							    cos0, x, y))
				  + M2_NeutrinoPseudoMeson( Neutrino::Helicity(_opts),
						   	    cos0, x, y);
}

double Amplitude::M2_NeutrinoPseudoMeson(int hel, double cos0, double x, double y)
{
	return M2_LeptonPseudoMeson(hel, x, y, cos0) / 2.0;
}

double Amplitude::M2_LeptonVectorMeson(double cos0, double x, double y)
{
		// added only if Majorana
	return (IsMajorana() ? 0. : M2_LeptonVectorMeson(-Neutrino::Helicity(_opts), cos0, x, y))
				  + M2_LeptonVectorMeson( Neutrino::Helicity(_opts), cos0, x, y);
}

//					lepton		meson	
double Amplitude::M2_LeptonVectorMeson(int hel, double cos0, double x, double y)	//must be divided by vector meson mass
{
	return Const::GF2 * std::pow(_m_parent, 4)
	     * (pow(1 - x, 2) + y * (1 + x) - 2 * y*y
		- (1 - x - 2*y) *  hel * std::sqrt(Kallen(1, x, y)) * cos0);
}

double Amplitude::M2_NeutrinoVectorMeson(double cos0, double x, double y)
{
		// added only if Majorana
	return (IsMajorana() ? 0. : M2_NeutrinoVectorMeson(-Neutrino::Helicity(_opts),
							   cos0, x, y))
				  + M2_NeutrinoVectorMeson( Neutrino::Helicity(_opts),
						   	   cos0, x, y);
}

//					lepton		meson	angle
double Amplitude::M2_NeutrinoVectorMeson(int hel, double cos0, double x, double y)
{
	return M2_LeptonVectorMeson(hel, cos0, x, y) / 2.0;
}

double Amplitude::M2_WW(double s, double cos0, double x, double y, double z)	//gL^2 + gR^2
{
		// added only if Majorana
	return (IsMajorana() ? 0. : M2_WW(-Neutrino::Helicity(_opts), s, cos0, x, y, z))
				  + M2_WW( Neutrino::Helicity(_opts), s, cos0, x, y, z);
}
//			lepton energy is s = (p0-p2)² and cos0s the angle wrt z-axis
//					       neutrino, letpon,   lepton
double Amplitude::M2_WW(int hel, double s, double cos0s, double x, double y, double z)	//gL^2 + gR^2
{
	return 16 * Const::GF2 * std::pow(_m_parent, 4) *
	       (s - x - y) * (1 + z - s - hel * std::sqrt(Kallen(1, z, s)) * cos0s);
}

double Amplitude::M2_WZ(double u, double cos0u, double x, double y, double z)	//gL^2 + gR^2
{
		// added only if Majorana
	return (IsMajorana() ? 0. : M2_WZ(-Neutrino::Helicity(_opts), u, cos0u, x, y, z))
				  + M2_WZ( Neutrino::Helicity(_opts), u, cos0u, x, y, z);
}
//			lepton energiy is u = (p0-p1)² and cos0u the angle wrt z-axis
//			       neutrino, letpon,   lepton
double Amplitude::M2_WZ(int hel, double u, double cos0u, double x, double y, double z)	//2gL*gR
{
	return 16 * Const::GF2 * std::pow(_m_parent, 4) *
		sqrt(y * z) * (1 + x - u - hel * std::sqrt(Kallen(1, x, u)) * cos0u);
}

double Amplitude::M2_WZ(double s, double t, double cos0s, double cos0t, double x, double y, double z)	//2gL*gR
{
		// added only if Majorana
	return (IsMajorana() ? 0. : M2_WZ(-Neutrino::Helicity(_opts), s, t,
					  cos0s, cos0t, x, y, z))
				  + M2_WZ( Neutrino::Helicity(_opts), s, t,
					  cos0s, cos0t, x, y, z);
}
//			lepton energies are s = (p0-p2)², t = (p0-p3)² and cos0s,t the angles wrt z-axis
//			       neutrino, letpon,   lepton
double Amplitude::M2_WZ(int hel, double s, double t, double cos0s, double cos0t, double x, double y, double z)	//2gL*gR
{
	double u = 1 + x + y + z - s - t;
	double cos0u = (std::sqrt(Kallen(1, y, s)) * cos0s + std::sqrt(Kallen(1, z, t)) * cos0t)
		     / std::sqrt(Kallen(1, x, u));

	return M2_WZ(hel, u, cos0u, x, y, z);
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
//{
//	return (GetFermion() ? 0.0 : M2_LeptonNeutrino(-Helicity(), u, x, y, z)) +	//added only if majorana
//				     M2_LeptonNeutrino( Helicity(), u, x, y, z);
//}
//
//double Amplitude::M2_LeptonNeutrino(int hel, double u, double x, double y, double z)
{
	return 16 * Const::GF2 * std::pow(_m_parent, 4) *
		(1 + x - u) * (u - y - z - Neutrino::Helicity(_opts) * std::sqrt(Kallen(u, y, z));
}

//	This amplitude is to be used if the mixing comes from the flavour in final state
//				     neutrino  lepton    neutrino  neutrino  angle betw. lepton and neutr
//	production is from antilepeton
double Amplitude::M2_AntiLeptonNeutrino(double s, double x, double y, double z)
//{
//	return (GetFermion() ? 0.0 : M2_AntiLeptonNeutrino(-Helicity(_opts), s, x, y, z)) +	//added only if majorana
//				     M2_AntiLeptonNeutrino( Helicity(_opts), s, x, y, z);
//}
//
//double Amplitude::M2_AntiLeptonNeutrino(int hel, double s, double x, double y, double z)
{
	return 16 * Const::GF2 * std::pow(_m_parent, 4) * 
		(s - x - y) * (1 + z - s - Neutrino::Helicity(_opts) * std::sqrt(Kallen(1, s, z));
}

double Amplitude::M2_LeptonTwo(double x, double y)
//{
//	return (GetFermion() ? 0.0 : M2_LeptonTwo(-Helicity(_opts), x, y)) +	//added only if majorana
//				     M2_LeptonTwo( Helicity(_opts), x, y);
//}
////					      neutrino	meson
//double Amplitude::M2_LeptonTwo(int hel, double x, double y)	//y is the meson
{
	return Const::GF2 * std::pow(_m_parent, 4) * 
		(std::pow(1 - x, 2) - y * (1 + x)
		 - (1 - x) * Neutrino::Helicity(_opts) * std::sqrt(Kallen(1, x, y)));
}

double Amplitude::M2_LeptonThree(double x, double y, double z)
//{
//	return (GetFermion() ? 0.0 : M2_LeptonThree(-Helicity(_opts), x, y, z)) +	//added only if majorana
//				     M2_LeptonThree( Helicity(_opts), x, y, z);
//}
////	not implemented
//double Amplitude::M2_LeptonThree(int hel, double x, double y, double z)	
{
	return Const::GF2 * std::pow(_m_parent, 4);
}

double Amplitude::M2_MesonTwo(double x, double y)
//{
//	return (GetFermion() ? 0.0 : M2_MesonTwo(-Helicity(_opts), x, y)) +	//added only if majorana
//				     M2_MesonTwo( Helicity(_opts), x, y);
//}
////					   neutrino  lepton
//double Amplitude::M2_MesonTwo(int hel, double x, double y)
{
	return Const::GF2 * std::pow(_m_parent, 4) * 
		(x + y - std::pow(x - y, 2)
		 - Neutrino::Helicity(_opts) * (y - x) * std::sqrt(Kallen(1, x, y)));
}

double Amplitude::M2_MesonThree(double s, double t, double x, double y, double z, double L_, double L0)
//{
//	return (GetFermion() ? 0.0 : M2_MesonThree(-Helicity(_opts), s, t, x, y, z, L_, L0)) +	//added only if majorana
//				     M2_MesonThree( Helicity(_opts), s, t, x, y, z, L_, L0);
//}
////Jackson frame??
//double Amplitude::M2_MesonThree(int hel, double s, double t, double x, double y, double z, double L_, double L0)
{
	double u = 1 + x + y + z - s - t;

	double F = 2 * (1 + L_ * u / x);
	double G = (1 + L_ * u / x) - (L_ - L0) * (1 + 1 / x);

	double A = (1 + y - t)
		 * (1 + z - s - Neutrino::Helicity(_opts) * std::sqrt(Kallen(1, z, s)))
		 - (u - y - z - Neutrino::Helicity(_opts) * std::sqrt(Kallen(u, y, z)));
	double B = (y + z) * (u - y - z) + 4 * y * z -
	           (y - z) * Neutrino::Helicity(_opts) * std::sqrt(Kallen(u, y, z));
	double C = (1 + y - t) * 2*z + (1 + z - s)
		 * (2*y + Neutrino::Helicity(_opts) * std::sqrt(Kallen(u, y, z)))
		 - Neutrino::Helicity(_opts) * (u - z + y) * std::sqrt(Kallen(1, z, s));

	return Const::GF2 * std::pow(_m_parent, 4) * ( (F*F) * A + (G*G) * B - (F*G) * C );
}

void Amplitude::SetNeutrino(const Neutrino &N)
{
	if (N != _N)
		Reset();

	_N = N;
}
