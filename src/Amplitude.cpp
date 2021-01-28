/* 
 * unpolarised amplitudes for processes involving N
 * also other util procedures for derived classes
 * and definition of decay width from unpolarised amplitude
 */

#include "physics/Amplitude.h"

Amplitude::Amplitude(Neutrino N) :
	_m_parent(0),
	_channel(Channel::undefined),
	_N(std::move(N))
{
}

void Amplitude::LoadMass(Channel::Name chan)
{
	switch(chan)
	{
		//masses in channelname order
		case Channel::ALL:
		case Channel::nnn:
			_masspdg = {{Const::MNeutrino, 12},
				   {Const::MNeutrino, -12},
				   {Const::MNeutrino, 12}};
			break;
		case Channel::nGAMMA:
			_masspdg = {{Const::MNeutrino, 12},
				   {Const::MPhoton, 22}};
			break;
		case Channel::ExpALL:
			//neutrino lepton lepton AA
		case Channel::nEE:
			_masspdg = {{Const::MNeutrino, 12},
				   {Const::MElectron, -11},
				   {Const::MElectron, 11}};
			break;
		case Channel::nMM:
			_masspdg = {{Const::MNeutrino, 12},
				   {Const::MMuon, -13},
				   {Const::MMuon, 13}};
			break;
			//neutrino lepton lepton AB
		case Channel::nEM:
			_masspdg = {{Const::MNeutrino, 12},
				   {Const::MElectron, -11},
				   {Const::MMuon, 13}};
			break;
		case Channel::nET:
			_masspdg = {{Const::MNeutrino, 12},
				   {Const::MElectron, -11},
				   {Const::MTau, 15}};
			break;
		case Channel::nMT:
			_masspdg = {{Const::MNeutrino, 12},
				   {Const::MMuon, -13},
				   {Const::MTau, 15}};
			break;
			//neutrino psuedomeson
		case Channel::nPI0:
			_masspdg = {{Const::MNeutrino, 12},
				   {Const::MPion0, 111}};
			break;
		case Channel::nETA:
			_masspdg = {{Const::MNeutrino, 12},
				   {Const::MEta, 221}};
			break;
		case Channel::nETAi:
			_masspdg = {{Const::MNeutrino, 12},
				   {Const::MEtai, 331}};
			break;
			//lepton psuedomeson
		case Channel::EPI:
			_masspdg = {{Const::MElectron, 11},
				   {Const::MPion, 211}};
			break;
		case Channel::MPI:
			_masspdg = {{Const::MMuon, 13},
				   {Const::MPion, 211}};
			break;
		case Channel::TPI:
			_masspdg = {{Const::MTau, 15},
				   {Const::MPion, 211}};
			break;
		case Channel::EKA:
			_masspdg = {{Const::MElectron, 11},
				   {Const::MKaon, 321}};
			break;
		case Channel::MKA:
			_masspdg = {{Const::MMuon, 13},
				   {Const::MKaon, 321}};
			break;
		case Channel::ECHARM:
			_masspdg = {{Const::MElectron, 12},
				   {Const::MD, 411}};
			break;
			//neutrino vectormeson
		case Channel::nRHO0:
			_masspdg = {{Const::MNeutrino, 12},
				   {Const::MRho0, 113}};
			break;
		case Channel::nOMEGA:
			_masspdg = {{Const::MNeutrino, 12},
				   {Const::MOmega, 223}};
			break;
		case Channel::nPHI:
			_masspdg = {{Const::MNeutrino, 12},
				   {Const::MPhi, 333}};
			break;
			//lepton vectormeson
		case Channel::ERHO:
			_masspdg = {{Const::MElectron, 11},
				   {Const::MRho, 213}};
			break;
		case Channel::MRHO:
			_masspdg = {{Const::MMuon, 13},
				   {Const::MRho0, 213}};
			break;
		case Channel::EKAx:
			_masspdg = {{Const::MElectron, 11},
				   {Const::MKaonx, 9000321}};
			break;
		case Channel::MKAx:
			_masspdg = {{Const::MMuon, 13},
				   {Const::MKaonx, 9000321}};
			break;
		//PRODUCTION
		//Parent first
		case Channel::MuonE:
			_masspdg = {{Const::MMuon, 13},
				   {Const::MElectron, 11},
				   {Const::MNeutrino, 14}};
		case Channel::MuonM:
			_masspdg = {{Const::MMuon, 13},
				   {Const::MElectron, 11},
				   {Const::MNeutrino, -12}};
			break;
		case Channel::TauEE:
			_masspdg = {{Const::MTau, 15},
				   {Const::MElectron, 11},
				   {Const::MNeutrino, 16}};
		case Channel::TauET:
			_masspdg = {{Const::MTau, 15},
				   {Const::MElectron, 11},
				   {Const::MNeutrino, -12}};
			break;
		case Channel::TauMM:
			_masspdg = {{Const::MTau, 15},
				   {Const::MMuon, 13},
				   {Const::MNeutrino, 16}};
		case Channel::TauMT:
			_masspdg = {{Const::MTau, 15},
				   {Const::MMuon, 13},
				   {Const::MNeutrino, -14}};
			break;
		case Channel::TauPI:
			_masspdg = {{Const::MTau, 15},
				   {Const::MPion, 211}};
			break;
		case Channel::Tau2PI:
			_masspdg = {{Const::MTau, 15},
				   {Const::MPion, -211},
				   {Const::MPion0, 111}};
			break;
		case Channel::PionE:
			_masspdg = {{Const::MPion, 211},
				   {Const::MElectron, 11}};
			break;
		case Channel::PionM:
			_masspdg = {{Const::MPion, 211},
				   {Const::MMuon, 13}};
			break;
		case Channel::KaonE:
			_masspdg = {{Const::MKaon, 321},
				   {Const::MElectron, 11}};
			break;
		case Channel::KaonM:
			_masspdg = {{Const::MKaon, 321},
				   {Const::MMuon, 13}};
			break;
		case Channel::CharmE:
			_masspdg = {{Const::MDs, 431},
				   {Const::MElectron, 11}};
			break;
		case Channel::CharmM:
			_masspdg = {{Const::MDs, 431},
				   {Const::MMuon, 13}};
			break;
		case Channel::CharmT:
			_masspdg = {{Const::MDs, 431},
				   {Const::MTau, 15}};
			break;
		case Channel::Kaon0E:
			_masspdg = {{Const::MKaon0, 130},
				   {Const::MPion, 211},
				   {Const::MElectron, -11}};
			break;
		case Channel::Kaon0M:
			_masspdg = {{Const::MKaon0, 130},
				   {Const::MPion, 211},
				   {Const::MMuon, -13}};
			break;
		case Channel::KaonCE:
			_masspdg = {{Const::MKaon, 321},
				   {Const::MPion0, 111},
				   {Const::MElectron, -11}};
			break;
		case Channel::KaonCM:
			_masspdg = {{Const::MKaon, 321},
				   {Const::MPion0, 111},
				   {Const::MMuon, -13}};
			break;
		default:
			throw std::invalid_argument("Channel " + Channel::toString(chan) + " unknown");
	}
}

double Amplitude::MassThreshold(Channel::Name chan)
{
	if (_channel != chan) {
		LoadMass(chan);
		_channel = chan;
	}

	switch (Channel::whichType(chan)) {
		case Channel::Type::decayrates:
			return std::accumulate(_masspdg.begin(), _masspdg.end(), 0.,
					[](double sum, const std::pair<double, int> mp) {
						return sum + mp.first; });
		case Channel::Type::production:
			return std::accumulate(_masspdg.begin()+1, _masspdg.end(),
								_masspdg.front().first,
					[](double sum, const std::pair<double, int> mp) {
						return sum - mp.first; });
		default:
			return 0.;
	}
}

//Kinematic function
double Amplitude::Kallen(double x, double y, double z)
{
	return std::max(0., x*x + y*y + z*z - 2*(x*y + x*z + y*z));
}

/////////////////////
//Diff decay widths//
////////////////////

//	to be integrated over ds/m0², dt/m0², da, dcosb, dc
//	ds = -2M dE
double Amplitude::dGammad5_3B()
{
	return _m_parent / (2048. * Const::pi5);	//integration over dsigma dtau
}

//	after integrating angular depend., to be integrated over ds, dt
double Amplitude::dGammad2_3B()
{
	return 8. * Const::pi2 * dGammad5_3B();
}

//	to be integrated over da, dcosb		(m1/m0)²  (m2/m0)²
double Amplitude::dGammad2_2B(double x, double y)
{
	return std::sqrt(Kallen(1, x, y)) / (64. * Const::pi2 * _m_parent);
}

//	constant, after integrating angular depend.
double Amplitude::dGammad0_2B(double x, double y)
{
	return 4. * Const::pi * dGammad2_2B(x, y);
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
	double slow = x  + y + 2. * std::sqrt(x * y);
	double ssup = 1. + z - 2. * std::sqrt(z);
	s = slow + (ssup - slow) * s;	//this is s

	if (s == 0.) {
		t = z + x;
		return 0.;
	}

	double kal = std::sqrt(Kallen(s, y, x) * Kallen(1., s, z));
	double tlow = z + x + ((1. - s - z) * (s - y + x) - kal) / (2. * s);
	//double tsup = z + x + ((1. - s - z) * (s - y + x) + kal) / (2. * s);
	t = tlow + kal * t / s;	// this is t

	return (ssup - slow) * kal / s;
}

//////////////////////
//Generic amplitudes//
//////////////////////

//					      angle lep    lepton    meson	
double Amplitude::M2_LeptonPseudoMeson(double cos0, double x, double y)
{
		// added only if Majorana
	int hel = _N.IsParticle() ? _N.Helicity() : -_N.Helicity();
	return (!_N.IsMajorana() ? 0. : M2_LeptonPseudoMeson(-hel, cos0, x, y))
				      + M2_LeptonPseudoMeson( hel, cos0, x, y);
}

double Amplitude::M2_LeptonPseudoMeson(int hel, double cos0, double x, double y)
{
	return Const::GF2 * std::pow(_m_parent, 4) * 
		((1 - x) * (1 - x) - y * (1 + x)
		 - (1 - x) * hel * std::sqrt(Kallen(1, x, y)) * cos0);
}

//						angle lep    lepton    meson	
double Amplitude::M2_NeutrinoPseudoMeson(double cos0, double x, double y)
{
		// added only if Majorana
	int hel = _N.IsParticle() ? _N.Helicity() : -_N.Helicity();
	return (!_N.IsMajorana() ? 0. : M2_NeutrinoPseudoMeson(-hel, cos0, x, y))
				      + M2_NeutrinoPseudoMeson( hel, cos0, x, y);
}

double Amplitude::M2_NeutrinoPseudoMeson(int hel, double cos0, double x, double y)
{
	return M2_LeptonPseudoMeson(hel, x, y, cos0) / 2.0;
}

double Amplitude::M2_LeptonVectorMeson(double cos0, double x, double y)
{
		// added only if Majorana
	int hel = _N.IsParticle() ? _N.Helicity() : -_N.Helicity();
	return (!_N.IsMajorana() ? 0. : M2_LeptonVectorMeson(-hel, cos0, x, y))
				      + M2_LeptonVectorMeson( hel, cos0, x, y);
}

//					lepton		meson	
double Amplitude::M2_LeptonVectorMeson(int hel, double cos0, double x, double y)	//must be divided by vector meson mass
{
	return Const::GF2 * std::pow(_m_parent, 4)
	     * ((1 - x) * (1 - x) + y * (1 + x) - 2 * y * y
		- (1 - x - 2 * y) *  hel * std::sqrt(Kallen(1, x, y)) * cos0);
}

double Amplitude::M2_NeutrinoVectorMeson(double cos0, double x, double y)
{
	int hel = _N.IsParticle() ? _N.Helicity() : -_N.Helicity();
	return (!_N.IsMajorana() ? 0. : M2_NeutrinoVectorMeson(-hel, cos0, x, y))
				      + M2_NeutrinoVectorMeson( hel, cos0, x, y);
}

//					lepton		meson	angle
double Amplitude::M2_NeutrinoVectorMeson(int hel, double cos0, double x, double y)
{
	return M2_LeptonVectorMeson(hel, cos0, x, y) / 2.0;
}

double Amplitude::M2_WW(double s, double cos0, double x, double y, double z)	//gL^2 + gR^2
{
	int hel = _N.IsParticle() ? _N.Helicity() : -_N.Helicity();
	return (!_N.IsMajorana() ? 0. : M2_WW(-hel, s, cos0, x, y, z))
				      + M2_WW( hel, s, cos0, x, y, z);
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
	int hel = _N.IsParticle() ? _N.Helicity() : -_N.Helicity();
	return (!_N.IsMajorana() ? 0. : M2_WZ(-hel, u, cos0u, x, y, z))
				      + M2_WZ( hel, u, cos0u, x, y, z);
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
	int hel = _N.IsParticle() ? _N.Helicity() : -_N.Helicity();
	return (!_N.IsMajorana() ? 0. : M2_WZ(-hel, s, t, cos0s, cos0t, x, y, z))
				      + M2_WZ( hel, s, t, cos0s, cos0t, x, y, z);
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
{
	int hel = _N.IsParticle() ? _N.Helicity() : -_N.Helicity();
	return (hel ? 1. : 2.) * M2_LeptonNeutrino(hel, u, x, y, z);
}

double Amplitude::M2_LeptonNeutrino(int hel, double u, double x, double y, double z)
{
	return 16. * Const::GF2 * std::pow(_m_parent, 4) *
		(1. + x - u) * (u - y - z - hel * std::sqrt(Kallen(u, y, z)));
}

//	This amplitude is to be used if the mixing comes from the flavour in final state
//	production is from antilepeton
//			     neutrino  lepton    neutrino  neutrino  angle betw. lepton and neutr
double Amplitude::M2_AntileptonNeutrino(double s, double x, double y, double z)
{
	int hel = _N.IsParticle() ? _N.Helicity() : -_N.Helicity();
	return (hel ? 1. : 2.) * M2_AntileptonNeutrino(hel, s, x, y, z);
}

double Amplitude::M2_AntileptonNeutrino(int hel, double s, double x, double y, double z)
{
	return 16. * Const::GF2 * std::pow(_m_parent, 4.) * 
		(s - x - y) * (1. + z - s - hel * std::sqrt(Kallen(1., s, z)));
}

double Amplitude::M2_LeptonTwo(double x, double y)
{
	int hel = _N.IsParticle() ? _N.Helicity() : -_N.Helicity();
	return (hel ? 1. : 2.) * M2_LeptonTwo(hel, x, y);
}

double Amplitude::M2_LeptonTwo(int hel, double x, double y)
{
	return Const::GF2 * std::pow(_m_parent, 4) * 
		((1 - x) * (1 - x) - y * (1 + x) - (1 - x) * hel * std::sqrt(Kallen(1, x, y)));
}

double Amplitude::M2_LeptonThree(double x, double y, double z)
{
	int hel = _N.IsParticle() ? _N.Helicity() : -_N.Helicity();
	return (hel ? 1. : 2.) * M2_LeptonThree(hel, x, y);
}

double Amplitude::M2_LeptonThree(int hel, double x, double y, double z)
{
	return Const::GF2 * std::pow(_m_parent, 4);
}

double Amplitude::M2_MesonTwo(double x, double y)
{
	int hel = _N.IsParticle() ? _N.Helicity() : -_N.Helicity();
	return (hel ? 1. : 2.) * M2_MesonTwo(hel, x, y);
}

double Amplitude::M2_MesonTwo(int hel, double x, double y)
{
	return Const::GF2 * std::pow(_m_parent, 4) * 
		(x + y - std::pow(x - y, 2) + hel * (x - y) * std::sqrt(Kallen(1, x, y)));
}

double Amplitude::M2_MesonThree(double s, double t, double x, double y, double z, double L_, double L0)
{
	int hel = _N.IsParticle() ? _N.Helicity() : -_N.Helicity();
	return (hel ? 1. : 2.) * M2_MesonThree(hel, s, t, x, y, z, L_, L0);
}

double Amplitude::M2_MesonThree(int hel, double s, double t, double x, double y, double z, double L_, double L0)
{
	double u = 1. + x + y + z - s - t;

	double F = 2. * (1 + L_ * u / x);
	double G = (1. + L_ * u / x) - (L_ - L0) * (1. + 1. / x);

	double A = (1. + y - t) * (1. + z - s - hel * std::sqrt(Kallen(1., z, s)))
		 	         - (u - y - z - hel * std::sqrt(Kallen(u, y, z)));
	double B = (y + z) * (u - y - z) + 4. * y * z
		 - hel * (y - z) * std::sqrt(Kallen(u, y, z));
	double C = 2 * z * (1. + y - t)
		 + (1. + z - s) * (2. * y + hel * std::sqrt(Kallen(u, y, z)))
		 - hel * (u - z + y) * std::sqrt(Kallen(1., z, s));

	return Const::GF2 * std::pow(_m_parent, 4) * ( (F*F) * A + (G*G) * B - (F*G) * C );
}

void Amplitude::SetNeutrino(const Neutrino &N)
{
	if (N != _N)
		Reset();

	_N = N;
}
