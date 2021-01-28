#include "physics/PhaseSpace.h"

// using external MT generator
PhaseSpace::PhaseSpace(Neutrino N) : Amplitude(std::move(N)),
	_genps(new TGenPhaseSpace())
{
}

PhaseSpace::~PhaseSpace()
{
	delete _genps;
}

std::pair<std::vector<Particle>, double> PhaseSpace::Generate(Channel::Name chan,
					const Particle &part, const Mixing &mix) {
	TLorentzVector frame = static_cast<TLorentzVector>(part);
	return Generate(chan, frame, mix);

}

std::pair<std::vector<Particle>, double> PhaseSpace::Generate(Channel::Name chan,
				const TLorentzVector &frame, const Mixing &mix) {
	std::vector<Particle> parts;
	if (!SetDecay(chan)) // decay not allowed
		return std::make_pair(parts, 0.); 	// return empty vector

	size_t cc = 0;
	std::uniform_real_distribution<> rndm;
	double weight = _genps->Generate();	// generate phase space
	double val = Gamma(chan, mix); 
	if (val <= 1e-12 || val - 1e-12 > 1.)	// phase space impossible
		return std::make_pair(parts, 0.);

	while (rndm(RNG::_mt) > val && cc++ < 1e4) {
		weight = _genps->Generate();	// generate phase space
		val = Gamma(chan, mix);
	}

	if (cc > 1e4)
		return std::make_pair(parts, 0.);

	// create parent frame
	TLorentzVector vec = *(_genps->GetDecay(0));
	vec.Boost(frame.BoostVector());

	// allocate space for particles
	auto psgs = Process::Pdg(chan);
	parts.reserve(pdgs.size());

	// reverse pdg sign if particle is "anti"
	int sign = 1;
	switch (Channel::whichType(chan))
	{
		case Channel::Type::decayrates:
			// assign sign during decay
			if (_N.IsMajorana()) {
				std::bernoulli_distribution b(0.5);
				sign = b(RNG::_mt) ? 1 : -1;
			}
			else // is dirac
				sign = _N.IsParticle() ? 1 : -1;
			parts.push_back(Particle(sign * pdgs.front(), vec));
			break;
		// don't care about sign in production
		case Channel::Type::production:
			parts.push_back(Particle(sign * _N.Pdg(), vec));
			break;
		default:
			throw std::invalid_argument("Channel " + Channel::toString(chan) + " unknown");
	}

	for (size_t i = 1; i < pdgs.size(); ++i) {
		TLorentzVector vec = *(_genps->GetDecay(i));
		vec.Boost(frame.BoostVector());
		parts.push_back(Particle(sign * pdgs[i], vec));
	}

	return std::make_pair(parts, weight);
}

bool PhaseSpace::SetDecay(Channel::Name chan)
{
	auto mass = Process::Mass(chan);
	double ps_masses[10];
	switch (Channel::whichType(chan))
	{
		case Channel::Type::decayrates:
			ps_masses[0] = mass.front();
			_m_parent = _N.M();
			break;
		case Channel::Type::production:
			ps_masses[0] = _N.M();
			_m_parent = mass.front();
			break;
		default:
			throw std::invalid_argument("Channel " + Channel::toString(chan) + " unknown");
	}

	for (size_t i = 1; i < _masspdg.size(); ++i)
		ps_masses[i] = mass[i];

	// compute decay in rest frame, boost later
	TLorentzVector rest(0, 0, 0, _m_parent);

	return _genps->SetDecay(rest, mass.size(), ps_masses);
}

double PhaseSpace::Gamma(Channel::Name chan, const Mixing &mix)
{
	using GammaR = double (PhaseSpace::*)(const Mixing &mix);
	GammaR gr;

	switch (chan)
	{
		case Channel::ALL:
		case Channel::nnn:
		case Channel::nGAMMA:
			return 1.;
		case Channel::nEE:
			gr = &PhaseSpace::nEE;
			break;
		case Channel::nEM:
			gr = &PhaseSpace::nEM;
			break;
		case Channel::nMM:
			gr = &PhaseSpace::nMM;
			break;
		case Channel::nET:
			gr = &PhaseSpace::nET;
			break;
		case Channel::nMT:
			gr = &PhaseSpace::nMT;
			break;
		case Channel::nPI0:
			gr = &PhaseSpace::nPI0;
			break;
		case Channel::EPI:
			gr = &PhaseSpace::EPI;
			break;
		case Channel::MPI:
			gr = &PhaseSpace::MPI;
			break;
		case Channel::TPI:
			gr = &PhaseSpace::TPI;
			break;
		case Channel::EKA:
			gr = &PhaseSpace::EKA;
			break;
		case Channel::MKA:
			gr = &PhaseSpace::MKA;
			break;
		case Channel::nRHO0:
			gr = &PhaseSpace::nRHO0;
			break;
		case Channel::ERHO:
			gr = &PhaseSpace::ERHO;
			break;
		case Channel::MRHO:
			gr = &PhaseSpace::MRHO;
			break;
		case Channel::EKAx:
			gr = &PhaseSpace::EKAx;
			break;
		case Channel::MKAx:
			gr = &PhaseSpace::MKAx;
			break;
		case Channel::nOMEGA:
			gr = &PhaseSpace::nOMEGA;
			break;
		case Channel::nETA:
			gr = &PhaseSpace::nETA;
			break;
		case Channel::nETAi:
			gr = &PhaseSpace::nETAi;
			break;
		case Channel::nPHI:
			gr = &PhaseSpace::nPHI;
			break;
		case Channel::ECHARM:
			gr = &PhaseSpace::ECHARM;
			break;
		case Channel::MuonE:	//these are needed as well
			gr = &PhaseSpace::MuonE;
			break;
			return 1.;
		case Channel::MuonM:
			gr = &PhaseSpace::MuonM;
			break;
		case Channel::TauEE:
			gr = &PhaseSpace::TauEE;
			break;
		case Channel::TauET:
			gr = &PhaseSpace::TauET;
			break;
		case Channel::TauMM:
			gr = &PhaseSpace::TauMM;
			break;
		case Channel::TauMT:
			gr = &PhaseSpace::TauMT;
			break;
		case Channel::TauPI:				//1
			//gr = &PhaseSpace::TauPI;
			//break;
			return 1.;
		case Channel::Tau2PI:				//1
			//gr = &PhaseSpace::Tau2PI;
			//break;
			return 1.;
		case Channel::PionE:				//1
			//gr = &PhaseSpace::PionE;
			//break;
			return 1.;
		case Channel::PionM:				//1
			//gr = &PhaseSpace::PionM;
			//break;
			return 1.;
		case Channel::KaonE:				//1
			//gr = &PhaseSpace::KaonE;
			//break;
			return 1.;
		case Channel::KaonM:				//1
			//gr = &PhaseSpace::KaonM;
			//break;
			return 1.;
		case Channel::CharmE:				//1
			//gr = &PhaseSpace::CharmE;
			//break;
			return 1.;
		case Channel::CharmM:				//1
			//gr = &PhaseSpace::CharmM;
			//break;
			return 1.;
		case Channel::CharmT:				//1
			//gr = &PhaseSpace::CharmT;
			//break;
			return 1.;
		case Channel::Kaon0E:
			gr = &PhaseSpace::Kaon0E;
			break;
		case Channel::Kaon0M:
			gr = &PhaseSpace::Kaon0M;
			break;
		case Channel::KaonCE:
			gr = &PhaseSpace::KaonCE;
			break;
		case Channel::KaonCM:
			gr = &PhaseSpace::KaonCM;
			break;
		default:
			throw std::invalid_argument("Channel " + Channel::toString(chan) + " unknown");
			break;
	}
	
	return (this->*gr)(mix);
}

//Neutrino LeptonLepton AA

double PhaseSpace::nEE(const Mixing &mix)
{
	auto gamma = NeutrinoLeptonAA(Channel::nEE, Channel::nEE_o, _masspdg[0].first, _masspdg[1].first);

	if (gamma.first < 0.)
		gamma.first = 0.;
	if (gamma.second < 0.)
		gamma.second = 0.;

	//return (maxName_u > 1e-25 || maxName_o > 1e-25) ? 
	return ( gamma.first * mix.Ue(2) + gamma.second * (mix.Um(2) + mix.Ut(2)))
	     / (_table[Channel::nEE] * mix.Ue(2) + _table[Channel::nEE_o] * (mix.Um(2) + mix.Ut(2)));
}

double PhaseSpace::nMM(const Mixing &mix)
{
	auto gamma = NeutrinoLeptonAA(Channel::nMM, Channel::nMM_o, _masspdg[0].first, _masspdg[1].first);

	if (gamma.first < 0.)
		gamma.first = 0.;

	if (gamma.second < 0.)
		gamma.second = 0.;

	//return (maxName_u > 1e-25 || maxName_o > 1e-25) ? 
	return ( gamma.first * mix.Um(2) + gamma.second * (mix.Ue(2) + mix.Ut(2)))
	     / (_table[Channel::nMM] * mix.Um(2) + _table[Channel::nMM_o] * (mix.Ue(2) + mix.Ut(2)));
}


std::pair<double, double> PhaseSpace::NeutrinoLeptonAA(Channel::Name chan, Channel::Name chan_o, double m_neut, double m_lepton)
{
	_m_parent = _N.M();
	double dMn2 = std::pow(m_neut / _m_parent, 2);
	double dML2 = std::pow(m_lepton / _m_parent, 2);

	//neutrino flavour is the same as leptons -> both Z and W interaction
	//double gL_CC = -0.5 + Const::sin2W + (2*GetParticle() - 1);	//times U(lepton flavour)
	double gL_CC = -0.5 + Const::sin2W + 1;	//times U(lepton flavour)
	double gR_CC = Const::sin2W;

	//neutrino flavour is different from leptons -> Z interaction
	double gL_NC = -0.5 + Const::sin2W;	//times U(neutrino flavour)
	double gR_NC = Const::sin2W;

	// compute max first, if needed
	if (!_table.count(chan)) {
		_table[chan]   = NeutrinoLeptonLepton_max(dMn2, dML2, dML2, gL_CC, gR_CC);
		_table[chan_o] = NeutrinoLeptonLepton_max(dMn2, dML2, dML2, gL_NC, gR_NC);
	}

	std::array<double, 6> kine = Kinematic_3B();
						   // s      u       cos0s     cos0u
	return std::make_pair(NeutrinoLeptonLepton(kine[0], kine[2], kine[3], kine[5], dMn2, dML2, dML2, gL_CC, gR_CC),
			      NeutrinoLeptonLepton(kine[0], kine[2], kine[3], kine[5], dMn2, dML2, dML2, gL_NC, gR_NC));
}



//Neutrino LeptonLepton AB

double PhaseSpace::nEM(const Mixing &mix)
{
	auto gamma = NeutrinoLeptonAB(Channel::nEM, Channel::nEM_o, _masspdg[0].first, _masspdg[1].first, _masspdg[2].first);

	if (gamma.first < 0.)
		gamma.first = 0.;
	if (gamma.second < 0.)
		gamma.second = 0.;

	//return (maxName_u > 1e-25 || maxName_o > 1e-25) ? 
	return ( gamma.first * mix.Ue(2) + gamma.second * mix.Um(2))
	     / (_table[Channel::nEM] * mix.Ue(2) + _table[Channel::nEM_o] * mix.Um(2));
}

double PhaseSpace::nET(const Mixing &mix)
{
	auto gamma = NeutrinoLeptonAB(Channel::nET, Channel::nET_o, _masspdg[0].first, _masspdg[1].first, _masspdg[2].first);

	if (gamma.first < 0.)
		gamma.first = 0.;
	if (gamma.second < 0.)
		gamma.second = 0.;

	//return (maxName_u > 1e-25 || maxName_o > 1e-25) ? 
	return ( gamma.first * mix.Ue(2) + gamma.second * mix.Ut(2))
	     / (_table[Channel::nET] * mix.Ue(2) + _table[Channel::nET_o] * mix.Ut(2));
}

double PhaseSpace::nMT(const Mixing &mix)
{
	// computes max too
	auto gamma = NeutrinoLeptonAB(Channel::nMT, Channel::nMT_o, _masspdg[0].first, _masspdg[1].first, _masspdg[2].first);

	if (gamma.first < 0.)
		gamma.first = 0.;
	if (gamma.second < 0.)
		gamma.second = 0.;

	//return (maxName_u > 1e-25 || maxName_o > 1e-25) ? 
	return (gamma.first * mix.Um(2) + gamma.second * mix.Ut(2))
	     / (_table[Channel::nMT] * mix.Um(2) + _table[Channel::nMT_o] * mix.Ut(2));
}


std::pair<double, double> PhaseSpace::NeutrinoLeptonAB(Channel::Name chan, Channel::Name chan_o, double m_neut, double m_leptonA, double m_leptonB)
{
	_m_parent = _N.M();
	double dMn2 = std::pow(m_neut / _m_parent, 2);
	double dMA2 = std::pow(m_leptonA / _m_parent, 2);
	double dMB2 = std::pow(m_leptonB / _m_parent, 2);

	//double gL = 1.0;
	//double gR = 0.0;

	if (!_table.count(chan)) {
		_table[chan]   = NeutrinoLeptonLepton_max(dMn2, dMA2, dMB2, 1., 0.);
		_table[chan_o] = NeutrinoLeptonLepton_max(dMn2, dMB2, dMA2, 1., 0.);
	}

	std::array<double, 6> kine = Kinematic_3B();
						   // s      u       cos0s     cos0u
	return std::make_pair(NeutrinoLeptonLepton(kine[0], kine[2], kine[3], kine[5], dMn2, dMA2, dMB2, 1., 0.),
			      NeutrinoLeptonLepton(kine[2], kine[0], kine[5], kine[3], dMn2, dMB2, dMA2, 1., 0.));
}


// common max functions

double PhaseSpace::NeutrinoLeptonLepton_max(double x, double y, double z, double gL, double gR)
{
	_vars = {x, y, z, gL, gR};

	auto func = [&](double p[]) { return F_NeutrinoLeptonLepton_max(p); };
	//std::bind(&PhaseSpace::F_NeutrinoLeptonLepton_max, this, std::placeholders::_1);
	auto res = Optimization::NelderMead<double>(func, 4); //-1 to invert function
	return -func(&res[0]);
}

// returns a reversed sign for minimization algorithm
double PhaseSpace::F_NeutrinoLeptonLepton_max(double p[])
{
	double x  = _vars[0];
	double y  = _vars[1];
	double z  = _vars[2];
	double gL = _vars[3];
	double gR = _vars[4];

	double s_   = p[0];
	double u_   = p[1];
	double cos0 = p[2];
	double cos1 = p[3];

	if (std::abs(cos0) <= 1 && std::abs(cos1) <= 1)
		return -NeutrinoLeptonLepton(s_, u_, cos0, cos1, x, y, z, gL, gR);
	return 0.;
}

double PhaseSpace::NeutrinoLeptonLepton(double s, double u, double cos0, double cos1,
					double x, double y, double z,
					double gL, double gR)
{
	return dGammad5_3B() * (gL * gL + gR * gR) * M2_WW(s, cos0, x, y, z) +
		                     (2 * gL * gR) * M2_WZ(u, cos1, x, y, z);
}




//neutrino psuedomeson

double PhaseSpace::nPI0(const Mixing &mix)
{
	return NeutrinoPseudoMeson(Channel::nPI0, _masspdg[0].first, _masspdg[1].first) / _table[Channel::nPI0];
}

double PhaseSpace::nETA(const Mixing &mix)
{
	return NeutrinoPseudoMeson(Channel::nETA, _masspdg[0].first, _masspdg[1].first) / _table[Channel::nETA];
}

double PhaseSpace::nETAi(const Mixing &mix)
{
	return NeutrinoPseudoMeson(Channel::nETAi, _masspdg[0].first, _masspdg[1].first) / _table[Channel::nETAi];
}

//
double PhaseSpace::NeutrinoPseudoMeson(Channel::Name chan, double m_neut, double m_meson)
{
	_m_parent = _N.M();
	double dMn2 = std::pow(m_neut / _m_parent, 2);
	double dMM2 = std::pow(m_meson / _m_parent, 2);

	_M2_F = &Amplitude::M2_NeutrinoPseudoMeson;

	if (!_table.count(chan)) {
		_table[chan] = ToPseudoMeson_max(dMn2, dMM2);
	}

	double cos0_ = Kinematic_2B();

	return ToPseudoMeson(cos0_, dMn2, dMM2);
}	//     2 is factor from decay constant which is sqrt(2) wrt to charged meson



//lepton psuedomeson

double PhaseSpace::EPI(const Mixing &mix)
{
	return LeptonPseudoMeson(Channel::EPI, _masspdg[0].first, _masspdg[1].first) / _table[Channel::EPI];
}

double PhaseSpace::MPI(const Mixing &mix)
{
	return LeptonPseudoMeson(Channel::MPI, _masspdg[0].first, _masspdg[1].first) / _table[Channel::MPI];
}

double PhaseSpace::TPI(const Mixing &mix)
{
	return LeptonPseudoMeson(Channel::TPI, _masspdg[0].first, _masspdg[1].first) / _table[Channel::TPI];
}

double PhaseSpace::EKA(const Mixing &mix)
{
	return LeptonPseudoMeson(Channel::EKA, _masspdg[0].first, _masspdg[1].first) / _table[Channel::EKA];
}

double PhaseSpace::MKA(const Mixing &mix)
{
	return LeptonPseudoMeson(Channel::MKA, _masspdg[0].first, _masspdg[1].first) / _table[Channel::MKA];
}

double PhaseSpace::ECHARM(const Mixing &mix)
{
	return LeptonPseudoMeson(Channel::ECHARM, _masspdg[0].first, _masspdg[1].first) / _table[Channel::ECHARM];
}

double PhaseSpace::LeptonPseudoMeson(Channel::Name chan, double m_lepton, double m_meson)
{
	_m_parent = _N.M();
	double dML2 = std::pow(m_lepton / _m_parent, 2);
	double dMM2 = std::pow(m_meson / _m_parent, 2);

	_M2_F = &Amplitude::M2_LeptonPseudoMeson;

	if (!_table.count(chan))
		_table[chan] = ToPseudoMeson_max(dML2, dMM2);

	double cos0_ = Kinematic_2B();

	return ToPseudoMeson(cos0_, dML2, dMM2);
}


// common for decays into pseudomeson + lepton

double PhaseSpace::ToPseudoMeson_max(double x, double y)
{
	_vars = {x, y};

	auto func = std::bind(&PhaseSpace::ToPseudoMeson_cos0_max, this, std::placeholders::_1);
	double res = Optimization::GoldenRatio<double>(func); //-1 to invert function
	return -func(res);
}

double PhaseSpace::ToPseudoMeson_cos0_max(double cos0)
{
	double x = _vars[0];
	double y = _vars[1];
	double cos0_ = -1 + 2 * cos0;

	return -ToPseudoMeson(cos0_, x, y);
}


double PhaseSpace::ToPseudoMeson(double cos0, double x, double y)
{	//either M2_NeutrinoPseudoMeson or M2_LeptonPseudoMeson
	return dGammad2_2B(x, y) * (this->*_M2_F)(cos0, x, y);
}



//neutrino vector meson

double PhaseSpace::nRHO0(const Mixing &mix)
{
	return NeutrinoVectorMeson(Channel::nRHO0, _masspdg[0].first, _masspdg[1].first) / _table[Channel::nRHO0];
}

double PhaseSpace::nOMEGA(const Mixing &mix)
{
	return NeutrinoVectorMeson(Channel::nOMEGA, _masspdg[0].first, _masspdg[1].first) / _table[Channel::nOMEGA];
}

double PhaseSpace::nPHI(const Mixing &mix)
{
	return NeutrinoVectorMeson(Channel::nOMEGA, _masspdg[0].first, _masspdg[1].first) / _table[Channel::nPHI];
}

double PhaseSpace::NeutrinoVectorMeson(Channel::Name chan, double m_neut, double m_meson)
{
	_m_parent = _N.M();
	double dMn2 = std::pow(m_neut / _m_parent, 2);
	double dMM2 = std::pow(m_meson / _m_parent, 2);

	_M2_F = &Amplitude::M2_NeutrinoVectorMeson;

	if (!_table.count(chan))
		_table[chan] = ToVectorMeson_max(dMn2, dMM2);

	double cos0 = Kinematic_2B();

	_M2_F = &Amplitude::M2_NeutrinoVectorMeson;
	return ToVectorMeson(cos0, dMn2, dMM2);
}

double PhaseSpace::ERHO(const Mixing &mix)
{
	return LeptonVectorMeson(Channel::ERHO, _masspdg[0].first, _masspdg[1].first) / _table[Channel::ERHO];
}

double PhaseSpace::MRHO(const Mixing &mix)
{
	return LeptonVectorMeson(Channel::MRHO, _masspdg[0].first, _masspdg[1].first) / _table[Channel::MRHO];
}

double PhaseSpace::EKAx(const Mixing &mix)
{
	return LeptonVectorMeson(Channel::EKAx, _masspdg[0].first, _masspdg[1].first) / _table[Channel::EKAx];
}

double PhaseSpace::MKAx(const Mixing &mix)
{
	return LeptonVectorMeson(Channel::MKAx, _masspdg[0].first, _masspdg[1].first) / _table[Channel::MKAx];
}

double PhaseSpace::LeptonVectorMeson(Channel::Name chan, double m_lepton, double m_meson)
{
	_m_parent = _N.M();
	double dML2 = std::pow(m_lepton / _m_parent, 2);
	double dMM2 = std::pow(m_meson / _m_parent, 2);

	_M2_F = &Amplitude::M2_LeptonVectorMeson;

	if (!_table.count(chan))
		_table[chan] = ToVectorMeson_max(dML2, dMM2);

	double cos0 = Kinematic_2B();

	return ToVectorMeson(cos0, dML2, dMM2);
}

// common for decays into vector meson + lepton

double PhaseSpace::ToVectorMeson_max(double x, double y)
{
	_vars = {x, y};

	auto func = std::bind(&PhaseSpace::ToVectorMeson_cos0_max, this, std::placeholders::_1);
	double res = Optimization::GoldenRatio<double>(func); //-1 to invert function
	return -func(res);
}

double PhaseSpace::ToVectorMeson_cos0_max(double cos0)
{
	double x = _vars[0];
	double y = _vars[1];
	double cos0_ = -1 + 2 * cos0;

	return -ToVectorMeson(cos0_, x, y);
}

double PhaseSpace::ToVectorMeson(double cos0, double x, double y)
{
	return dGammad2_2B(x, y) * (this->*_M2_F)(cos0, x, y);
}


//// PRODUCTION modes

//pure leptonic decays
//
double PhaseSpace::MuonE(const Mixing &mix)
{
	return AntileptonNeutrino(Channel::MuonE, _masspdg[0].first, _masspdg[1].first, _masspdg[2].first) / _table[Channel::MuonE];
}

double PhaseSpace::MuonM(const Mixing &mix)
{
	return LeptonNeutrino(Channel::MuonM, _masspdg[0].first, _masspdg[1].first, _masspdg[2].first) / _table[Channel::MuonM];
}

double PhaseSpace::TauEE(const Mixing &mix)
{
	return AntileptonNeutrino(Channel::TauEE, _masspdg[0].first, _masspdg[1].first, _masspdg[2].first) / _table[Channel::TauEE];
}

double PhaseSpace::TauET(const Mixing &mix)
{
	return LeptonNeutrino(Channel::TauET, _masspdg[0].first, _masspdg[1].first, _masspdg[2].first) / _table[Channel::TauET];
}

double PhaseSpace::TauMM(const Mixing &mix)
{
	return AntileptonNeutrino(Channel::TauMM, _masspdg[0].first, _masspdg[1].first, _masspdg[2].first) / _table[Channel::TauMM];
}

double PhaseSpace::TauMT(const Mixing &mix)
{
	return LeptonNeutrino(Channel::TauMT, _masspdg[0].first, _masspdg[1].first, _masspdg[2].first) / _table[Channel::TauMT];
}

//	p1,u = N, p2,t = L, p3,s = n
double PhaseSpace::LeptonNeutrino(Channel::Name chan, double m_lepton0, double m_lepton, double m_neut)
{
	_m_parent = m_lepton0;
	double dML2 = std::pow(m_lepton / _m_parent, 2);
	double dMn2 = std::pow(m_neut / _m_parent, 2);
	double dMN2 = std::pow(_N.M() / _m_parent, 2);

	if (!_table.count(chan))
		_table[chan] = dGammad5_3B() * LeptonNeutrino_max(dMn2, dML2, dMN2);

	std::array<double, 6> kine = Kinematic_3B();// = s
	return dGammad5_3B() * M2_LeptonNeutrino(kine[2], dMn2, dML2, dMN2);
}

double PhaseSpace::LeptonNeutrino_max(double x, double y, double z)
{
	_vars = {x, y, z};
	
	auto func = std::bind(&PhaseSpace::LeptonNeutrino_u_max, this, std::placeholders::_1);
	double res = Optimization::GoldenRatio<double>(func); //-1 to invert function
	return -func(res);
}

double PhaseSpace::LeptonNeutrino_u_max(double u)	//vars is s 
{
	double x = _vars[0];
	double y = _vars[1];
	double z = _vars[2];

	// Limit changes u
	Limit(u, y, z, x);

	return -M2_LeptonNeutrino(u, x, y, z);
}

double PhaseSpace::AntileptonNeutrino(Channel::Name chan, double m_lepton0, double m_lepton, double m_neut)
{
	_m_parent = m_lepton0;
	double dML2 = std::pow(m_lepton / _m_parent, 2);
	double dMn2 = std::pow(m_neut / _m_parent, 2);
	double dMN2 = std::pow(_N.M() / _m_parent, 2);

	if (!_table.count(chan))
		_table[chan] = dGammad5_3B() * AntileptonNeutrino_max(dMn2, dML2, dMN2);

	std::array<double, 6> kine = Kinematic_3B();// = s
	return dGammad5_3B() * M2_AntileptonNeutrino(kine[0], dMn2, dML2, dMN2);
}

double PhaseSpace::AntileptonNeutrino_max(double x, double y, double z)	//vars is s 
{
	_vars = {x, y, z};

	auto func = std::bind(&PhaseSpace::AntileptonNeutrino_s_max, this, std::placeholders::_1);
	double res = Optimization::GoldenRatio<double>(func); //-1 to invert function
	return -func(res);
}

double PhaseSpace::AntileptonNeutrino_s_max(double s)	//vars is s 
{
	double x = _vars[0];
	double y = _vars[1];
	double z = _vars[2];

	// limit changes s
	Limit(s, x, y, z);

	return -M2_AntileptonNeutrino(s, x, y, z);
}



//lepton decay into meson stuff
//
double PhaseSpace::TauPI(const Mixing &mix)
{
	return LeptonMeson(Channel::TauPI, _masspdg[0].first, _masspdg[1].first) / _table[Channel::TauPI];
}

/*
double PhaseSpace::Tau2PI_ratio()
{
	return LeptonMeson(Channel::Tau2PI, _masspdg[0].first, _masspgs[1].first) / _table[Channel::Tau2PI];
}
*/


double PhaseSpace::LeptonMeson(Channel::Name chan, double m_lepton, double m_meson)
{
	_m_parent = m_lepton;
	double dMN2 = std::pow(_N.M() / _m_parent, 2);
	double dMM2 = std::pow(m_meson / _m_parent, 2);

	if (!_table.count(chan))
		_table[chan] = dGammad2_2B(dMN2, dMM2) * M2_LeptonTwo(dMN2, dMM2);

	return dGammad2_2B(dMN2, dMM2) * M2_LeptonTwo(dMN2, dMM2);
}

//meson two body decay
//
double PhaseSpace::PionE(const Mixing &mix)
{
	return MesonTwo(Channel::PionE, _masspdg[0].first, _masspdg[1].first) / _table[Channel::PionE];
}

double PhaseSpace::PionM(const Mixing &mix)
{
	return MesonTwo(Channel::PionM, _masspdg[0].first, _masspdg[1].first) / _table[Channel::PionM];
}

double PhaseSpace::KaonE(const Mixing &mix)
{
	return MesonTwo(Channel::KaonE, _masspdg[0].first, _masspdg[1].first) / _table[Channel::KaonE];
}

double PhaseSpace::KaonM(const Mixing &mix)
{
	return MesonTwo(Channel::KaonM, _masspdg[0].first, _masspdg[1].first) / _table[Channel::KaonM];
}

double PhaseSpace::CharmE(const Mixing &mix)
{
	return MesonTwo(Channel::CharmE, _masspdg[0].first, _masspdg[1].first) / _table[Channel::CharmE];
}

double PhaseSpace::CharmM(const Mixing &mix)
{
	return MesonTwo(Channel::CharmM, _masspdg[0].first, _masspdg[1].first) / _table[Channel::CharmM];
}

double PhaseSpace::CharmT(const Mixing &mix)
{
	return MesonTwo(Channel::CharmT, _masspdg[0].first, _masspdg[1].first) / _table[Channel::CharmT];
}


// is the ratio just 1?
double PhaseSpace::MesonTwo(Channel::Name chan, double m_meson, double m_lepton)
{
	_m_parent = m_meson;
	double dMN2 = std::pow(_N.M() / _m_parent, 2);
	double dML2 = std::pow(m_lepton / _m_parent, 2);

	if (!_table.count(chan))
		_table[chan] = dGammad2_2B(dMN2, dML2) * M2_MesonTwo(dMN2, dML2);

	return dGammad2_2B(dMN2, dML2) * M2_MesonTwo(dMN2, dML2);
}




//three body decays of meson
//
double PhaseSpace::Kaon0E(const Mixing &mix)
{
	return MesonThree(Channel::Kaon0E, _masspdg[0].first, _masspdg[1].first, _masspdg[2].first,
				Const::KCL_, Const::KCL0) / _table[Channel::Kaon0E];
}

double PhaseSpace::Kaon0M(const Mixing &mix)
{
	return MesonThree(Channel::Kaon0M, _masspdg[0].first, _masspdg[1].first, _masspdg[2].first,
				Const::KCL_, Const::KCL0) / _table[Channel::Kaon0M];
}

double PhaseSpace::KaonCE(const Mixing &mix)
{
	return MesonThree(Channel::KaonCE, _masspdg[0].first, _masspdg[1].first, _masspdg[2].first,
				Const::KCL_, Const::KCL0) / _table[Channel::KaonCE];
}

double PhaseSpace::KaonCM(const Mixing &mix)
{
	return MesonThree(Channel::KaonCM, _masspdg[0].first, _masspdg[1].first, _masspdg[2].first,
				Const::KCL_, Const::KCL0) / _table[Channel::KaonCM];
}


double PhaseSpace::MesonThree(Channel::Name chan, double m_meson0, double m_meson, double m_lepton,
				double L_, double L0)	//decay constant not important
{
	_m_parent = m_meson0;
	double dMM2 = std::pow(m_meson / _m_parent, 2);
	double dML2 = std::pow(m_lepton / _m_parent, 2);
	double dMN2 = std::pow(_N.M() / _m_parent, 2);	

	if (!_table.count(chan))
		_table[chan] = dGammad5_3B() * MesonThree_max(dMM2, dML2, dMN2, L_, L0);

	std::array<double, 6> kine = Kinematic_3B();	 // = s    = t
	return dGammad5_3B() * M2_MesonThree(kine[0], kine[1], dMM2, dML2, dMN2, L_, L0) ;
}

double PhaseSpace::MesonThree_max(double x, double y, double z, double L_, double L0)	//no angle dependence
{
	_vars = {x, y, z, L_, L0};

	auto func = std::bind(&PhaseSpace::F_MesonThree_max, this, std::placeholders::_1);
	auto res = Optimization::NelderMead<double>(func, 2); //-1 to invert function
	return -func(&res[0]);
}

// returns a reversed sign for minimization algorithm
double PhaseSpace::F_MesonThree_max(double p[])
{
	double x  = _vars[0];
	double y  = _vars[1];
	double z  = _vars[2];
	double gL = _vars[3];
	double gR = _vars[4];

	double s_ = p[0];
	double t_ = p[1];

	return -M2_MesonThree(s_, t_, x, y, z, gL, gR);
}





///KINEMATICS/////
/////obtaining kinematics from decay products
//
double PhaseSpace::Kinematic_2B()
{
	if (_masspdg.size() > 2)
		throw std::logic_error("Kinematic_2B error");

	TLorentzVector vec1 = *(_genps->GetDecay(1));

	return std::cos(vec1.Theta());
}

//the mandelstam variables refer to a specific vector
//s -> p3
//t -> p2
//u -> p1
std::array<double, 6> PhaseSpace::Kinematic_3B()
{
	if (_masspdg.size() > 3)
		throw std::logic_error("Kinematic_3B error");

	TLorentzVector vec1 = *(_genps->GetDecay(0));
	TLorentzVector vec2 = *(_genps->GetDecay(1));
	TLorentzVector vec3 = *(_genps->GetDecay(2));

	TLorentzVector rest(0, 0, 0, _m_parent);
	TLorentzVector vec_u = rest - vec1;
	TLorentzVector vec_t = rest - vec2;
	TLorentzVector vec_s = rest - vec3;

	std::array<double, 6> ret;

	ret[0] = vec_u.M2() / std::pow(_m_parent, 2); // s
	ret[1] = vec_t.M2() / std::pow(_m_parent, 2); // t
	ret[2] = vec_s.M2() / std::pow(_m_parent, 2); // u

	ret[3] = std::cos(vec1.Theta());	// cos0u 
	ret[4] = std::cos(vec2.Theta());	// cos0t
	ret[5] = std::cos(vec3.Theta());	// cos0s

	return ret;
}
