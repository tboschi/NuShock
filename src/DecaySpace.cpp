#include "physics/DecaySpace.h"

std::pair<std::vector<Particle>, double> DecaySpace::Generate(Decay::Channel chan,
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
	auto pdgs = Decay::Pdgs(chan);
	parts.reserve(pdgs.size());

	// reverse pdg sign if particle is "anti"
	int sign = 1;
	if (_N.IsMajorana()) {
		std::bernoulli_distribution b(0.5);
		sign = b(RNG::_mt) ? 1 : -1;
	}
	else // is dirac
		sign = _N.IsParticle() ? 1 : -1;

	parts.push_back(Particle(sign * pdgs.front(), vec));
	for (size_t i = 1; i < pdgs.size(); ++i) {
		TLorentzVector vec = *(_genps->GetDecay(i));
		vec.Boost(frame.BoostVector());
		parts.push_back(Particle(sign * pdgs[i], vec));
	}

	return std::make_pair(parts, weight);
}

bool DecaySpace::SetDecay(Decay::Channel chan)
{
	auto mass = Decay::Masses(chan);
	double *ps_masses = &mass[0];

	// compute decay in rest frame, boost later
	TLorentzVector rest(0, 0, 0, _N.M());
	return _genps->SetDecay(rest, mass.size(), ps_masses);
}

double DecaySpace::Gamma(Decay::Channel chan, const Mixing &mix)
{
	using GammaR = double (DecaySpace::*)(const Mixing &mix);
	GammaR gr;

	switch (chan)
	{
		//case Decay::Channel::ALL:
		case Decay::Channel::nnn:
		case Decay::Channel::nGAMMA:
			return 1.;
		case Decay::Channel::nEE:
			gr = &DecaySpace::nEE;
			break;
		case Decay::Channel::nEM:
			gr = &DecaySpace::nEM;
			break;
		case Decay::Channel::nMM:
			gr = &DecaySpace::nMM;
			break;
		case Decay::Channel::nET:
			gr = &DecaySpace::nET;
			break;
		case Decay::Channel::nMT:
			gr = &DecaySpace::nMT;
			break;
		case Decay::Channel::nPI0:
			gr = &DecaySpace::nPI0;
			break;
		case Decay::Channel::EPI:
			gr = &DecaySpace::EPI;
			break;
		case Decay::Channel::MPI:
			gr = &DecaySpace::MPI;
			break;
		case Decay::Channel::TPI:
			gr = &DecaySpace::TPI;
			break;
		case Decay::Channel::EKA:
			gr = &DecaySpace::EKA;
			break;
		case Decay::Channel::MKA:
			gr = &DecaySpace::MKA;
			break;
		case Decay::Channel::nRHO0:
			gr = &DecaySpace::nRHO0;
			break;
		case Decay::Channel::ERHO:
			gr = &DecaySpace::ERHO;
			break;
		case Decay::Channel::MRHO:
			gr = &DecaySpace::MRHO;
			break;
		case Decay::Channel::EKAx:
			gr = &DecaySpace::EKAx;
			break;
		case Decay::Channel::MKAx:
			gr = &DecaySpace::MKAx;
			break;
		case Decay::Channel::nOMEGA:
			gr = &DecaySpace::nOMEGA;
			break;
		case Decay::Channel::nETA:
			gr = &DecaySpace::nETA;
			break;
		case Decay::Channel::nETAi:
			gr = &DecaySpace::nETAi;
			break;
		case Decay::Channel::nPHI:
			gr = &DecaySpace::nPHI;
			break;
		case Decay::Channel::EDs:
			gr = &DecaySpace::EDs;
			break;
		default:
			throw std::invalid_argument("Decay channel "
						+ Decay::toString(chan) + " is unknown");
	}
	
	return (this->*gr)(mix);
}

//Neutrino LeptonLepton AA

double DecaySpace::nEE(const Mixing &mix)
{
	auto gamma = NeutrinoLeptonAA(Decay::Channel::nEE, Decay::Channel::nEE_o,
				Const::MNeutrino, Const::MElectron);

	if (gamma.first < 0.)
		gamma.first = 0.;
	if (gamma.second < 0.)
		gamma.second = 0.;

	return ( gamma.first * mix.Ue(2) + gamma.second * (mix.Um(2) + mix.Ut(2)))
	     / ( _table[Decay::Channel::nEE] * mix.Ue(2)
	       + _table[Decay::Channel::nEE_o] * (mix.Um(2) + mix.Ut(2)) );
}

double DecaySpace::nMM(const Mixing &mix)
{
	auto gamma = NeutrinoLeptonAA(Decay::Channel::nMM, Decay::Channel::nMM_o,
					Const::MNeutrino, Const::MMuon);

	if (gamma.first < 0.)
		gamma.first = 0.;

	if (gamma.second < 0.)
		gamma.second = 0.;

	return ( gamma.first * mix.Um(2) + gamma.second * (mix.Ue(2) + mix.Ut(2)))
	     / ( _table[Decay::Channel::nMM] * mix.Um(2)
	       + _table[Decay::Channel::nMM_o] * (mix.Ue(2) + mix.Ut(2)));
}


std::pair<double, double> DecaySpace::NeutrinoLeptonAA(Decay::Channel chan, Decay::Channel chan_o, double m_neut, double m_lepton)
{
	_m_parent = _N.M();
	double x = std::pow(m_neut / _m_parent, 2);
	double y = std::pow(m_lepton / _m_parent, 2);

	//neutrino flavour is the same as leptons -> both Z and W interaction
	double gLCC = +0.5 + Const::sin2W;	//times U(lepton flavour)
	double gRCC = Const::sin2W;

	//neutrino flavour is different from leptons -> Z interaction
	double gLNC = -0.5 + Const::sin2W;	//times U(neutrino flavour)
	double gRNC = Const::sin2W;

	// compute max first, if needed
	if (!_table.count(chan)) {
		auto func = [&](double p[]) { return NeutrinoLeptonAA_max(p); };

		_vars = {x, y, gLCC, gRCC};
		auto max_ = Optimization::NelderMead<double>(func, 4); //-1 to invert function

		_vars = {x, y, gLNC, gRNC};
		auto max_o = Optimization::NelderMead<double>(func, 4); //-1 to invert function
		_table[chan]   = dGammad5_3B() * (-func(&max_[0])); // NeutrinoLeptonLepton_max(x, y, y, gL_CC, gR_CC);
		_table[chan_o] = dGammad5_3B() * (-func(&max_o[0])); // NeutrinoLeptonLepton_max(x, y, y, gL_NC, gR_NC);
	}

	//return std::make_pair(NeutrinoLeptonLepton(kine[0], kine[2], kine[3], kine[5], dMn2, dML2, dML2, gL_CC, gR_CC),
			      //NeutrinoLeptonLepton(kine[0], kine[2], kine[3], kine[5], dMn2, dML2, dML2, gL_NC, gR_NC));
	std::array<double, 6> k = Kinematic3B(); // s      u       cos0s     cos0u
	double m2_ww = M2_WW(k[0], k[3], x, y, y);
	double m2_wz = M2_WZ(k[2], k[5], x, y, y);
	return std::make_pair(dGammad5_3B() * (gLCC * gLCC + gRCC * gRCC) * m2_ww
						      + (2 * gLCC * gRCC) * m2_wz,
			      dGammad5_3B() * (gLNC * gLNC + gRNC * gRNC) * m2_ww
						      + (2 * gLNC * gRNC) * m2_wz);
}


// returns a reversed sign for minimization algorithm
double DecaySpace::NeutrinoLeptonAA_max(double p[])
{
	double x  = _vars[0];
	double y  = _vars[1];
	double gL = _vars[2];
	double gR = _vars[3];

	double s   = p[0];
	double u   = p[1];
	double cos0 = p[2];
	double cos1 = p[3];

	// Limits...?

	if (std::abs(cos0) <= 1 && std::abs(cos1) <= 1)
		return (gL * gL + gR * gR) * M2_WW(s, cos0, x, y, y)
			   + (2 * gL * gR) * M2_WZ(u, cos1, x, y, y);
	return 0.;
}


//Neutrino LeptonLepton AB

double DecaySpace::nEM(const Mixing &mix)
{
	auto gamma = NeutrinoLeptonAB(Decay::Channel::nEM, Decay::Channel::nEM_o,
					Const::MNeutrino, Const::MElectron, Const::MMuon);

	if (gamma.first < 0.)
		gamma.first = 0.;
	if (gamma.second < 0.)
		gamma.second = 0.;

	//return (maxName_u > 1e-25 || maxName_o > 1e-25) ? 
	return ( gamma.first * mix.Ue(2) + gamma.second * mix.Um(2))
	     / ( _table[Decay::Channel::nEM] * mix.Ue(2)
	       + _table[Decay::Channel::nEM_o] * mix.Um(2));
}

double DecaySpace::nET(const Mixing &mix)
{
	auto gamma = NeutrinoLeptonAB(Decay::Channel::nET, Decay::Channel::nET_o,
				Const::MNeutrino, Const::MElectron, Const::MTau);

	if (gamma.first < 0.)
		gamma.first = 0.;
	if (gamma.second < 0.)
		gamma.second = 0.;

	//return (maxName_u > 1e-25 || maxName_o > 1e-25) ? 
	return ( gamma.first * mix.Ue(2) + gamma.second * mix.Ut(2))
	     / ( _table[Decay::Channel::nET] * mix.Ue(2)
	       + _table[Decay::Channel::nET_o] * mix.Ut(2));
}

double DecaySpace::nMT(const Mixing &mix)
{
	// computes max too
	auto gamma = NeutrinoLeptonAB(Decay::Channel::nMT, Decay::Channel::nMT_o,
				Const::MNeutrino, Const::MMuon, Const::MTau);

	if (gamma.first < 0.)
		gamma.first = 0.;
	if (gamma.second < 0.)
		gamma.second = 0.;

	//return (maxName_u > 1e-25 || maxName_o > 1e-25) ? 
	return (gamma.first * mix.Um(2) + gamma.second * mix.Ut(2))
	     / ( _table[Decay::Channel::nMT] * mix.Um(2)
	       + _table[Decay::Channel::nMT_o] * mix.Ut(2));
}


std::pair<double, double> DecaySpace::NeutrinoLeptonAB(Decay::Channel chan, Decay::Channel chan_o, double m_neut, double m_leptonA, double m_leptonB)
{
	_m_parent = _N.M();
	double x = std::pow(m_neut / _m_parent, 2);
	double y = std::pow(m_leptonA / _m_parent, 2);
	double z = std::pow(m_leptonB / _m_parent, 2);

	//double gL = 1.0;
	//double gR = 0.0;

	if (!_table.count(chan)) {
		auto func = [&](double p[]) { return NeutrinoLeptonAB_max(p); };

		_vars = {x, z, y};
		auto max_  = Optimization::NelderMead<double>(func, 2); //-1 to invert function

		_vars = {x, z, y};
		auto max_o = Optimization::NelderMead<double>(func, 2); //-1 to invert function
		_table[chan]   = dGammad5_3B() * (-func(&max_[0])); //NeutrinoLeptonLepton_max(x, z, y, 1., 0.);
		_table[chan_o] = dGammad5_3B() * (-func(&max_o[0])); //NeutrinoLeptonLepton_max(x, y, z, 1., 0.);
	}

	std::array<double, 6> k = Kinematic3B();		// s      u       cos0s     cos0u
	//return std::make_pair(dGammad5_3B() * NeutrinoLeptonLepton(k[0], k[2], k[3], k[5], x, y, z, 1., 0.),
			      //dGammad5_3B() * NeutrinoLeptonLepton(k[2], k[0], k[5], k[3], x, y, z, 1., 0.));
	return std::make_pair(dGammad5_3B() * M2_WW(k[0], k[3], x, y, z), // passing s and cos0s
			      dGammad5_3B() * M2_WW(k[2], k[5], x, y, z)); // passing u and cos0u
}

// returns a reversed sign for minimization algorithm
double DecaySpace::NeutrinoLeptonAB_max(double p[])
{
	double x  = _vars[0];
	double y  = _vars[1];
	double z  = _vars[2];
	//double gL = _vars[3];	 == 1
	//double gR = _vars[4];	 == 0

	double s   = p[0];
	//double u   = p[1];
	double cos0 = p[1];
	//double cos1 = p[3];
	// LIMIT??

	if (std::abs(cos0) <= 1) // && std::abs(cos1) <= 1)
		return M2_WW(s, cos0, x, y, z);
	return 0.;
}

// common max functions

/*
double DecaySpace::NeutrinoLeptonLepton_max(double x, double y, double z, double gL, double gR)
{
	_vars = {x, y, z, gL, gR};

	auto func = [&](double p[]) { return F_NeutrinoLeptonLepton_max(p); };
	auto max_ = Optimization::NelderMead<double>(func, 4); //-1 to invert function
	return -func(&max_[0]);
}
*/

/*
double DecaySpace::NeutrinoLeptonLepton(double s, double u, double cos0, double cos1,
					double x, double y, double z,
					double gL, double gR)
{
	return dGammad5_3B() * (gL * gL + gR * gR) * M2_WW(s, cos0, x, y, z) +
		                     (2 * gL * gR) * M2_WZ(u, cos1, x, y, z);
}
*/




//neutrino psuedomeson

double DecaySpace::nPI0(const Mixing &mix)
{
	return NeutrinoPseudoMeson(Decay::Channel::nPI0, Const::MNeutrino, Const::MPion0)
		/ _table[Decay::Channel::nPI0];
}

double DecaySpace::nETA(const Mixing &mix)
{
	return NeutrinoPseudoMeson(Decay::Channel::nETA, Const::MNeutrino, Const::MEta)
		/ _table[Decay::Channel::nETA];
}

double DecaySpace::nETAi(const Mixing &mix)
{
	return NeutrinoPseudoMeson(Decay::Channel::nETAi, Const::MNeutrino, Const::MEtai)
		/ _table[Decay::Channel::nETAi];
}

//
double DecaySpace::NeutrinoPseudoMeson(Decay::Channel chan, double m_neut, double m_meson)
{
	_m_parent = _N.M();
	double x = std::pow(m_neut / _m_parent, 2);
	double y = std::pow(m_meson / _m_parent, 2);

	if (!_table.count(chan)) {
		_vars = {x, y};

		auto func = [&](double cos0) { return NeutrinoPseudoMeson_max(cos0); };
		double max_ = Optimization::GoldenRatio<double>(func); //-1 to invert function
		_table[chan] = dGammad2_2B(x, y) * (-func(max_));
	}

	double cos0 = Kinematic2B();

	return dGammad2_2B(x, y) * M2_NeutrinoPseudoMeson(cos0, x, y);
	//return ToPseudoMeson(cos0_, dMn2, dMM2);
}	//     2 is factor from decay constant which is sqrt(2) wrt to charged meson


double DecaySpace::NeutrinoPseudoMeson_max(double cos0)
{
	double x = _vars[0];
	double y = _vars[1];
	cos0 = -1 + 2 * cos0;

	//return -ToPseudoMeson(cos0_, x, y);
	return - M2_NeutrinoPseudoMeson(cos0, x, y);
}


//lepton psuedomeson

double DecaySpace::EPI(const Mixing &mix)
{
	return LeptonPseudoMeson(Decay::Channel::EPI, Const::MElectron, Const::MPion)
		/ _table[Decay::Channel::EPI];
}

double DecaySpace::MPI(const Mixing &mix)
{
	return LeptonPseudoMeson(Decay::Channel::MPI, Const::MMuon, Const::MPion)
		/ _table[Decay::Channel::MPI];
}

double DecaySpace::TPI(const Mixing &mix)
{
	return LeptonPseudoMeson(Decay::Channel::TPI, Const::MTau, Const::MPion)
		/ _table[Decay::Channel::TPI];
}

double DecaySpace::EKA(const Mixing &mix)
{
	return LeptonPseudoMeson(Decay::Channel::EKA, Const::MElectron, Const::MKaon)
		/ _table[Decay::Channel::EKA];
}

double DecaySpace::MKA(const Mixing &mix)
{
	return LeptonPseudoMeson(Decay::Channel::MKA, Const::MMuon, Const::MKaon)
		/ _table[Decay::Channel::MKA];
}

double DecaySpace::EDs(const Mixing &mix)
{
	return LeptonPseudoMeson(Decay::Channel::EDs, Const::MElectron, Const::MDs)
		/ _table[Decay::Channel::EDs];
}

double DecaySpace::LeptonPseudoMeson(Decay::Channel chan, double m_lepton, double m_meson)
{
	_m_parent = _N.M();
	double x = std::pow(m_lepton / _m_parent, 2);
	double y = std::pow(m_meson / _m_parent, 2);

	if (!_table.count(chan)) {
		_vars = {x, y};

		auto func = [&](double cos0) { return LeptonPseudoMeson_max(cos0); };
		double max_ = Optimization::GoldenRatio<double>(func); //-1 to invert function
		_table[chan] = dGammad2_2B(x, y) * (-func(max_));
	}

	double cos0 = Kinematic2B();

	return dGammad2_2B(x, y) * M2_LeptonPseudoMeson(cos0, x, y);
	//return ToPseudoMeson(cos0_, dML2, dMM2);
}


// common for decays into pseudomeson + lepton

/*
double DecaySpace::ToPseudoMeson_max(double x, double y)
{
	_vars = {x, y};

	auto func = std::bind(&DecaySpace::ToPseudoMeson_cos0_max, this, std::placeholders::_1);
	double res = Optimization::GoldenRatio<double>(func); //-1 to invert function
	return -func(res);
}
*/

double DecaySpace::LeptonPseudoMeson_max(double cos0)
{
	double x = _vars[0];
	double y = _vars[1];
	cos0 = -1 + 2 * cos0;

	//return -ToPseudoMeson(cos0_, x, y);
	return - M2_LeptonPseudoMeson(cos0, x, y);
}


/*
double DecaySpace::ToPseudoMeson(double cos0, double x, double y)
{	//either M2_NeutrinoPseudoMeson or M2_LeptonPseudoMeson
	return dGammad2_2B(x, y) * (this->*_M2_F)(cos0, x, y);
}
*/



//neutrino vector meson

double DecaySpace::nRHO0(const Mixing &mix)
{
	return NeutrinoVectorMeson(Decay::Channel::nRHO0, Const::MNeutrino, Const::MRho)
		/ _table[Decay::Channel::nRHO0];
}

double DecaySpace::nOMEGA(const Mixing &mix)
{
	return NeutrinoVectorMeson(Decay::Channel::nOMEGA, Const::MNeutrino, Const::MOmega)
		/ _table[Decay::Channel::nOMEGA];
}

double DecaySpace::nPHI(const Mixing &mix)
{
	return NeutrinoVectorMeson(Decay::Channel::nOMEGA, Const::MNeutrino, Const::MPhi)
		/ _table[Decay::Channel::nPHI];
}

double DecaySpace::NeutrinoVectorMeson(Decay::Channel chan, double m_neut, double m_meson)
{
	_m_parent = _N.M();
	double x = std::pow(m_neut / _m_parent, 2);
	double y = std::pow(m_meson / _m_parent, 2);

	if (!_table.count(chan)) {
		// find max
		_vars = {x, y};

		auto func = [&](double cos0) { return NeutrinoVectorMeson_max(cos0); };
		double max_ = Optimization::GoldenRatio<double>(func); //-1 to invert function
		_table[chan] = dGammad2_2B(x, y) * (-func(max_));
	}

	double cos0 = Kinematic2B();

	return dGammad2_2B(x, y) * M2_NeutrinoVectorMeson(cos0, x, y);
	//return ToVectorMeson(cos0, dMn2, dMM2);
}

double DecaySpace::NeutrinoVectorMeson_max(double cos0)
{
	double x = _vars[0];
	double y = _vars[1];
	cos0 = -1 + 2 * cos0;

	//return -ToVectorMeson(cos0_, x, y);
	return - M2_NeutrinoVectorMeson(cos0, x, y);
}


double DecaySpace::ERHO(const Mixing &mix)
{
	return LeptonVectorMeson(Decay::Channel::ERHO, Const::MElectron, Const::MRho)
		/ _table[Decay::Channel::ERHO];
}

double DecaySpace::MRHO(const Mixing &mix)
{
	return LeptonVectorMeson(Decay::Channel::MRHO, Const::MMuon, Const::MRho)
		/ _table[Decay::Channel::MRHO];
}

double DecaySpace::EKAx(const Mixing &mix)
{
	return LeptonVectorMeson(Decay::Channel::EKAx, Const::MElectron, Const::MKaonx)
		/ _table[Decay::Channel::EKAx];
}

double DecaySpace::MKAx(const Mixing &mix)
{
	return LeptonVectorMeson(Decay::Channel::MKAx, Const::MMuon, Const::MKaonx)
		/ _table[Decay::Channel::MKAx];
}

double DecaySpace::LeptonVectorMeson(Decay::Channel chan, double m_lepton, double m_meson)
{
	_m_parent = _N.M();
	double x = std::pow(m_lepton / _m_parent, 2);
	double y = std::pow(m_meson / _m_parent, 2);

	if (!_table.count(chan)) {
		// find max
		_vars = {x, y};

		auto func = [&](double cos0) { return LeptonVectorMeson_max(cos0); };
		double max_ = Optimization::GoldenRatio<double>(func); //-1 to invert function
		_table[chan] = dGammad2_2B(x, y) * (-func(max_));
	}

	double cos0 = Kinematic2B();

	return dGammad2_2B(x, y) * M2_LeptonVectorMeson(cos0, x, y);
	//return ToVectorMeson(cos0, dML2, dMM2);
}

// common for decays into vector meson + lepton

/*
double DecaySpace::ToVectorMeson_max(double x, double y)
{
	_vars = {x, y};

	auto func = std::bind(&DecaySpace::ToVectorMeson_cos0_max, this, std::placeholders::_1);
	double res = Optimization::GoldenRatio<double>(func); //-1 to invert function
	return -func(res);
}
*/

double DecaySpace::LeptonVectorMeson_max(double cos0)
{
	double x = _vars[0];
	double y = _vars[1];
	cos0 = -1 + 2 * cos0;

	//return -ToVectorMeson(cos0_, x, y);
	return - M2_LeptonVectorMeson(cos0, x, y);
}

/*
double DecaySpace::ToVectorMeson(double cos0, double x, double y)
{
	return dGammad2_2B(x, y) * (this->*_M2_F)(cos0, x, y);
}
*/

void DecaySpace::Reset() {
	_table.clear();
}
