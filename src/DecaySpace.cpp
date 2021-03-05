#include "physics/DecaySpace.h"

std::pair<std::vector<Particle>, double> DecaySpace::Generate(Decay::Channel chan,
				const TLorentzVector &frame, const Mixing &mix) {
	std::vector<Particle> parts;
	if (!SetDecay(chan))	// decay not allowed
		return std::make_pair(parts, 0.); 	// return empty vector

	size_t cc = 0;
	std::uniform_real_distribution<> rndm;
	double weight = _genps->Generate();	// generate phase space
	double val = Gamma(chan, mix); 
	if (val <= 1e-12 || val - 1e-12 > 1.)	// phase space impossible
		return std::make_pair(parts, val);

	while (rndm(RNG::_mt) > val && cc++ < 1e4) {
		weight = _genps->Generate();	// generate phase space
		val = Gamma(chan, mix);
	}

	if (cc > 1e4)
		return std::make_pair(parts, val);

	// create parent frame
	//TLorentzVector vec = *(_genps->GetDecay(0));
	//vec.Boost(frame.BoostVector());

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

	//parts.push_back(Particle(sign * pdgs.front(), vec));
	for (size_t i = 0; i < pdgs.size(); ++i) {
		TLorentzVector vec = *(_genps->GetDecay(i));
		vec.Boost(frame.BoostVector());
		parts.emplace_back(sign * pdgs[i], vec);
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

//Neutrino LeptonLepton AA -- phase space depends on the mixing!

double DecaySpace::nEE(const Mixing &mix)
{
	return NeutrinoLeptonAA(Decay::Channel::nEE, Const::MNeutrino, Const::MElectron,
			mix.Ue(2), mix.Um(2) + mix.Ut(2)) / _table[Decay::Channel::nEE];
}

double DecaySpace::nMM(const Mixing &mix)
{
	return NeutrinoLeptonAA(Decay::Channel::nMM, Const::MNeutrino, Const::MMuon,
			mix.Um(2), mix.Ue(2) + mix.Ut(2)) / _table[Decay::Channel::nMM];
}


double DecaySpace::NeutrinoLeptonAA(Decay::Channel chan, double m_neut, double m_lepton,
				    double uu, double uo)
{
	_m_parent = _N.M();
	double x = std::pow(m_neut / _m_parent, 2);
	double y = std::pow(m_lepton / _m_parent, 2);

	//neutrino flavour is the same as leptons -> both Z and W interaction
	double gLcc = +0.5 + Const::sin2W;	//times U(lepton flavour)
	double gRcc = Const::sin2W;

	//neutrino flavour is different from leptons -> Z interaction
	double gLnc = -0.5 + Const::sin2W;	//times U(neutrino flavour)
	double gRnc = Const::sin2W;

	double mixWW = uu * (gLcc * gLcc + gRcc * gRcc)
		     + uo * (gLnc * gLnc + gRnc * gRnc);
	double mixWZ = uu * (2 * gLcc * gRcc) + uo * (2 * gLnc * gRnc);


	// compute max first, if needed
	if (!_table.count(chan)) {
		auto func = [&](double p[]) {	    // s     u     cos0  cos1
			return - NeutrinoLeptonAA_max(p[0], p[1], p[2], p[3]); };

		_vars = {x, y, mixWW, mixWZ};
		auto max_ = Optimization::NelderMead<double>(func, 4);

		_table[chan] = - dGammad5_3B() * func(&max_[0]);
	}

	std::array<double, 6> k = Kinematic3B(); // s      u       cos0s     cos0u
	return dGammad5_3B() * (mixWW * M2_WW(k[0], k[2], x, y, y)
			      + mixWZ * M2_WZ(k[1], k[3], x, y, y) );
}


// returns a reversed sign for minimization algorithm
double DecaySpace::NeutrinoLeptonAA_max(double s, double u, double cos0, double cos1)
{
	if (std::abs(cos0) > 1 || std::abs(cos1) > 1)
		return 0.;
	if (s < 0 || u < 0)
		return 0.;

	const double &x  = _vars[0];
	const double &y  = _vars[1];

	if (s + u > 1 + x + 2 * y)
		return 0.;

	const double &ww = _vars[2];
	const double &wz = _vars[3];

	return ww * M2_WW(s, cos0, x, y, y) + wz * M2_WZ(u, cos1, x, y, y);
}


//Neutrino LeptonLepton AB -- depends on mixing

double DecaySpace::nEM(const Mixing &mix)
{
	return NeutrinoLeptonAB(Decay::Channel::nEM, Const::MNeutrino,
			        Const::MElectron, Const::MMuon, mix.Ue(2), mix.Um(2))
			/ _table[Decay::Channel::nEM];
}

double DecaySpace::nET(const Mixing &mix)
{
	return NeutrinoLeptonAB(Decay::Channel::nET, Const::MNeutrino,
			        Const::MElectron, Const::MTau, mix.Ue(2), mix.Ut(2))
			/ _table[Decay::Channel::nET];
}

double DecaySpace::nMT(const Mixing &mix)
{
	return NeutrinoLeptonAB(Decay::Channel::nMT, Const::MNeutrino,
			        Const::MMuon, Const::MTau, mix.Um(2), mix.Ut(2))
			/ _table[Decay::Channel::nMT];
}


double DecaySpace::NeutrinoLeptonAB(Decay::Channel chan,
		double m_neut, double m_leptonA, double m_leptonB,
		double uu, double uo)
{
	_m_parent = _N.M();
	double x = std::pow(m_neut / _m_parent, 2);
	double y = std::pow(m_leptonA / _m_parent, 2);
	double z = std::pow(m_leptonB / _m_parent, 2);

	//double gL = 1.0;
	//double gR = 0.0;

	if (!_table.count(chan)) {
		// compute maximum of just one component
		auto func = [&](double p[]) {
			return - NeutrinoLeptonAB_max(p[0], p[1], p[2], p[3]); };

		_vars = {x, z, y, uu, uo};
		auto max_ = Optimization::NelderMead<double>(func, 4);
		_table[chan] = - dGammad5_3B() * func(&max_[0]);
		//NeutrinoLeptonLepton_max(x, z, y, 1., 0.);
	}

	// combine two components here with mixing
	std::array<double, 6> k = Kinematic3B();// s    u    cos0s   cos0u
	return dGammad5_3B() * (uu * M2_WW(k[0], k[3], x, y, z)
			      + uo * M2_WW(k[2], k[5], x, y, z) );
}

// returns a reversed sign for minimization algorithm
double DecaySpace::NeutrinoLeptonAB_max(double s, double u, double cos0, double cos1)
{
	if (std::abs(cos0) > 1 || std::abs(cos1) > 1)
		return 0.;
	if (s < 0 || u < 0)
		return 0.;

	const double &x  = _vars[0];
	const double &y  = _vars[1];
	const double &z  = _vars[2];

	if (s + u > 1 + x + y + z)
		return 0.;

	const double &mu = _vars[3];
	const double &mo = _vars[4];

	return (mu + mo) * M2_WW(s, cos0, x, y, z); // + mo * M2_WW(u, cos1, x, y, z);
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

		auto func = [&](double cos0) {
			return - NeutrinoPseudoMeson_max(cos0); };
		double max_ = Optimization::GoldenRatio<double>(func);
		_table[chan] = - dGammad2_2B(x, y) * func(max_);
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
	return M2_NeutrinoPseudoMeson(cos0, x, y);
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

		auto func = [&](double cos0) {
			return - LeptonPseudoMeson_max(cos0); };
		double max_ = Optimization::GoldenRatio<double>(func);
		_table[chan] = - dGammad2_2B(x, y) * func(max_);
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

	return M2_LeptonPseudoMeson(cos0, x, y);
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

		auto func = [&](double cos0) {
			return - NeutrinoVectorMeson_max(cos0); };
		double max_ = Optimization::GoldenRatio<double>(func);
		_table[chan] = - dGammad2_2B(x, y) * func(max_);
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
	return M2_NeutrinoVectorMeson(cos0, x, y);
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

		auto func = [&](double cos0) {
			return - LeptonVectorMeson_max(cos0); };
		double max_ = Optimization::GoldenRatio<double>(func);
		_table[chan] = - dGammad2_2B(x, y) * func(max_);
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

	return M2_LeptonVectorMeson(cos0, x, y);
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
