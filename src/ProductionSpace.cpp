#include "physics/ProductionSpace.h"

std::pair<std::vector<Particle>, double> ProductionSpace::Generate(Production::Channel chan,
				const TLorentzVector &frame, const Mixing &mix) {
	std::vector<Particle> parts;
	if (!SetProduction(chan)) // decay not allowed
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
	auto pdgs = Production::Pdgs(chan);
	parts.reserve(pdgs.size());

	parts.push_back(Particle(_N.Pdg(), vec));

	for (size_t i = 1; i < pdgs.size(); ++i) {
		TLorentzVector vec = *(_genps->GetDecay(i));
		vec.Boost(frame.BoostVector());
		parts.push_back(Particle(pdgs[i], vec));
	}

	return std::make_pair(parts, weight);
}


bool ProductionSpace::SetProduction(Production::Channel chan)
{
	auto mass = Production::Masses(chan);
	double *ps_masses = &mass[0];
	ps_masses[0] = _N.M();

	// compute decay in rest frame, boost later
	TLorentzVector rest(0, 0, 0, mass.front());
	return _genps->SetDecay(rest, mass.size(), ps_masses);
}

double ProductionSpace::Gamma(Production::Channel chan, const Mixing &mix)
{
	using GammaR = double (ProductionSpace::*)(const Mixing &mix);
	GammaR gr;

	switch (chan)
	{
		case Production::Channel::MuonE:	//these are needed as well
			gr = &ProductionSpace::MuonE;
			break;
			return 1.;
		case Production::Channel::MuonM:
			gr = &ProductionSpace::MuonM;
			break;
		case Production::Channel::TauEE:
			gr = &ProductionSpace::TauEE;
			break;
		case Production::Channel::TauET:
			gr = &ProductionSpace::TauET;
			break;
		case Production::Channel::TauMM:
			gr = &ProductionSpace::TauMM;
			break;
		case Production::Channel::TauMT:
			gr = &ProductionSpace::TauMT;
			break;
		case Production::Channel::TauPI:				//1
			//gr = &ProductionSpace::TauPI;
			//break;
			return 1.;
		case Production::Channel::Tau2PI:				//1
			//gr = &ProductionSpace::Tau2PI;
			//break;
			return 1.;
		case Production::Channel::PionE:				//1
			//gr = &ProductionSpace::PionE;
			//break;
			return 1.;
		case Production::Channel::PionM:				//1
			//gr = &ProductionSpace::PionM;
			//break;
			return 1.;
		case Production::Channel::KaonE:				//1
			//gr = &ProductionSpace::KaonE;
			//break;
			return 1.;
		case Production::Channel::KaonM:				//1
			//gr = &ProductionSpace::KaonM;
			//break;
			return 1.;
		case Production::Channel::CharmE:				//1
			//gr = &ProductionSpace::CharmE;
			//break;
			return 1.;
		case Production::Channel::CharmM:				//1
			//gr = &ProductionSpace::CharmM;
			//break;
			return 1.;
		case Production::Channel::CharmT:				//1
			//gr = &ProductionSpace::CharmT;
			//break;
			return 1.;
		case Production::Channel::Kaon0E:
			gr = &ProductionSpace::Kaon0E;
			break;
		case Production::Channel::Kaon0M:
			gr = &ProductionSpace::Kaon0M;
			break;
		case Production::Channel::KaonCE:
			gr = &ProductionSpace::KaonCE;
			break;
		case Production::Channel::KaonCM:
			gr = &ProductionSpace::KaonCM;
			break;
		default:
			throw std::invalid_argument("Production channel "
					+ Production::toString(chan) + " is unknown");
	}
	
	return (this->*gr)(mix);
}


//// PRODUCTION modes

//pure leptonic decays
//
double ProductionSpace::MuonE(const Mixing &mix)
{
	return AntileptonNeutrino(Production::Channel::MuonE, Const::MMuon, Const::MElectron, Const::MNeutrino)
		/ _table[Production::Channel::MuonE];
}

double ProductionSpace::MuonM(const Mixing &mix)
{
	return LeptonNeutrino(Production::Channel::MuonM, Const::MMuon, Const::MElectron, Const::MNeutrino)
		/ _table[Production::Channel::MuonM];
}

double ProductionSpace::TauEE(const Mixing &mix)
{
	return AntileptonNeutrino(Production::Channel::TauEE, Const::MTau, Const::MElectron, Const::MNeutrino)
		/ _table[Production::Channel::TauEE];
}

double ProductionSpace::TauET(const Mixing &mix)
{
	return LeptonNeutrino(Production::Channel::TauET, Const::MTau, Const::MElectron, Const::MNeutrino)
		/ _table[Production::Channel::TauET];
}

double ProductionSpace::TauMM(const Mixing &mix)
{
	return AntileptonNeutrino(Production::Channel::TauMM, Const::MTau, Const::MMuon, Const::MNeutrino)
		/ _table[Production::Channel::TauMM];
}

double ProductionSpace::TauMT(const Mixing &mix)
{
	return LeptonNeutrino(Production::Channel::TauMT, Const::MTau, Const::MMuon, Const::MNeutrino)
		/ _table[Production::Channel::TauMT];
}

//	p1,u = N, p2,t = L, p3,s = n
double ProductionSpace::LeptonNeutrino(Production::Channel chan, double m_lepton0, double m_lepton, double m_neut)
{
	_m_parent = m_lepton0;
	double x = std::pow(m_lepton / _m_parent, 2);
	double y = std::pow(m_neut / _m_parent, 2);
	double z = std::pow(_N.M() / _m_parent, 2);

	if (!_table.count(chan)) {
		// compute max
		_vars = {x, y, z};
		
		auto func = std::bind(&ProductionSpace::LeptonNeutrino_max, this, std::placeholders::_1);
		double max_ = Optimization::GoldenRatio<double>(func); //-1 to invert function
		_table[chan] = dGammad5_3B() * (- func(max_));
	}

	std::array<double, 6> kine = Kinematic3B();// = s
	return dGammad5_3B() * M2_LeptonNeutrino(kine[2], x, y, z);
}

/*
double ProductionSpace::LeptonNeutrino_max(double x, double y, double z)
{
	_vars = {x, y, z};
	
	auto func = std::bind(&ProductionSpace::LeptonNeutrino_u_max, this, std::placeholders::_1);
	double res = Optimization::GoldenRatio<double>(func); //-1 to invert function
	return -func(res);
}
*/

double ProductionSpace::LeptonNeutrino_max(double u)	//vars is s 
{
	double x = _vars[0];
	double y = _vars[1];
	double z = _vars[2];

	// Limit changes u
	double fc = Limit(u, y, z, x);

	return - fc * M2_LeptonNeutrino(u, x, y, z);
}

double ProductionSpace::AntileptonNeutrino(Production::Channel chan, double m_lepton0, double m_lepton, double m_neut)
{
	_m_parent = m_lepton0;
	double x = std::pow(m_lepton / _m_parent, 2);
	double y = std::pow(m_neut / _m_parent, 2);
	double z = std::pow(_N.M() / _m_parent, 2);

	if (!_table.count(chan)) {
		// compute max
		_vars = {x, y, z};

		auto func = std::bind(&ProductionSpace::AntileptonNeutrino_max, this, std::placeholders::_1);
		double max_ = Optimization::GoldenRatio<double>(func); //-1 to invert function
		_table[chan] = dGammad5_3B() * (-func(max_));
	}

	std::array<double, 6> kine = Kinematic3B();// = s
	return dGammad5_3B() * M2_AntileptonNeutrino(kine[0], x, y, z);
}

/*
double ProductionSpace::AntileptonNeutrino_max(double x, double y, double z)	//vars is s 
{
	_vars = {x, y, z};

	auto func = std::bind(&ProductionSpace::AntileptonNeutrino_s_max, this, std::placeholders::_1);
	double res = Optimization::GoldenRatio<double>(func); //-1 to invert function
	return -func(res);
}
*/

double ProductionSpace::AntileptonNeutrino_max(double s)	//vars is s 
{
	double x = _vars[0];
	double y = _vars[1];
	double z = _vars[2];

	// limit changes s
	double fc = Limit(s, x, y, z);

	return - fc * M2_AntileptonNeutrino(s, x, y, z);
}


//lepton decay into meson stuff : production of 2 body decays is easy
/*
double ProductionSpace::TauPI(const Mixing &mix)
{
	return LeptonMeson(Production::Channel::TauPI, _masspdg[0].first, _masspdg[1].first) / _table[Production::Channel::TauPI];
}

double ProductionSpace::Tau2PI_ratio()
{
	return LeptonMeson(Production::Channel::Tau2PI, _masspdg[0].first, _masspgs[1].first) / _table[Production::Channel::Tau2PI];
}


double ProductionSpace::LeptonMeson(Production::Channel chan, double m_lepton, double m_meson)
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
double ProductionSpace::PionE(const Mixing &mix)
{
	return MesonTwo(Production::Channel::PionE, _masspdg[0].first, _masspdg[1].first) / _table[Production::Channel::PionE];
}

double ProductionSpace::PionM(const Mixing &mix)
{
	return MesonTwo(Production::Channel::PionM, _masspdg[0].first, _masspdg[1].first) / _table[Production::Channel::PionM];
}

double ProductionSpace::KaonE(const Mixing &mix)
{
	return MesonTwo(Production::Channel::KaonE, _masspdg[0].first, _masspdg[1].first) / _table[Production::Channel::KaonE];
}

double ProductionSpace::KaonM(const Mixing &mix)
{
	return MesonTwo(Production::Channel::KaonM, _masspdg[0].first, _masspdg[1].first) / _table[Production::Channel::KaonM];
}

double ProductionSpace::CharmE(const Mixing &mix)
{
	return MesonTwo(Production::Channel::CharmE, _masspdg[0].first, _masspdg[1].first) / _table[Production::Channel::CharmE];
}

double ProductionSpace::CharmM(const Mixing &mix)
{
	return MesonTwo(Production::Channel::CharmM, _masspdg[0].first, _masspdg[1].first) / _table[Production::Channel::CharmM];
}

double ProductionSpace::CharmT(const Mixing &mix)
{
	return MesonTwo(Production::Channel::CharmT, _masspdg[0].first, _masspdg[1].first) / _table[Production::Channel::CharmT];
}


// is the ratio just 1?
double ProductionSpace::MesonTwo(Production::Channel chan, double m_meson, double m_lepton)
{
	_m_parent = m_meson;
	double dMN2 = std::pow(_N.M() / _m_parent, 2);
	double dML2 = std::pow(m_lepton / _m_parent, 2);

	if (!_table.count(chan))
		_table[chan] = dGammad2_2B(dMN2, dML2) * M2_MesonTwo(dMN2, dML2);

	return dGammad2_2B(dMN2, dML2) * M2_MesonTwo(dMN2, dML2);
}
*/



//three body decays of meson
//
double ProductionSpace::Kaon0E(const Mixing &mix)
{
	return MesonThree(Production::Channel::Kaon0E, Const::MKaon0, Const::MPion, Const::MElectron,
				Const::KCL_, Const::KCL0) / _table[Production::Channel::Kaon0E];
}

double ProductionSpace::Kaon0M(const Mixing &mix)
{
	return MesonThree(Production::Channel::Kaon0M, Const::MKaon0, Const::MPion, Const::MMuon,
				Const::KCL_, Const::KCL0) / _table[Production::Channel::Kaon0M];
}

double ProductionSpace::KaonCE(const Mixing &mix)
{
	return MesonThree(Production::Channel::KaonCE, Const::MKaon, Const::MPion0, Const::MElectron,
				Const::KCL_, Const::KCL0) / _table[Production::Channel::KaonCE];
}

double ProductionSpace::KaonCM(const Mixing &mix)
{
	return MesonThree(Production::Channel::KaonCM, Const::MKaon, Const::MPion0, Const::MMuon,
				Const::KCL_, Const::KCL0) / _table[Production::Channel::KaonCM];
}


double ProductionSpace::MesonThree(Production::Channel chan, double m_meson0, double m_meson, double m_lepton,
				double L_, double L0)	//decay constant not important
{
	_m_parent = m_meson0;
	double x = std::pow(m_meson / _m_parent, 2);
	double y = std::pow(m_lepton / _m_parent, 2);
	double z = std::pow(_N.M() / _m_parent, 2);	

	if (!_table.count(chan)) {
		// find max
		_vars = {x, y, z, L_, L0};

		auto func = std::bind(&ProductionSpace::MesonThree_max, this, std::placeholders::_1);
		auto max_ = Optimization::NelderMead<double>(func, 2); //-1 to invert function
		_table[chan] = dGammad5_3B() * -func(&max_[0]);
	}

	std::array<double, 6> kine = Kinematic3B();	 // = s    = t
	return dGammad5_3B() * M2_MesonThree(kine[0], kine[1], x, y, z, L_, L0) ;
}

/*
double ProductionSpace::MesonThree_max(double x, double y, double z, double L_, double L0)	//no angle dependence
{
	_vars = {x, y, z, L_, L0};

	auto func = std::bind(&ProductionSpace::F_MesonThree_max, this, std::placeholders::_1);
	auto res = Optimization::NelderMead<double>(func, 2); //-1 to invert function
	return -func(&res[0]);
}
*/

// returns a reversed sign for minimization algorithm
double ProductionSpace::MesonThree_max(double p[])
{
	double x  = _vars[0];
	double y  = _vars[1];
	double z  = _vars[2];
	double gL = _vars[3];
	double gR = _vars[4];

	double s = p[0];
	double t = p[1];

	//double fc = Limit(s, t, x, y, x); ???

	return - M2_MesonThree(s, t, x, y, z, gL, gR);
}

void ProductionSpace::Reset() {
	_table.clear();
}
