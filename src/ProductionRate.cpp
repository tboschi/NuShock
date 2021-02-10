#include "physics/ProductionRate.h"

ProductionRate ProductionRate::_sm = ProductionRate(Neutrino());

ProductionRate::ProductionRate(Neutrino N) : Amplitude(std::move(N))
{
	// check if default argument is being used
	kSM = (_N == Neutrino());
}

double ProductionRate::Scale(Production::Channel chan, const Mixing &mix)
{
	if (kSM)	// save computation time
		return 1.;
	return Gamma(chan, mix) / _sm.Gamma(chan, mix);
}

double ProductionRate::Gamma(Production::Channel chan, const Mixing &mix)
{
	if (!IsAllowed(chan))
		return 0.;

	using GammaF = double (ProductionRate::*)(const Mixing &mix);
	GammaF gf;

	switch(chan)
	{
		//case Production::Channel::ALL:
		//	gf = &ProductionRate::Total;
		//	break;
		case Production::Channel::MuonE:
			gf = &ProductionRate::MuonE;
			break;
		case Production::Channel::MuonM:
			gf = &ProductionRate::MuonM;
			break;
		case Production::Channel::TauEE:
			gf = &ProductionRate::TauEE;
			break;
		case Production::Channel::TauET:
			gf = &ProductionRate::TauET;
			break;
		case Production::Channel::TauMM:
			gf = &ProductionRate::TauMM;
			break;
		case Production::Channel::TauMT:
			gf = &ProductionRate::TauMT;
			break;
		case Production::Channel::TauPI:
			gf = &ProductionRate::TauPI;
			break;
		case Production::Channel::Tau2PI:
			gf = &ProductionRate::Tau2PI;
			break;
		case Production::Channel::PionE:
			gf = &ProductionRate::PionE;
			break;
		case Production::Channel::PionM:
			gf = &ProductionRate::PionM;
			break;
		case Production::Channel::KaonE:
			gf = &ProductionRate::KaonE;
			break;
		case Production::Channel::KaonM:
			gf = &ProductionRate::KaonM;
			break;
		case Production::Channel::CharmE:
			gf = &ProductionRate::CharmE;
			break;
		case Production::Channel::CharmM:
			gf = &ProductionRate::CharmM;
			break;
		case Production::Channel::CharmT:
			gf = &ProductionRate::CharmT;
			break;
		case Production::Channel::Kaon0E:
			gf = &ProductionRate::Kaon0E;
			break;
		case Production::Channel::Kaon0M:
			gf = &ProductionRate::Kaon0M;
			break;
		case Production::Channel::KaonCE:
			gf = &ProductionRate::KaonCE;
			break;
		case Production::Channel::KaonCM:
			gf = &ProductionRate::KaonCM;
			break;
		default:
			throw std::invalid_argument("Production channel "
					+ Production::toString(chan) + " is unknown");
	}

	//return (Helicity() ? 1.0 : 2.0) * Result;
	return (this->*gf)(mix);
}

double ProductionRate::Total(const Mixing &mix)
{
	auto prods = Production::Channels();
	return std::accumulate(prods.begin(), prods.end(), 0.,
			[=](double sum, Production::Channel chan) {
				return sum + Gamma(chan, mix); });
	// check capture here
}


double ProductionRate::MuonE(const Mixing &mix)
{
	if (!_table.count(Production::Channel::MuonE))
		_table[Production::Channel::MuonE] = AntileptonNeutrinoDecay(Const::MMuon, Const::MElectron,
							Const::MNeutrino);

	return _table[Production::Channel::MuonE] * mix.Ue(2);
}

double ProductionRate::MuonM(const Mixing &mix)
{
	if (!_table.count(Production::Channel::MuonM))
		_table[Production::Channel::MuonM] = LeptonNeutrinoDecay(Const::MMuon, Const::MElectron,
						    Const::MNeutrino);

	return _table[Production::Channel::MuonM] * mix.Um(2);
}

double ProductionRate::TauEE(const Mixing &mix)
{

	if (!_table.count(Production::Channel::TauEE))
		_table[Production::Channel::TauEE] = AntileptonNeutrinoDecay(Const::MTau, Const::MElectron,
							Const::MNeutrino);

	return _table[Production::Channel::TauEE] * mix.Ue(2);
}

double ProductionRate::TauET(const Mixing &mix)
{

	if (!_table.count(Production::Channel::TauET))
		_table[Production::Channel::TauET] = LeptonNeutrinoDecay(Const::MTau, Const::MElectron,
						    Const::MNeutrino);

	return _table[Production::Channel::TauET] * mix.Ut(2);
}

double ProductionRate::TauMM(const Mixing &mix)
{

	if (!_table.count(Production::Channel::TauMM))
		_table[Production::Channel::TauMM] = AntileptonNeutrinoDecay(Const::MTau, Const::MMuon,
							Const::MNeutrino);

	return _table[Production::Channel::TauMM] * mix.Um(2);
}

double ProductionRate::TauMT(const Mixing &mix)
{

	if (!_table.count(Production::Channel::TauMT))
		_table[Production::Channel::TauMT] = LeptonNeutrinoDecay(Const::MTau, Const::MMuon,
						    Const::MNeutrino);

	return _table[Production::Channel::TauMT] * mix.Ut(2);
}

double ProductionRate::TauPI(const Mixing &mix)
{

	if (!_table.count(Production::Channel::TauPI))
		_table[Production::Channel::TauPI] = std::pow(Const::U_ud * Const::DPion, 2)
				* LeptonTwoDecay(Const::MTau, Const::MPion);

	return _table[Production::Channel::TauPI] * mix.Ut(2);
}

double ProductionRate::Tau2PI(const Mixing &mix)
{
	if (!_table.count(Production::Channel::Tau2PI))
		_table[Production::Channel::Tau2PI] = std::pow(Const::U_ud * Const::DPion, 2)
				* LeptonThreeDecay(Const::MTau, Const::MPion, Const::MPion0);

	return _table[Production::Channel::Tau2PI] * mix.Ut(2);
}

double ProductionRate::PionE(const Mixing &mix)
{
	if (!_table.count(Production::Channel::PionE))
		_table[Production::Channel::PionE] = std::pow(Const::U_ud * Const::DPion, 2)
				* MesonTwoDecay(Const::MPion, Const::MElectron);

	return _table[Production::Channel::PionE] * mix.Ue(2);
}

double ProductionRate::PionM(const Mixing &mix)
{
	if (!_table.count(Production::Channel::PionM))
		_table[Production::Channel::PionM] = std::pow(Const::U_ud * Const::DPion, 2)
				* MesonTwoDecay(Const::MPion, Const::MMuon);

	return _table[Production::Channel::PionM] * mix.Um(2);
}

double ProductionRate::KaonE(const Mixing &mix)
{
	if (!_table.count(Production::Channel::KaonE))
		_table[Production::Channel::KaonE] = std::pow(Const::U_us * Const::DKaon, 2)
				* MesonTwoDecay(Const::MKaon, Const::MElectron);

	return _table[Production::Channel::KaonE] * mix.Ue(2);
}

double ProductionRate::KaonM(const Mixing &mix)
{
	if (!_table.count(Production::Channel::KaonM))
		_table[Production::Channel::KaonM] = std::pow(Const::U_us * Const::DKaon, 2)
				* MesonTwoDecay(Const::MKaon, Const::MMuon);

	return _table[Production::Channel::KaonM] * mix.Um(2);
}

double ProductionRate::CharmE(const Mixing &mix)
{
	if (!_table.count(Production::Channel::CharmE))
		_table[Production::Channel::CharmE] = std::pow(Const::U_cs * Const::DCharm, 2)
				* MesonTwoDecay(Const::MDs, Const::MElectron);

	return _table[Production::Channel::CharmE] * mix.Ue(2);
}

double ProductionRate::CharmM(const Mixing &mix)
{
	if (!_table.count(Production::Channel::CharmM))
		_table[Production::Channel::CharmM] = std::pow(Const::U_cs * Const::DCharm, 2)
				* MesonTwoDecay(Const::MDs, Const::MMuon);

	return _table[Production::Channel::CharmM] * mix.Um(2);
}

double ProductionRate::CharmT(const Mixing &mix)
{
	if (!_table.count(Production::Channel::CharmT))
		_table[Production::Channel::CharmT] = std::pow(Const::U_cs * Const::DCharm, 2)
				* MesonTwoDecay(Const::MDs, Const::MTau);

	return _table[Production::Channel::CharmT] * mix.Ut(2);
}

double ProductionRate::Kaon0E(const Mixing &mix)
{
	if (!_table.count(Production::Channel::Kaon0E))// || IsChanged())
		_table[Production::Channel::Kaon0E] = std::pow(Const::U_us * Const::KaPi, 2)
			* MesonThreeDecay(Const::MKaon0, Const::MPion, Const::MElectron,
					  Const::K0L_, Const::K0L0);

	return _table[Production::Channel::Kaon0E] * mix.Ue(2);
}

double ProductionRate::Kaon0M(const Mixing &mix)
{
	if (!_table.count(Production::Channel::Kaon0M))// || IsChanged())
		_table[Production::Channel::Kaon0M] = std::pow(Const::U_us * Const::KaPi, 2)
			* MesonThreeDecay(Const::MKaon0, Const::MPion, Const::MMuon,
					  Const::K0L_, Const::K0L0);

	return _table[Production::Channel::Kaon0M] * mix.Um(2);
}

double ProductionRate::KaonCE(const Mixing &mix)
{
	if (!_table.count(Production::Channel::KaonCE))
		_table[Production::Channel::KaonCE] = std::pow(Const::U_us * Const::KaPi, 2) / 2.
			* MesonThreeDecay(Const::MKaon, Const::MPion0, Const::MElectron,
					  Const::KCL_, Const::KCL0);

	return _table[Production::Channel::KaonCE] * mix.Ue(2);
}

double ProductionRate::KaonCM(const Mixing &mix)
{
	if (!_table.count(Production::Channel::KaonCM))
		_table[Production::Channel::KaonCM] = std::pow(Const::U_us * Const::KaPi, 2) / 2.
			* MesonThreeDecay(Const::MKaon, Const::MPion0, Const::MMuon,
					  Const::KCL_, Const::KCL0);

	return _table[Production::Channel::KaonCM] * mix.Um(2);
}

/////////////////
//Generic decay//
/////////////////
//
//production from lepton -> neutrino is from same leptonic line
double ProductionRate::LeptonNeutrinoDecay(double m_lepton0, double m_lepton, double m_neut)
{
	_m_parent = m_lepton0;
	double x = std::pow(m_neut / _m_parent, 2);
	double y = std::pow(m_lepton / _m_parent, 2);
	double z = std::pow(_N.M() / _m_parent, 2);

	_vars = {x, y, z};
	auto func = std::bind(&ProductionRate::LeptonNeutrino_u, this, std::placeholders::_1);
	return Integration::Boole<1, double>(func); 
}

double ProductionRate::LeptonNeutrino_u(double u)	//fixing one variable
{
	//shouldnt use alias for clarity
	double x = _vars[0];	//light neutrino	u
	double y = _vars[1];	//lepton
	double z = _vars[2];	//heavy neutrino
	//const double &cos0 = F_var.at(3);

	double fc = Limit(u, y, z, x);

	//double cos0 = theta < 0 ? 0.0 : cos(theta);	//theta < 0 means integration over theta
	//double fc = theta < 0 ? 2.0 : 1.0;		//so an overall factor of 2

	return fc * M2_LeptonNeutrino(u, x, y, z);
}


//production from antilepton -> neutrino is from opposite leptonic line
double ProductionRate::AntileptonNeutrinoDecay(double m_lepton0, double m_lepton, double m_neut)
{
	_m_parent = m_lepton0;

	double x = std::pow(m_neut / _m_parent, 2);
	double y = std::pow(m_lepton / _m_parent, 2);
	double z = std::pow(_N.M() / _m_parent, 2);

	_vars = {x, y, z};
	auto func = std::bind(&ProductionRate::AntileptonNeutrino_s, this, std::placeholders::_1);
	return dGammad2_3B() * Integration::Boole<1, double>(func); 
}

double ProductionRate::AntileptonNeutrino_s(double s)//the term is written for a neutrino production
{						   //therefore with helicities inverted
	//aliases for clarity
	double x = _vars[0];
	double y = _vars[1];
	double z = _vars[2];
	//const double &cos0 = F_var.at(5);

	double fc = Limit(s, x, y, z);

	//double cos0 = theta < 0 ? 0.0 : cos(theta);	//theta < 0 means integration over theta
	//double fc = theta < 0 ? 2.0 : 1.0;		//so an overall factor of 2

	return fc * M2_AntileptonNeutrino(s, x, y, z);
}

double ProductionRate::LeptonTwoDecay(double m_lepton, double m_meson)
{
	_m_parent = m_lepton;
	double x = std::pow(_N.M() / _m_parent, 2);
	double y = std::pow(m_meson / _m_parent, 2);
			   
	return dGammad0_2B(x, y) * M2_LeptonTwo(x, y);
}

double ProductionRate::LeptonThreeDecay(double m_lepton, double m_meson0, double m_meson)
{
	_m_parent = m_lepton;
	double x = std::pow(_N.M() / _m_parent, 2);
	double y = std::pow(m_meson / _m_parent, 2);
	double z = std::pow(m_meson0 / _m_parent, 2);

	_vars = {x, y, z};

	auto func = std::bind(&ProductionRate::LeptonThree_s, this, std::placeholders::_1);
	return dGammad2_3B() * Integration::Boole<1, double>(func); 
}

double ProductionRate::LeptonThree_s(double s)
{
	double x = _vars[0];
	double y = _vars[1];
	double z = _vars[2];
	double fc = Limit(s, x, y, z);

	return fc * M2_LeptonThree(x, y, z);
}

double ProductionRate::MesonTwoDecay(double m_meson, double m_lepton)
{
	_m_parent = m_meson;
	double x = std::pow(_N.M() / _m_parent, 2);
	double y = std::pow(m_lepton / _m_parent, 2);
	return dGammad0_2B(x, y) * M2_MesonTwo(x, y);
}

double ProductionRate::MesonThreeDecay(double m_meson0, double m_meson, double m_lepton,
					double L_, double L0)	//decay constant not important
{
	_m_parent = m_meson0;
	double x = std::pow(m_meson / _m_parent, 2);
	double y = std::pow(m_lepton / _m_parent, 2);
	double z = std::pow(_N.M() / _m_parent, 2);
			    
	_vars = {x, y, z, L_, L0};

	auto func = std::bind(&ProductionRate::MesonThree_s_t, this,
			      std::placeholders::_1, std::placeholders::_2);
	return dGammad2_3B() * Integration::Boole<2, double>(func, 250); 
}

double ProductionRate::MesonThree_s_t(double s, double t)
{
	//aliases for clarity
	double x  = _vars[0];
	double y  = _vars[1];
	double z  = _vars[2];
	double L_ = _vars[3];
	double L0 = _vars[4];

	double fc = Limit(s, t, x, y, z);

	//double cos0 = theta < 0 ? 0.0 : cos(theta);	//theta < 0 means integration over theta
	//double fc = theta < 0 ? 2.0 : 1.0;              //so an overall factor of 2

	return fc * M2_MesonThree(s, t, x, y, z, L_, L0);
}

void ProductionRate::Reset() {
	_table.clear();
}
