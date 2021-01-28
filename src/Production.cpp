#include "physics/Production.h"

Production Production::_sm = Production(Neutrino());

Production::Production(Neutrino N) : Amplitude(std::move(N))
{
	// check if default argument is being used
	kSM = (_N == Neutrino());
}

bool Production::IsAllowed(Channel::Name chan)
{
	return (MassThreshold(chan) >= _N.M());
}

double Production::Scale(Channel::Name chan, const Mixing &mix)
{
	if (kSM)	// save computation time
		return 1.;
	return Gamma(chan, mix) / _sm.Gamma(chan, mix);
}

double Production::Gamma(Channel::Name chan, const Mixing &mix)
{
	if (!IsAllowed(chan))
		return 0.;

	using GammaF = double (Production::*)(const Mixing &mix);
	GammaF gf;

	switch(chan)
	{
		case Channel::ALL:
			gf = &Production::Total;
			break;
		case Channel::MuonE:
			gf = &Production::MuonE;
			break;
		case Channel::MuonM:
			gf = &Production::MuonM;
			break;
		case Channel::TauEE:
			gf = &Production::TauEE;
			break;
		case Channel::TauET:
			gf = &Production::TauET;
			break;
		case Channel::TauMM:
			gf = &Production::TauMM;
			break;
		case Channel::TauMT:
			gf = &Production::TauMT;
			break;
		case Channel::TauPI:
			gf = &Production::TauPI;
			break;
		case Channel::Tau2PI:
			gf = &Production::Tau2PI;
			break;
		case Channel::PionE:
			gf = &Production::PionE;
			break;
		case Channel::PionM:
			gf = &Production::PionM;
			break;
		case Channel::KaonE:
			gf = &Production::KaonE;
			break;
		case Channel::KaonM:
			gf = &Production::KaonM;
			break;
		case Channel::CharmE:
			gf = &Production::CharmE;
			break;
		case Channel::CharmM:
			gf = &Production::CharmM;
			break;
		case Channel::CharmT:
			gf = &Production::CharmT;
			break;
		case Channel::Kaon0E:
			gf = &Production::Kaon0E;
			break;
		case Channel::Kaon0M:
			gf = &Production::Kaon0M;
			break;
		case Channel::KaonCE:
			gf = &Production::KaonCE;
			break;
		case Channel::KaonCM:
			gf = &Production::KaonCM;
			break;
		default:
			throw std::invalid_argument("Channel " + Channel::toString(chan) + " unknown");
	}

	//return (Helicity() ? 1.0 : 2.0) * Result;
	return (this->*gf)(mix);
}

double Production::Total(const Mixing &mix)
{
	auto prods = Channel::Productions();
	return std::accumulate(prods.begin(), prods.end(), 0.,
			[=](double sum, const Channel::Name &chan) {
				return sum + Gamma(chan, mix); });
	// check capture here
}


double Production::MuonE(const Mixing &mix)
{
	if (!_table.count(Channel::MuonE))
		_table[Channel::MuonE] = AntileptonNeutrinoDecay(Const::MMuon, Const::MElectron,
							Const::MNeutrino);

	return _table[Channel::MuonE] * mix.Ue(2);
}

double Production::MuonM(const Mixing &mix)
{
	if (!_table.count(Channel::MuonM))
		_table[Channel::MuonM] = LeptonNeutrinoDecay(Const::MMuon, Const::MElectron,
						    Const::MNeutrino);

	return _table[Channel::MuonM] * mix.Um(2);
}

double Production::TauEE(const Mixing &mix)
{

	if (!_table.count(Channel::TauEE))
		_table[Channel::TauEE] = AntileptonNeutrinoDecay(Const::MTau, Const::MElectron,
							Const::MNeutrino);

	return _table[Channel::TauEE] * mix.Ue(2);
}

double Production::TauET(const Mixing &mix)
{

	if (!_table.count(Channel::TauET))
		_table[Channel::TauET] = LeptonNeutrinoDecay(Const::MTau, Const::MElectron,
						    Const::MNeutrino);

	return _table[Channel::TauET] * mix.Ut(2);
}

double Production::TauMM(const Mixing &mix)
{

	if (!_table.count(Channel::TauMM))
		_table[Channel::TauMM] = AntileptonNeutrinoDecay(Const::MTau, Const::MMuon,
							Const::MNeutrino);

	return _table[Channel::TauMM] * mix.Um(2);
}

double Production::TauMT(const Mixing &mix)
{

	if (!_table.count(Channel::TauMT))
		_table[Channel::TauMT] = LeptonNeutrinoDecay(Const::MTau, Const::MMuon,
						    Const::MNeutrino);

	return _table[Channel::TauMT] * mix.Ut(2);
}

double Production::TauPI(const Mixing &mix)
{

	if (!_table.count(Channel::TauPI))
		_table[Channel::TauPI] = std::pow(Const::U_ud * Const::DPion, 2)
				* LeptonTwoDecay(Const::MTau, Const::MPion);

	return _table[Channel::TauPI] * mix.Ut(2);
}

double Production::Tau2PI(const Mixing &mix)
{
	if (!_table.count(Channel::Tau2PI))
		_table[Channel::Tau2PI] = std::pow(Const::U_ud * Const::DPion, 2)
				* LeptonThreeDecay(Const::MTau, Const::MPion, Const::MPion0);

	return _table[Channel::Tau2PI] * mix.Ut(2);
}

double Production::PionE(const Mixing &mix)
{
	if (!_table.count(Channel::PionE))
		_table[Channel::PionE] = std::pow(Const::U_ud * Const::DPion, 2)
				* MesonTwoDecay(Const::MPion, Const::MElectron);

	return _table[Channel::PionE] * mix.Ue(2);
}

double Production::PionM(const Mixing &mix)
{
	if (!_table.count(Channel::PionM))
		_table[Channel::PionM] = std::pow(Const::U_ud * Const::DPion, 2)
				* MesonTwoDecay(Const::MPion, Const::MMuon);

	return _table[Channel::PionM] * mix.Um(2);
}

double Production::KaonE(const Mixing &mix)
{
	if (!_table.count(Channel::KaonE))
		_table[Channel::KaonE] = std::pow(Const::U_us * Const::DKaon, 2)
				* MesonTwoDecay(Const::MKaon, Const::MElectron);

	return _table[Channel::KaonE] * mix.Ue(2);
}

double Production::KaonM(const Mixing &mix)
{
	if (!_table.count(Channel::KaonM))
		_table[Channel::KaonM] = std::pow(Const::U_us * Const::DKaon, 2)
				* MesonTwoDecay(Const::MKaon, Const::MMuon);

	return _table[Channel::KaonM] * mix.Um(2);
}

double Production::CharmE(const Mixing &mix)
{
	if (!_table.count(Channel::CharmE))
		_table[Channel::CharmE] = std::pow(Const::U_cs * Const::DCharm, 2)
				* MesonTwoDecay(Const::MDs, Const::MElectron);

	return _table[Channel::CharmE] * mix.Ue(2);
}

double Production::CharmM(const Mixing &mix)
{
	if (!_table.count(Channel::CharmM))
		_table[Channel::CharmM] = std::pow(Const::U_cs * Const::DCharm, 2)
				* MesonTwoDecay(Const::MDs, Const::MMuon);

	return _table[Channel::CharmM] * mix.Um(2);
}

double Production::CharmT(const Mixing &mix)
{
	if (!_table.count(Channel::CharmT))
		_table[Channel::CharmT] = std::pow(Const::U_cs * Const::DCharm, 2)
				* MesonTwoDecay(Const::MDs, Const::MTau);

	return _table[Channel::CharmT] * mix.Ut(2);
}

double Production::Kaon0E(const Mixing &mix)
{
	if (!_table.count(Channel::Kaon0E))// || IsChanged())
		_table[Channel::Kaon0E] = std::pow(Const::U_us * Const::KaPi, 2)
			* MesonThreeDecay(Const::MKaon0, Const::MPion, Const::MElectron,
					  Const::K0L_, Const::K0L0);

	return _table[Channel::Kaon0E] * mix.Ue(2);
}

double Production::Kaon0M(const Mixing &mix)
{
	if (!_table.count(Channel::Kaon0M))// || IsChanged())
		_table[Channel::Kaon0M] = std::pow(Const::U_us * Const::KaPi, 2)
			* MesonThreeDecay(Const::MKaon0, Const::MPion, Const::MMuon,
					  Const::K0L_, Const::K0L0);

	return _table[Channel::Kaon0M] * mix.Um(2);
}

double Production::KaonCE(const Mixing &mix)
{
	if (!_table.count(Channel::KaonCE))
		_table[Channel::KaonCE] = std::pow(Const::U_us * Const::KaPi, 2) / 2.
			* MesonThreeDecay(Const::MKaon, Const::MPion0, Const::MElectron,
					  Const::KCL_, Const::KCL0);

	return _table[Channel::KaonCE] * mix.Ue(2);
}

double Production::KaonCM(const Mixing &mix)
{
	if (!_table.count(Channel::KaonCM))
		_table[Channel::KaonCM] = std::pow(Const::U_us * Const::KaPi, 2) / 2.
			* MesonThreeDecay(Const::MKaon, Const::MPion0, Const::MMuon,
					  Const::KCL_, Const::KCL0);

	return _table[Channel::KaonCM] * mix.Um(2);
}

/////////////////
//Generic decay//
/////////////////
//
//production from lepton -> neutrino is from same leptonic line
double Production::LeptonNeutrinoDecay(double m_lepton0, double m_lepton, double m_neut)
{
	_m_parent = m_lepton0;
	return I_LeptonNeutrino(std::pow(m_neut / _m_parent, 2),
				std::pow(m_lepton / _m_parent, 2),
				std::pow(_N.M() / _m_parent, 2));
}

//						c	  b	    a
double Production::I_LeptonNeutrino(double x, double y, double z)//, double theta)	//integrate first in y and then in x
{
	_vars = {x, y, z};
	//F_var.push_back(cos0);	//3	//theta

	auto func = std::bind(&Production::F_LeptonNeutrino_u, this, std::placeholders::_1);
	return Integration::Boole<1, double>(func); 
}

double Production::F_LeptonNeutrino_u(double u)	//fixing one variable
{
	//shouldnt use alias for clarity
	double x = _vars[0];	//light neutrino	u
	double y = _vars[1];	//lepton
	double z = _vars[2];	//heavy neutrino
	//const double &cos0 = F_var.at(3);

	double fc = Limit(u, y, z, x);

	//double cos0 = theta < 0 ? 0.0 : cos(theta);	//theta < 0 means integration over theta
	//double fc = theta < 0 ? 2.0 : 1.0;		//so an overall factor of 2

	return dGammad2_3B() * fc * M2_LeptonNeutrino(u, x, y, z);
}


//production from antilepton -> neutrino is from opposite leptonic line
double Production::AntileptonNeutrinoDecay(double m_lepton0, double m_lepton, double m_neut)
{
	_m_parent = m_lepton0;
	return I_AntileptonNeutrino(std::pow(m_neut / _m_parent, 2),
				    std::pow(m_lepton / _m_parent, 2),
				    std::pow(_N.M() / _m_parent, 2));
}
						//  c	      b		a
double Production::I_AntileptonNeutrino(double x, double y, double z)//, double theta)	//integrate first in y and then in x
{
	_vars = {x, y, z};
	//F_var.push_back(cos0);	//3	//theta

	auto func = std::bind(&Production::F_AntileptonNeutrino_s, this, std::placeholders::_1);
	return Integration::Boole<1, double>(func); 
}

double Production::F_AntileptonNeutrino_s(double s)//the term is written for a neutrino production
{						   //therefore with helicities inverted
	//aliases for clarity
	double x = _vars[0];
	double y = _vars[1];
	double z = _vars[2];
	//const double &cos0 = F_var.at(5);

	double fc = Limit(s, x, y, z);

	//double cos0 = theta < 0 ? 0.0 : cos(theta);	//theta < 0 means integration over theta
	//double fc = theta < 0 ? 2.0 : 1.0;		//so an overall factor of 2

	return dGammad2_3B() * fc * M2_AntileptonNeutrino(s, x, y, z);
}

double Production::LeptonTwoDecay(double m_lepton, double m_meson)
{
	_m_parent = m_lepton;
	return I_LeptonTwo(std::pow(_N.M() / _m_parent, 2),
			   std::pow(m_meson / _m_parent, 2));
}

double Production::I_LeptonTwo(double x, double y)
{
	return dGammad0_2B(x, y) * M2_LeptonTwo(x, y);
}

double Production::LeptonThreeDecay(double m_lepton, double m_meson0, double m_meson)
{
	_m_parent = m_lepton;
	return I_LeptonThree(std::pow(_N.M() / _m_parent, 2),
			     std::pow(m_meson / _m_parent, 2),
			     std::pow(m_meson0 / _m_parent, 2));
}

double Production::I_LeptonThree(double x, double y, double z)
{
	_vars = {x, y, z};

	auto func = std::bind(&Production::F_LeptonThree_s, this, std::placeholders::_1);
	return Integration::Boole<1, double>(func); 
}

double Production::F_LeptonThree_s(double s)
{
	double x = _vars[0];
	double y = _vars[1];
	double z = _vars[2];
	double fc = Limit(s, x, y, z);

	return dGammad2_3B() * fc * M2_LeptonThree(x, y, z);
}

double Production::MesonTwoDecay(double m_meson, double m_lepton)
{
	_m_parent = m_meson;
	return I_MesonTwo(std::pow(_N.M() / _m_parent, 2),
			  std::pow(m_lepton / _m_parent, 2));
}

double Production::I_MesonTwo(double x, double y)	//symetric in x and y
{
	return dGammad0_2B(x, y) * M2_MesonTwo(x, y);
}

double Production::MesonThreeDecay(double m_meson0, double m_meson, double m_lepton, double L_, double L0)	//decay constant not important
{
	_m_parent = m_meson0;
	return I_MesonThree(std::pow(m_meson / _m_parent, 2),
			    std::pow(m_lepton / _m_parent, 2),
			    std::pow(_N.M() / _m_parent, 2),
			    L_, L0);
}

//no angle dependence
double Production::I_MesonThree(double x, double y, double z, double L_, double L0)
{
	_vars = {x, y, z, L_, L0};

	auto func = std::bind(&Production::F_MesonThree_s_t, this,
			      std::placeholders::_1, std::placeholders::_2);
	return Integration::Boole<2, double>(func); // reduced precision, 16x faster
}

double Production::F_MesonThree_s_t(double s, double t)
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

	return dGammad2_3B() * fc * M2_MesonThree(s, t, x, y, z, L_, L0);
}
