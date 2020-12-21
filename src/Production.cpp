#include "Production.h"

Production::Production(const HNL &N) : Amplitude(N)
{
}

double Production::MassThreshold(Channel::Name chan)
{
	if (_channel != chan) {
		LoadMass(chan);
		_channel = chan;
	}

	return std::accumulate(_masspdg.begin()+1, _masspdg.end(), _masspdg.front().first,
			[](double sum, const std::pair<double, int> mp) {
				return sum - mp.first; });
}

bool Production::IsAllowed(Channel::Name chan)
{
	return (MassThreshold(chan) >= _N.M());
}

double Production::Scale(Channel::Name chan)
{
	return Gamma(chan, 1., 1., 1.);
}

double Production::Gamma(Channel::Name chan, std::arrya<double, 3> &mix) {
	return Gamma(chan, mix[0], mix[1], mix[2]);
}

double Production::Gamma(Channel::Name chan, double ue, double um, double ut)
{
	if (!IsAllowed(chan))
		return 0.;

	using GammaF = double (Production::*)(double, double, double);
	GammaF gf;

	switch(chan)
	{
		case Channel::_ALL:
			gf = &Production::Total;
			break;
		case Channel::_MuonE:
			gf = &Production::MuonE;
			break;
		case Channel::_MuonM:
			gf = &Production::MuonM;
			break;
		case Channel::_TauEE:
			gf = &Production::TauEE;
			break;
		case Channel::_TauET:
			gf = &Production::TauET;
			break;
		case Channel::_TauMM:
			gf = &Production::TauMM;
			break;
		case Channel::_TauMT:
			gf = &Production::TauMT;
			break;
		case Channel::_TauPI:
			gf = &Production::TauPI;
			break;
		case Channel::_Tau2PI:
			gf = &Production::Tau2PI;
			break;
		case Channel::_PionE:
			gf = &Production::PionE;
			break;
		case Channel::_PionM:
			gf = &Production::PionM;
			break;
		case Channel::_KaonE:
			gf = &Production::KaonE;
			break;
		case Channel::_KaonM:
			gf = &Production::KaonM;
			break;
		case Channel::_CharmE:
			gf = &Production::ChamE;
			break;
		case Channel::_CharmM:
			gf = &Production::ChamM;
			break;
		case Channel::_CharmT:
			gf = &Production::CharmT;
			break;
		case Channel::_Kaon0E:
			gf = &Production::Kaon0E;
			break;
		case Channel::_Kaon0M:
			gf = &Production::Kaon0M;
			break;
		case Channel::_KaonCE:
			gf = &Production::KaonCE;
			break;
		case Channel::_KaonCM:
			gf = &Production::KaonCM;
			break;
		default:
			throw std::invalid_argument("Channel " + toString(chan) + " unknown");
	}

	//return (Helicity() ? 1.0 : 2.0) * Result;
	return (this->*gf)(ue, um, ut);
}

double Production::Total(double ue, double um, double ut)
{
	return std::accumulate(std::begin(Channel::Productions), std::end(Channel::Productions),
			0., [=](double sum, const Channel::Name &chan) { return sum + Gamma(chan, ue, um, ut); });
	// check capture here
}


double Production::MuonE(double ue, double um, double ut)
{
	if (!fprod.count(Channel::_MuonE))
		fprod[Channel::_MuonE] = AntiLeptonNeutrinoDecay(Const::MMuon, Const::MElectron,
							Const::MNeutrino)

	return fprod[Channel::_MuonE] * ue*ue;
}

double Production::MuonM(double ue, double um, double ut)
{
	if (!fprod.count(Channel::_MuonM))
		fprod[Channel::_MuonM] = LeptonNeutrinoDecay(Const::MMuon, Const::MElectron,
						    Const::MNeutrino);

	return fprod[Channel::_MuonM] * um*um;
}

double Production::TauEE(double ue, double um, double ut)
{

	if (!fprod.count(Channel::_TauEE))
		fprod[Channel::_TauEE] = AntiLeptonNeutrinoDecay(Const::MTau, Const::MElectron,
							Const::MNeutrino);

	return fprod[Channel::_TauEE] * ue*ue;
}

double Production::TauET(double ue, double um, double ut)
{

	if (!fprod.count(Channel::_TauET))
		fprod[Channel::_TauET] = LeptonNeutrinoDecay(Const::MTau, Const::MElectron,
						    Const::MNeutrino);

	return fprod[Channel::_TauET] * ut*ut;
}

double Production::TauMM(double ue, double um, double ut)
{

	if (!fprod.count(Channel::_TauMM))
		fprod[Channel::_TauMM] = AntiLeptonNeutrinoDecay(Const::MTau, Const::MMuon,
							Const::MNeutrino);

	return fprod[Channel::_TauMM] * um*um;
}

double Production::TauMT(double ue, double um, double ut)
{

	if (!fprod.count(Channel::_TauMT))
		fprod[Channel::_TauMT] = LeptonNeutrinoDecay(Const::MTau, Const::MMuon,
						    Const::MNeutrino);

	return fprod[Channel::_TauMT] * ut*ut;
}

double Production::TauPI(double ue, double um, double ut)
{

	if (!fprod.count(Channel::_TauPI))
		fprod[Channel::_TauPI] = std::pow(Const::U_ud * Const::DPion, 2)
				* LeptonTwoDecay(Const::MTau, Const::MPion);

	return fprod[Channel::_TauPI] * ut*ut;
}

double Production::Tau2PI(double ue, double um, double ut)
{
	if (!fprod.count(Channel::_Tau2PI))
		fprod[Channel::_Tau2PI] = std::pow(Const::U_ud * Const::DPion, 2)
				* LeptonThreeDecay(Const::MTau, Const::MPion, Const::MPion0);

	return fprod[Channel::_Tau2PI] * ut*ut;
}

double Production::PionE(double ue, double um, double ut)
{
	if (!fprod.count(Channel::_PionE))
		fprod[Channel::_PionE] = std::pow(Const::U_ud * Const::DPion, 2)
				* MesonTwoDecay(Const::MPion, Const::MElectron);

	return fprod[Channel::_PionE] * ue*ue;
}

double Production::PionM(double ue, double um, double ut)
{

	if (!fprod.count(Channel::_PionM))
		fprod[Channel::_PionM] = std::pow(Const::U_ud * Const::DPion, 2)
				* MesonTwoDecay(Const::MPion, Const::MMuon);

	return fprod[Channel::_PionM] * um*um;
}

double Production::KaonE(double ue, double um, double ut)
{

	if (!fprod.count(Channel::_KaonE))
		fprod[Channel::_KaonE] = std::pow(Const::U_us * Const::DKaon, 2)
				* MesonTwoDecay(Const::MKaon, Const::MElectron);

	return fprod[Channel::_KaonE] * ue*ue;
}

double Production::KaonM(double ue, double um, double ut)
{

	if (!fprod.count(Channel::_KaonM))
		fprod[Channel::_KaonM] = std::pow(Const::U_us * Const::DKaon, 2)
				* MesonTwoDecay(Const::MKaon, Const::MMuon);

	return fprod[Channel::_KaonM] * um*um;
}

double Production::CharmE(double ue, double um, double ut)
{

	if (!fprod.count(Channel::_CharmE))
		fprod[Channel::_CharmE] = std::pow(Const::U_cs * Const::DCharm, 2)
				* MesonTwoDecay(Const::MDs, Const::MElectron);

	return fprod[Channel::_CharmE] * ue*ue;
}

double Production::CharmM(double ue, double um, double ut)
{

	if (!fprod.count(Channel::_CharmM))
		fprod[Channel::_CharmM] = std::pow(Const::U_cs * Const::DCharm, 2)
				* MesonTwoDecay(Const::MDs, Const::MMuon);

	return fprod[Channel::_CharmM] * um*um;
}

double Production::CharmT(double ue, double um, double ut)
{

	if (!fprod.count(Channel::_CharmT))
		fprod[Channel::_CharmT] = std::pow(Const::U_cs * Const::DCharm, 2)
				* MesonTwoDecay(Const::MDs, Const::MTau);

	return fprod[Channel::_CharmT] * ut*ut;
}

double Production::Kaon0E(double ue, double um, double ut)
{

	if (!fprod.count(Channel::_Kaon0E))// || IsChanged())
		fprod[Channel::_Kaon0E] = std::pow(Const::U_us * Const::KaPi, 2)
			* MesonThreeDecay(Const::MKaon0, Const::MPion, Const::MElectron,
					  Const::K0L_, Const::K0L0);

	return fprod[Channel::_Kaon0E] * ue*ue;
}

double Production::Kaon0M(double ue, double um, double ut)
{

	if (!fprod.count(Channel::_Kaon0M))// || IsChanged())
		fprod[Channel::_Kaon0M] = std::pow(Const::U_us * Const::KaPi, 2)
			* MesonThreeDecay(Const::MKaon0, Const::MPion, Const::MMuon,
					  Const::K0L_, Const::K0L0);

	return fprod[Channel::_Kaon0M] * um*um;
}

double Production::KaonCE(double ue, double um, double ut)
{
	if (!fprod.count(Channel::_KaonCE))
		fprod[Channel::_KaonCE] = std::pow(Const::U_us * Const::KaPi, 2) / 2.
			* MesonThreeDecay(Const::MKaon, Const::MPion0, Const::MElectron,
					  Const::KCL_, Const::KCL0);

	return fprod[Channel::_KaonCE] * ue*ue;
}

double Production::KaonCM(double ue, double um, double ut)
{
	if (!fprod.count(Channel::_KaonCM))
		fprod[Channel::_KaonCM] = std::pow(Const::U_us * Const::KaPi, 2) / 2.
			* MesonThreeDecay(Const::MKaon, Const::MPion0, Const::MMuon,
					  Const::KCL_, Const::KCL0);

	return fprod[Channel::_KaonCM] * um*um;
}

/////////////////
//Generic decay//
/////////////////
//
//production from lepton -> neutrino is from same leptonic line
double Production::LeptonNeutrinoDecay(double m_lepton0, double m_lepton, double m_neut)
{
	_m_parent = _m_lepton0;
	return I_LeptonNeutrino(std::pow(m_neut / _m_parent, 2),
				std::pow(m_lepton / _m_parent, 2),
				std::pow(_m_Nu / _m_parent, 2));
}

//						c	  b	    a
double Production::I_LeptonNeutrino(double x, double y, double z)//, double theta)	//integrate first in y and then in x
{
	_vars = {x, y, z};
	//F_var.push_back(cos0);	//3	//theta

	return BooleIntegration(&Production::I_LeptonNeutrino_u); 
}

double Production::I_LeptonNeutrino_u(double u)	//fixing one variable
{
	//shouldnt use alias for clarity
	double x = _vars[0];	//light neutrino	u
	double y = _vars[1];	//lepton
	double z = _vars[2];	//heavy neutrino
	//const double &cos0 = F_var.at(3);

	double u_ = u;
	double fc = Limit(u_, y, z, x);

	//double cos0 = theta < 0 ? 0.0 : cos(theta);	//theta < 0 means integration over theta
	//double fc = theta < 0 ? 2.0 : 1.0;		//so an overall factor of 2

	double M2 = fc * M2_LeptonNeutrino(u_, x, y, z);
	return dGammad2_3B(M2);
}


//production from antilepton -> neutrino is from opposite leptonic line
double Production::AntiLeptonNeutrinoDecay(double m_lepton0, double m_lepton, double m_neut)
{
	_m_parent = m_lepton0;
	return I_AntiLeptonNeutrino(std::pow(m_neut / _m_parent, 2),
				    std::pow(m_lepton / _m_parent, 2),
				    std::pow(_N.M() / _m_parent, 2));
}
						//  c	      b		a
double Production::I_AntiLeptonNeutrino(double x, double y, double z)//, double theta)	//integrate first in y and then in x
{
	_vars = {x, y, z};
	//F_var.push_back(cos0);	//3	//theta

	return BooleIntegration(&Production::I_AntiLeptonNeutrino_s);
}

double Production::I_AntiLeptonNeutrino_s(double s)	//the term is written for a neutrino production
{								//therefore with heliciies inverted
	//aliases for clarity
	double x = _vars[0];
	double y = _vars[1];
	double z = _vars[2];
	//const double &cos0 = F_var.at(5);

	//create S var
	double s_ = s;
	double fc = Limit(s_, x, y, z);

	//double cos0 = theta < 0 ? 0.0 : cos(theta);	//theta < 0 means integration over theta
	//double fc = theta < 0 ? 2.0 : 1.0;		//so an overall factor of 2

	double M2 = fc * M2_AntiLeptonNeutrino(s_, x, y, z);
	return dGammad2_3B(M2);
}

double Production::LeptonTwoDecay(double m_lepton, double m_meson)
{
	_m_parent = m_lepton;
	return I_LeptonTwo(std::pow(_N.M() / _m_parent, 2),
			   std::pow(m_meson / _m_parent, 2));
}

double Production::I_LeptonTwo(double x, double y)
{
	double M2 = M2_LeptonTwo(x, y);
	return dGammad0_2B(M2, x, y);
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

	return BooleIntegration(&Production::I_LeptonThree_s); 		//switch to Vega
}

double Production::I_LeptonThree_s(double s)
{
	double x = _vars[0];
	double y = _vars[1];
	double z = _vars[2];
	double s_ = s;
	double fc = Limit(s_, x, y, z);

	double M2 = fc * M2_LeptonThree(x, y, z);
	return dGammad2_3B(M2);
}

double Production::MesonTwoDecay(double m_meson, double m_lepton)
{
	_m_parent = m_meson;
	return I_MesonTwo(std::pow(_N.M() / _m_parent, 2),
			  std::pow(m_lepton / _m_parent, 2));
}

double Production::I_MesonTwo(double x, double y)	//symetric in x and y
{
	double M2 = M2_MesonTwo(x, y);
	return dGammad0_2B(M2, x, y);
}

double Production::MesonThreeDecay(double m_meson0, double m_meson, double m_lepton, double L_, double L0)	//decay constant not important
{
	_m_parent = m_meson0;
	return I_MesonThree(std::pow(m_meson / _m_parent, 2),
			    std::pow(m_lepton / _m_parent, 2),
			    std::pow(_N.M() / _m_parent, 2),
			    L_, L0);
}

double Production::I_MesonThree(double x, double y, double z, double L_, double L0)	//no angle dependence
{
	_vars = {x, y, z, L_, L0, 0};

	return BooleIntegration(&Production::I_MesonThree_s); 		//switch to Vega
}

// nested integration
double Production::I_MesonThree_s(double s)	//fixing one variable
{
	_vars[5] = s;
	return BooleIntegration(&Production::I_MesonThree_t);
}

double Production::I_MesonThree_t(double t)
{
	//aliases for clarity
	double x  = _vars[0];
	double y  = _vars[1];
	double z  = _vars[2];
	double L0 = _vars[3];
	double L_ = _vars[4];

	double s_ = _vars[5];
	double t_ = t;
	double fc = Limit(s_, t_, x, y, z);

	//double cos0 = theta < 0 ? 0.0 : cos(theta);	//theta < 0 means integration over theta
	//double fc = theta < 0 ? 2.0 : 1.0;              //so an overall factor of 2

	double M2 = fc * M2_MesonThree(s_, t_, x, y, z, L_, L0);
	return dGammad2_3B(M2);
}

double Reset() {
	fprod.clear();
}
