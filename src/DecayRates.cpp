#include "physics/DecayRates.h"

DecayRates::DecayRates(Neutrino N) : Amplitude(std::move(N))
{
}

bool DecayRates::IsAllowed(Channel::Name chan)
{
	return (_N.M() >= MassThreshold(chan));
}

//return the decay width (Gamma)
//
double DecayRates::Gamma(Channel::Name chan, const Mixing &mix)
{
	if (!IsAllowed(chan))
		return 0.;

	using GammaF = double (DecayRates::*)(const Mixing &mix);
	GammaF gf;

	switch(chan)
	{
		case Channel::ALL:
			gf = &DecayRates::Total;
			break;
		case Channel::nnn:
			gf = &DecayRates::nnn;
			break;
		case Channel::nGAMMA:
			gf = &DecayRates::nGAMMA;
			break;
		case Channel::nEE:
			gf = &DecayRates::nEE;
			break;
		case Channel::nEM:
			gf = &DecayRates::nEM;
			break;
		case Channel::nMM:
			gf = &DecayRates::nMM;
			break;
		case Channel::nET:
			gf = &DecayRates::nET;
			break;
		case Channel::nMT:
			gf = &DecayRates::nMT;
			break;
		case Channel::nPI0:
			gf = &DecayRates::nPI0;
			break;
		case Channel::EPI:
			gf = &DecayRates::EPI;
			break;
		case Channel::MPI:
			gf = &DecayRates::MPI;
			break;
		case Channel::TPI:
			gf = &DecayRates::TPI;
			break;
		case Channel::EKA:
			gf = &DecayRates::EKA;
			break;
		case Channel::MKA:
			gf = &DecayRates::MKA;
			break;
		case Channel::EKAx:
			gf = &DecayRates::EKAx;
			break;
		case Channel::MKAx:
			gf = &DecayRates::MKAx;
			break;
		case Channel::nRHO0:
			gf = &DecayRates::nRHO0;
			break;
		case Channel::ERHO:
			gf = &DecayRates::ERHO;
			break;
		case Channel::MRHO:
			gf = &DecayRates::MRHO;
			break;
		case Channel::nETA:
			gf = &DecayRates::nETA;
			break;
		case Channel::nETAi:
			gf = &DecayRates::nETAi;
			break;
		case Channel::nOMEGA:
			gf = &DecayRates::nOMEGA;
			break;
		case Channel::nPHI:
			gf = &DecayRates::nPHI;
			break;
		case Channel::ECHARM:
			gf = &DecayRates::ECHARM;
			break;
		case Channel::ExpALL:
			gf = &DecayRates::ExpALL;
			break;
		default:
			throw std::invalid_argument("Channel " + Channel::toString(chan) + " unknown");
	}

	return (this->*gf)(mix);
}

//Return Gamma_tot - Gamma of interest
//
double DecayRates::Other(Channel::Name chan, const Mixing &mix)
{
	return Gamma(Channel::ALL, mix) - Gamma(chan, mix);
}

//Return the branching ratio
//
double DecayRates::Branch(Channel::Name chan, const Mixing &mix)
{
	if (chan == Channel::ALL)
		return 1.;

	if (Gamma(chan, mix) <= 0.)
		return 0.;

	return Gamma(chan, mix)/Gamma(Channel::ALL, mix);
}

//total decay width
double DecayRates::Total(const Mixing &mix)
{
	auto decs = Channel::Decays();
	return std::accumulate(decs.begin(), decs.end(), 0.,
			[=](double sum, Channel::Name chan) { return sum + Gamma(chan, mix); });
	// check capture here
}

//special here
double DecayRates::ExpALL(const Mixing &mix)
{
	return (nEE(mix) + nEM(mix) + nMM(mix) +
		EPI(mix) + MPI(mix) + nPI0(mix) );
}

//individual decay channels
//all mixing factors are factorised out
double DecayRates::nGAMMA(const Mixing &mix)
{
	if (!_table.count(Channel::nGAMMA))
		_table[Channel::nGAMMA] = Const::GF2 * std::pow(_N.M(), 5)
					* (27.0 * Const::Aem) / (192.0 * 64.0 * Const::pi4);

	return (1.0 + _N.IsMajorana()) * _table[Channel::nGAMMA]
		* (mix.Ue(2) + mix.Um(2) + mix.Ut(2));
}

double DecayRates::nnn(const Mixing &mix)
{
	if (!_table.count(Channel::nnn))
		_table[Channel::nnn] = Const::GF2 * std::pow(_N.M(), 5) / (192.0 * Const::pi3);

	return (1.0 + _N.IsMajorana()) * _table[Channel::nnn]
		* (mix.Ue(2) + mix.Um(2) + mix.Ut(2));
}

//M_Sterile > 2 M_Electron (always)
double DecayRates::nEE(const Mixing &mix)
{
	if (!_table.count(Channel::nEE)) {
		auto gamma = NeutrinoLeptonAA(Const::MNeutrino, Const::MElectron);
		_table[Channel::nEE] = gamma.first;
		_table[Channel::nEE_o] = gamma.second;
	}

	return _table[Channel::nEE] * mix.Ue(2) + _table[Channel::nEE_o] * (mix.Um(2) + mix.Ut(2));
}

//M_Sterile > Const::MMuon + M_Electron
double DecayRates::nEM(const Mixing &mix)	//Antiparticle is Elec
{
	if (!_table.count(Channel::nEM)) {
		auto gamma = NeutrinoLeptonAB(Const::MNeutrino, Const::MElectron, Const::MMuon);
		_table[Channel::nEM] = gamma.first;
		_table[Channel::nEM_o] = gamma.second;
	}

	return _table[Channel::nEM] * mix.Ue(2) + _table[Channel::nEM_o] * mix.Um(2);
}

//M_Sterile > 2 Const::MMuon
double DecayRates::nMM(const Mixing &mix)
{
	if (!_table.count(Channel::nMM)) {
		auto gamma = NeutrinoLeptonAA(Const::MNeutrino, Const::MMuon);
		_table[Channel::nMM] = gamma.first;
		_table[Channel::nMM_o] = gamma.second;
	}

	return _table[Channel::nMM] * mix.Um(2) + _table[Channel::nMM_o] * (mix.Ue(2) + mix.Ut(2));
}

//M_Sterile > Const::MTau + M_Electron
double DecayRates::nET(const Mixing &mix)	//Antiparticle is Elec
{
	if (!_table.count(Channel::nET)) {
		auto gamma = NeutrinoLeptonAB(Const::MNeutrino, Const::MElectron, Const::MTau);
		_table[Channel::nET] = gamma.first;
		_table[Channel::nET_o] = gamma.second;
	}

	return _table[Channel::nET] * mix.Ue(2) + _table[Channel::nET_o] * mix.Ut(2);
}

//M_Sterile > Const::MTau + Const::MMuon
double DecayRates::nMT(const Mixing &mix)	//Antiparticle is Muon
{
	if (!_table.count(Channel::nMT)) {
		auto gamma = NeutrinoLeptonAB(Const::MNeutrino, Const::MMuon, Const::MTau);
		_table[Channel::nMT] = gamma.first;
		_table[Channel::nMT_o] = gamma.second;
	}

	return _table[Channel::nMT] * mix.Um(2) + _table[Channel::nMT_o] * mix.Ut(2);
}

//M_Sterile > Const::MPion0
double DecayRates::nPI0(const Mixing &mix)
{
	if (!_table.count(Channel::nPI0))
		_table[Channel::nPI0] = std::pow(Const::DPion, 2)
				* NeutrinoPseudoMeson(Const::MNeutrino, Const::MPion0);

	return _table[Channel::nPI0] * (mix.Ue(2) + mix.Um(2) + mix.Ut(2));
}

//M_Sterile > Const::MPion
double DecayRates::EPI(const Mixing &mix)
{
	if (!_table.count(Channel::EPI))
		_table[Channel::EPI] = std::pow(Const::U_ud * Const::DPion, 2)
				* LeptonPseudoMeson(Const::MElectron, Const::MPion);
	
	return _table[Channel::EPI] * mix.Ue(2);
}

//M_Sterile > Const::MPion + Const::MMuon
double DecayRates::MPI(const Mixing &mix)
{
	if (!_table.count(Channel::MPI))
		_table[Channel::MPI] = std::pow(Const::U_ud * Const::DPion, 2)
				* LeptonPseudoMeson(Const::MMuon, Const::MPion);
	
	return _table[Channel::MPI] * mix.Um(2);
}

//M_Sterile > Const::MTau + Const::MPion
double DecayRates::TPI(const Mixing &mix)
{
	if (!_table.count(Channel::TPI))
		_table[Channel::TPI] = std::pow(Const::U_ud * Const::DPion, 2)
				* LeptonPseudoMeson(Const::MTau, Const::MPion);
	
	return _table[Channel::TPI] * mix.Ut(2);
}

//M_Sterile > Const::MKaon + M_Electron
double DecayRates::EKA(const Mixing &mix)
{
	if (!_table.count(Channel::EKA))
		_table[Channel::EKA] = std::pow(Const::U_us * Const::DKaon, 2)
				* LeptonPseudoMeson(Const::MElectron, Const::MKaon);

	return _table[Channel::EKA] * mix.Ue(2);
}

//M_Sterile > Const::MKaon + Const::MMuon
double DecayRates::MKA(const Mixing &mix)
{
	if (!_table.count(Channel::MKA))
		_table[Channel::MKA] = std::pow(Const::U_us * Const::DKaon, 2)
				* LeptonPseudoMeson(Const::MMuon, Const::MKaon);

	return _table[Channel::MKA] * mix.Um(2);
}

//M_Sterile > Const::MRho
double DecayRates::nRHO0(const Mixing &mix)
{
	if (!_table.count(Channel::nRHO0))
		_table[Channel::nRHO0] = std::pow(Const::DRho * Const::VRho, 2)
				* NeutrinoVectorMeson(Const::MNeutrino, Const::MRho0);

	return _table[Channel::nRHO0] * (mix.Ue(2) + mix.Um(2) + mix.Ut(2));
}

//M_Sterile > Const::MRho + M_Electron 
double DecayRates::ERHO(const Mixing &mix)
{
	if (!_table.count(Channel::ERHO))
		_table[Channel::ERHO] = std::pow(Const::U_ud * Const::DRho, 2)
				* LeptonVectorMeson(Const::MElectron, Const::MRho);

	return _table[Channel::ERHO] * mix.Ue(2);
}

//M_Sterile > Const::MRho + Const::MMuon 
double DecayRates::MRHO(const Mixing &mix)
{
	if (!_table.count(Channel::MRHO))
		_table[Channel::MRHO] = std::pow(Const::U_ud * Const::DRho, 2)
				* LeptonVectorMeson(Const::MMuon, Const::MRho);

	return _table[Channel::MRHO] * mix.Um(2);
}

//M_Sterile > Const::MKaon* + M_Electron 
double DecayRates::EKAx(const Mixing &mix)
{
	if (!_table.count(Channel::EKAx))
		_table[Channel::EKAx] = std::pow(Const::U_us * Const::DKaonx, 2)
				* LeptonVectorMeson(Const::MElectron, Const::MKaonx);

	return _table[Channel::EKAx] * mix.Ue(2);
}

//M_Sterile > Const::MKaon* + Const::MMuon 
double DecayRates::MKAx(const Mixing &mix)
{
	if (!_table.count(Channel::MKAx))
		_table[Channel::MKAx] = std::pow(Const::U_us * Const::DKaonx, 2)
				* LeptonVectorMeson(Const::MMuon, Const::MKaonx);

	return _table[Channel::MKAx] * mix.Um(2);
}

//M_Sterile > Const::MEta
double DecayRates::nETA(const Mixing &mix)
{
	if (!_table.count(Channel::nETA))
		_table[Channel::nETA] = std::pow(Const::DEta, 2)
				* NeutrinoPseudoMeson(Const::MNeutrino, Const::MEta);

	return _table[Channel::nETA] * (mix.Ue(2) + mix.Um(2) + mix.Ut(2));
}

//M_Sterile > Const::MEta'
double DecayRates::nETAi(const Mixing &mix)
{
	if (!_table.count(Channel::nETAi))
		_table[Channel::nETAi] = std::pow(Const::DEta, 2)
				* NeutrinoPseudoMeson(Const::MNeutrino, Const::MEtai);

	return _table[Channel::nETAi] * (mix.Ue(2) + mix.Um(2) + mix.Ut(2));
}

//M_Sterile > Const::MOmega 
double DecayRates::nOMEGA(const Mixing &mix)
{
	if (!_table.count(Channel::nOMEGA))
		_table[Channel::nOMEGA] = std::pow(Const::DOmega * Const::VOmega, 2)
				* NeutrinoVectorMeson(Const::MNeutrino, Const::MOmega);

	return _table[Channel::nOMEGA] * (mix.Ue(2) + mix.Um(2) + mix.Ut(2));
}

//M_Sterile > Const::MPhi 
double DecayRates::nPHI(const Mixing &mix)
{
	if (!_table.count(Channel::nPHI))
		_table[Channel::nPHI] = std::pow(Const::DOmega * Const::VPhi, 2)
				* NeutrinoVectorMeson(Const::MNeutrino, Const::MPhi);

	return _table[Channel::nPHI] * (mix.Ue(2) + mix.Um(2) + mix.Ut(2));
}

//M_Sterile > Const::MD + M_Electron
double DecayRates::ECHARM(const Mixing &mix)
{
	if (!_table.count(Channel::ECHARM))
		_table[Channel::ECHARM] = std::pow(Const::U_cd * Const::DCharm, 2)
				* LeptonPseudoMeson(Const::MElectron, Const::MD);

	return _table[Channel::ECHARM] * mix.Ue(2);
}

/////////////////
//Generic decay//
/////////////////
//
//CC version possible
double DecayRates::LeptonPseudoMeson(double m_lepton, double m_meson)
{
	_m_parent = _N.M();
	return I_LeptonPseudoMeson(std::pow(m_lepton / _m_parent, 2),
				   std::pow(m_meson / _m_parent, 2));
}

//integrated over angular dep.
double DecayRates::I_LeptonPseudoMeson(double x, double y)
{
	//double cos0 = theta < 0 ? 0.0 : cos(theta);
	//double fc = theta < 0.0 ? 2.0 : 1.0;

	return dGammad0_2B(x, y) * M2_LeptonPseudoMeson(0., x, y);
}

double DecayRates::NeutrinoPseudoMeson(double m_neut, double m_meson)
{
	_m_parent = _N.M();
	return I_NeutrinoPseudoMeson(std::pow(m_neut / _m_parent, 2),
				     std::pow(m_meson / _m_parent, 2));
}

//integrated over angular dep.
double DecayRates::I_NeutrinoPseudoMeson(double x, double y)
{
	//double cos0 = theta < 0 ? 0.0 : cos(theta);
	//double fc = theta < 0.0 ? 2.0 : 1.0;

	return dGammad0_2B(x, y) * M2_NeutrinoPseudoMeson(0., x, y);
}

//CC version possible
double DecayRates::LeptonVectorMeson(double m_lepton, double m_meson)
{
	_m_parent = _N.M();
	return I_LeptonVectorMeson(std::pow(m_lepton / _m_parent, 2),
				   std::pow(m_meson / _m_parent, 2));
}

//integrated over angular dep.
double DecayRates::I_LeptonVectorMeson(double x, double y)
{
	//double cos0 = theta < 0 ? 0.0 : cos(theta);
	//double fc = theta < 0.0 ? 2.0 : 1.0;

	return dGammad0_2B(x, y) * M2_LeptonVectorMeson(0., x, y);
}

double DecayRates::NeutrinoVectorMeson(double m_neut, double m_meson)
{
	_m_parent = _N.M();
	return  I_NeutrinoVectorMeson(std::pow(m_neut / _m_parent, 2),
				      std::pow(m_meson / _m_parent, 2));
}

//integrated over angular dep.
double DecayRates::I_NeutrinoVectorMeson(double x, double y)
{
	//double cos0 = theta < 0 ? 0.0 : cos(theta);
	//double fc = theta < 0.0 ? 2.0 : 1.0;

	return dGammad0_2B(x, y) * M2_NeutrinoVectorMeson(0., x, y);
}

std::pair<double, double> DecayRates::NeutrinoLeptonAA(double m_neut, double m_lepton)
{
	_m_parent = _N.M();
	double dMn2 = std::pow(m_neut / _m_parent, 2);
	double dML2 = std::pow(m_lepton / _m_parent, 2);

	//neutrino flavour is the same as leptons -> both Z and W interaction
	//double gL_CC = -0.5 + Const::sin2W + (2*GetParticle() - 1);	//times U(lepton flavour)
	//double gL_CC = -0.5 + Const::sin2W + 1;	//times U(lepton flavour)
	//double gR_CC = Const::sin2W;

	//neutrino flavour is different from leptons -> Z interaction
	double gL = -0.5 + Const::sin2W;	//times U(neutrino flavour)
	double gR = Const::sin2W;

	return std::make_pair(NeutrinoLeptonLepton(dMn2, dML2, dML2, gL + 1., gR),
                              NeutrinoLeptonLepton(dMn2, dML2, dML2, gL,      gR));
}

//CC version also available
std::pair<double, double> DecayRates::NeutrinoLeptonAB(double m_neut,
					double m_leptonA, double m_leptonB)
{
	_m_parent = _N.M();
	double dMn2 = std::pow(m_neut / _m_parent, 2);
	double dMA2 = std::pow(m_leptonA / _m_parent, 2);
	double dMB2 = std::pow(m_leptonB / _m_parent, 2);

	//for CC
	//  gL = 1.0;
	//  gR = 0.0;

	return std::make_pair(NeutrinoLeptonLepton(dMn2, dMA2, dMB2, 1., 0.),
			      NeutrinoLeptonLepton(dMn2, dMB2, dMA2, 1., 0.));
}

double DecayRates::NeutrinoLeptonLepton(double x, double y, double z, double gL, double gR)
{
	return dGammad2_3B() * I_NeutrinoLeptonLepton(x, y, z, gL, gR);
}

double DecayRates::I_NeutrinoLeptonLepton(double x, double y, double z, double gL, double gR)//, double theta)
{
	_vars = {x, y, z, gL, gR};
	//F_var.push_back(theta);	//3

	auto func = std::bind(&DecayRates::F_NeutrinoLeptonLepton_s, this, std::placeholders::_1);
	return Integration::Boole<1, double>(func); 
}

double DecayRates::F_NeutrinoLeptonLepton_s(double s)
{
	double x  = _vars[0];
	double y  = _vars[1];
	double z  = _vars[2];
	double gL = _vars[3];
	double gR = _vars[4];


	double t = s, u = s;
	double fcs = Limit(s, x, y, z);
	double fct = Limit(t, z, x, y);
	double fcu = Limit(u, y, z, x);

	if (gR > 0.)
		return gL * gL * fcs * M2_WW(s, 0, x, y, z) +
		       gR * gR * fct * M2_WW(t, 0, x, y, z) +
	   	   2 * gL * gR * fcu * M2_WZ(u, 0, x, y, z);
	return gL * gL * fcs * M2_WW(s, 0, x, y, z);
}
