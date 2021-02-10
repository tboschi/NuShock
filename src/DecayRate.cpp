#include "physics/DecayRate.h"

//return the decay width (Gamma)
//
double DecayRate::Gamma(Decay::Channel chan, const Mixing &mix)
{
	if (!IsAllowed(chan))
		return 0.;

	using GammaF = double (DecayRate::*)(const Mixing &mix);
	GammaF gf;

	switch(chan)
	{
		//case Channel::ALL:
			//gf = &DecayRate::Total;
			//break;
		case Decay::Channel::nnn:
			gf = &DecayRate::nnn;
			break;
		case Decay::Channel::nGAMMA:
			gf = &DecayRate::nGAMMA;
			break;
		case Decay::Channel::nEE:
			gf = &DecayRate::nEE;
			break;
		case Decay::Channel::nEM:
			gf = &DecayRate::nEM;
			break;
		case Decay::Channel::nMM:
			gf = &DecayRate::nMM;
			break;
		case Decay::Channel::nET:
			gf = &DecayRate::nET;
			break;
		case Decay::Channel::nMT:
			gf = &DecayRate::nMT;
			break;
		case Decay::Channel::nPI0:
			gf = &DecayRate::nPI0;
			break;
		case Decay::Channel::EPI:
			gf = &DecayRate::EPI;
			break;
		case Decay::Channel::MPI:
			gf = &DecayRate::MPI;
			break;
		case Decay::Channel::TPI:
			gf = &DecayRate::TPI;
			break;
		case Decay::Channel::EKA:
			gf = &DecayRate::EKA;
			break;
		case Decay::Channel::MKA:
			gf = &DecayRate::MKA;
			break;
		case Decay::Channel::EKAx:
			gf = &DecayRate::EKAx;
			break;
		case Decay::Channel::MKAx:
			gf = &DecayRate::MKAx;
			break;
		case Decay::Channel::nRHO0:
			gf = &DecayRate::nRHO0;
			break;
		case Decay::Channel::ERHO:
			gf = &DecayRate::ERHO;
			break;
		case Decay::Channel::MRHO:
			gf = &DecayRate::MRHO;
			break;
		case Decay::Channel::nETA:
			gf = &DecayRate::nETA;
			break;
		case Decay::Channel::nETAi:
			gf = &DecayRate::nETAi;
			break;
		case Decay::Channel::nOMEGA:
			gf = &DecayRate::nOMEGA;
			break;
		case Decay::Channel::nPHI:
			gf = &DecayRate::nPHI;
			break;
		case Decay::Channel::EDs:
			gf = &DecayRate::EDs;
			break;
		case Decay::Channel::ExpALL:
			gf = &DecayRate::ExpALL;
			break;
		default:
			throw std::invalid_argument("Decay channel " + Decay::toString(chan)
						+ " is unknown");
	}

	return (this->*gf)(mix);
}

//Return Gamma_tot - Gamma of interest
//
double DecayRate::Other(Decay::Channel chan, const Mixing &mix)
{
	return Total(mix) - Gamma(chan, mix);
}

//Return the branching ratio
//
double DecayRate::Branch(Decay::Channel chan, const Mixing &mix)
{
	if (Gamma(chan, mix) <= 0.)
		return 0.;

	return Gamma(chan, mix)/Total(mix);
}

//total decay width
double DecayRate::Total(const Mixing &mix)
{
	auto decs = Decay::Channels();
	return std::accumulate(decs.begin(), decs.end(), 0.,
			[&](double sum, Decay::Channel chan) { return sum + Gamma(chan, mix); });
	// check capture here
}

//special here
double DecayRate::ExpALL(const Mixing &mix)
{
	auto decs = Decay::Detections();
	return std::accumulate(decs.begin(), decs.end(), 0.,
			[&](double sum, Decay::Channel chan) { return sum + Gamma(chan, mix); });
}

//individual decay channels
//all mixing factors are factorised out
double DecayRate::nGAMMA(const Mixing &mix)
{
	if (!_table.count(Decay::Channel::nGAMMA))
		_table[Decay::Channel::nGAMMA] = Const::GF2 * std::pow(_N.M(), 5)
					* (27.0 * Const::Aem) / (192.0 * 64.0 * Const::pi4);

	return (1.0 + _N.IsMajorana()) * _table[Decay::Channel::nGAMMA]
		* (mix.Ue(2) + mix.Um(2) + mix.Ut(2));
}

double DecayRate::nnn(const Mixing &mix)
{
	if (!_table.count(Decay::Channel::nnn))
		_table[Decay::Channel::nnn] = Const::GF2 * std::pow(_N.M(), 5) / (192.0 * Const::pi3);

	return (1.0 + _N.IsMajorana()) * _table[Decay::Channel::nnn]
		* (mix.Ue(2) + mix.Um(2) + mix.Ut(2));
}

//M_Sterile > 2 M_Electron (always)
double DecayRate::nEE(const Mixing &mix)
{
	if (!_table.count(Decay::Channel::nEE)) {
		auto gamma = NeutrinoLeptonAA(Const::MNeutrino, Const::MElectron);
		_table[Decay::Channel::nEE] = gamma.first;
		_table[Decay::Channel::nEE_o] = gamma.second;
	}

	return _table[Decay::Channel::nEE] * mix.Ue(2) + _table[Decay::Channel::nEE_o] * (mix.Um(2) + mix.Ut(2));
}

//M_Sterile > Const::MMuon + M_Electron
double DecayRate::nEM(const Mixing &mix)	//Antiparticle is Elec
{
	if (!_table.count(Decay::Channel::nEM)) {
		auto gamma = NeutrinoLeptonAB(Const::MNeutrino, Const::MElectron, Const::MMuon);
		_table[Decay::Channel::nEM] = gamma.first;
		_table[Decay::Channel::nEM_o] = gamma.second;
	}

	return _table[Decay::Channel::nEM] * mix.Ue(2) + _table[Decay::Channel::nEM_o] * mix.Um(2);
}

//M_Sterile > 2 Const::MMuon
double DecayRate::nMM(const Mixing &mix)
{
	if (!_table.count(Decay::Channel::nMM)) {
		auto gamma = NeutrinoLeptonAA(Const::MNeutrino, Const::MMuon);
		_table[Decay::Channel::nMM] = gamma.first;
		_table[Decay::Channel::nMM_o] = gamma.second;
	}

	return _table[Decay::Channel::nMM] * mix.Um(2) + _table[Decay::Channel::nMM_o] * (mix.Ue(2) + mix.Ut(2));
}

//M_Sterile > Const::MTau + M_Electron
double DecayRate::nET(const Mixing &mix)	//Antiparticle is Elec
{
	if (!_table.count(Decay::Channel::nET)) {
		auto gamma = NeutrinoLeptonAB(Const::MNeutrino, Const::MElectron, Const::MTau);
		_table[Decay::Channel::nET] = gamma.first;
		_table[Decay::Channel::nET_o] = gamma.second;
	}

	return _table[Decay::Channel::nET] * mix.Ue(2) + _table[Decay::Channel::nET_o] * mix.Ut(2);
}

//M_Sterile > Const::MTau + Const::MMuon
double DecayRate::nMT(const Mixing &mix)	//Antiparticle is Muon
{
	if (!_table.count(Decay::Channel::nMT)) {
		auto gamma = NeutrinoLeptonAB(Const::MNeutrino, Const::MMuon, Const::MTau);
		_table[Decay::Channel::nMT] = gamma.first;
		_table[Decay::Channel::nMT_o] = gamma.second;
	}

	return _table[Decay::Channel::nMT] * mix.Um(2) + _table[Decay::Channel::nMT_o] * mix.Ut(2);
}

//M_Sterile > Const::MPion0
double DecayRate::nPI0(const Mixing &mix)
{
	if (!_table.count(Decay::Channel::nPI0))
		_table[Decay::Channel::nPI0] = std::pow(Const::DPion, 2)
				* NeutrinoPseudoMeson(Const::MNeutrino, Const::MPion0);

	return _table[Decay::Channel::nPI0] * (mix.Ue(2) + mix.Um(2) + mix.Ut(2));
}

//M_Sterile > Const::MPion
double DecayRate::EPI(const Mixing &mix)
{
	if (!_table.count(Decay::Channel::EPI))
		_table[Decay::Channel::EPI] = std::pow(Const::U_ud * Const::DPion, 2)
				* LeptonPseudoMeson(Const::MElectron, Const::MPion);
	
	return _table[Decay::Channel::EPI] * mix.Ue(2);
}

//M_Sterile > Const::MPion + Const::MMuon
double DecayRate::MPI(const Mixing &mix)
{
	if (!_table.count(Decay::Channel::MPI))
		_table[Decay::Channel::MPI] = std::pow(Const::U_ud * Const::DPion, 2)
				* LeptonPseudoMeson(Const::MMuon, Const::MPion);
	
	return _table[Decay::Channel::MPI] * mix.Um(2);
}

//M_Sterile > Const::MTau + Const::MPion
double DecayRate::TPI(const Mixing &mix)
{
	if (!_table.count(Decay::Channel::TPI))
		_table[Decay::Channel::TPI] = std::pow(Const::U_ud * Const::DPion, 2)
				* LeptonPseudoMeson(Const::MTau, Const::MPion);
	
	return _table[Decay::Channel::TPI] * mix.Ut(2);
}

//M_Sterile > Const::MKaon + M_Electron
double DecayRate::EKA(const Mixing &mix)
{
	if (!_table.count(Decay::Channel::EKA))
		_table[Decay::Channel::EKA] = std::pow(Const::U_us * Const::DKaon, 2)
				* LeptonPseudoMeson(Const::MElectron, Const::MKaon);

	return _table[Decay::Channel::EKA] * mix.Ue(2);
}

//M_Sterile > Const::MKaon + Const::MMuon
double DecayRate::MKA(const Mixing &mix)
{
	if (!_table.count(Decay::Channel::MKA))
		_table[Decay::Channel::MKA] = std::pow(Const::U_us * Const::DKaon, 2)
				* LeptonPseudoMeson(Const::MMuon, Const::MKaon);

	return _table[Decay::Channel::MKA] * mix.Um(2);
}

//M_Sterile > Const::MRho
double DecayRate::nRHO0(const Mixing &mix)
{
	if (!_table.count(Decay::Channel::nRHO0))
		_table[Decay::Channel::nRHO0] = std::pow(Const::DRho * Const::VRho, 2)
				* NeutrinoVectorMeson(Const::MNeutrino, Const::MRho0);

	return _table[Decay::Channel::nRHO0] * (mix.Ue(2) + mix.Um(2) + mix.Ut(2));
}

//M_Sterile > Const::MRho + M_Electron 
double DecayRate::ERHO(const Mixing &mix)
{
	if (!_table.count(Decay::Channel::ERHO))
		_table[Decay::Channel::ERHO] = std::pow(Const::U_ud * Const::DRho, 2)
				* LeptonVectorMeson(Const::MElectron, Const::MRho);

	return _table[Decay::Channel::ERHO] * mix.Ue(2);
}

//M_Sterile > Const::MRho + Const::MMuon 
double DecayRate::MRHO(const Mixing &mix)
{
	if (!_table.count(Decay::Channel::MRHO))
		_table[Decay::Channel::MRHO] = std::pow(Const::U_ud * Const::DRho, 2)
				* LeptonVectorMeson(Const::MMuon, Const::MRho);

	return _table[Decay::Channel::MRHO] * mix.Um(2);
}

//M_Sterile > Const::MKaon* + M_Electron 
double DecayRate::EKAx(const Mixing &mix)
{
	if (!_table.count(Decay::Channel::EKAx))
		_table[Decay::Channel::EKAx] = std::pow(Const::U_us * Const::DKaonx, 2)
				* LeptonVectorMeson(Const::MElectron, Const::MKaonx);

	return _table[Decay::Channel::EKAx] * mix.Ue(2);
}

//M_Sterile > Const::MKaon* + Const::MMuon 
double DecayRate::MKAx(const Mixing &mix)
{
	if (!_table.count(Decay::Channel::MKAx))
		_table[Decay::Channel::MKAx] = std::pow(Const::U_us * Const::DKaonx, 2)
				* LeptonVectorMeson(Const::MMuon, Const::MKaonx);

	return _table[Decay::Channel::MKAx] * mix.Um(2);
}

//M_Sterile > Const::MEta
double DecayRate::nETA(const Mixing &mix)
{
	if (!_table.count(Decay::Channel::nETA))
		_table[Decay::Channel::nETA] = std::pow(Const::DEta, 2)
				* NeutrinoPseudoMeson(Const::MNeutrino, Const::MEta);

	return _table[Decay::Channel::nETA] * (mix.Ue(2) + mix.Um(2) + mix.Ut(2));
}

//M_Sterile > Const::MEta'
double DecayRate::nETAi(const Mixing &mix)
{
	if (!_table.count(Decay::Channel::nETAi))
		_table[Decay::Channel::nETAi] = std::pow(Const::DEta, 2)
				* NeutrinoPseudoMeson(Const::MNeutrino, Const::MEtai);

	return _table[Decay::Channel::nETAi] * (mix.Ue(2) + mix.Um(2) + mix.Ut(2));
}

//M_Sterile > Const::MOmega 
double DecayRate::nOMEGA(const Mixing &mix)
{
	if (!_table.count(Decay::Channel::nOMEGA))
		_table[Decay::Channel::nOMEGA] = std::pow(Const::DOmega * Const::VOmega, 2)
				* NeutrinoVectorMeson(Const::MNeutrino, Const::MOmega);

	return _table[Decay::Channel::nOMEGA] * (mix.Ue(2) + mix.Um(2) + mix.Ut(2));
}

//M_Sterile > Const::MPhi 
double DecayRate::nPHI(const Mixing &mix)
{
	if (!_table.count(Decay::Channel::nPHI))
		_table[Decay::Channel::nPHI] = std::pow(Const::DOmega * Const::VPhi, 2)
				* NeutrinoVectorMeson(Const::MNeutrino, Const::MPhi);

	return _table[Decay::Channel::nPHI] * (mix.Ue(2) + mix.Um(2) + mix.Ut(2));
}

//M_Sterile > Const::MD + M_Electron
double DecayRate::EDs(const Mixing &mix)
{
	if (!_table.count(Decay::Channel::EDs))
		_table[Decay::Channel::EDs] = std::pow(Const::U_cd * Const::DCharm, 2)
				* LeptonPseudoMeson(Const::MElectron, Const::MD);

	return _table[Decay::Channel::EDs] * mix.Ue(2);
}

/////////////////
//Generic decay//
/////////////////
//
//CC version possible
double DecayRate::LeptonPseudoMeson(double m_lepton, double m_meson)
{
	_m_parent = _N.M();
	double x = std::pow(m_lepton / _m_parent, 2);
        double y = std::pow(m_meson / _m_parent, 2);

	// integrated over angular dependencies
	return dGammad0_2B(x, y) * M2_LeptonPseudoMeson(0., x, y);
}

/*
//integrated over angular dep.
double DecayRate::I_LeptonPseudoMeson(double x, double y)
{
	//double cos0 = theta < 0 ? 0.0 : cos(theta);
	//double fc = theta < 0.0 ? 2.0 : 1.0;

}
*/

double DecayRate::NeutrinoPseudoMeson(double m_neut, double m_meson)
{
	_m_parent = _N.M();
	double x = std::pow(m_neut / _m_parent, 2);
        double y = std::pow(m_meson / _m_parent, 2);

	// integrated over angular dependencies
	return dGammad0_2B(x, y) * M2_NeutrinoPseudoMeson(0., x, y);
}

/*
//integrated over angular dep.
double DecayRate::I_NeutrinoPseudoMeson(double x, double y)
{
	//double cos0 = theta < 0 ? 0.0 : cos(theta);
	//double fc = theta < 0.0 ? 2.0 : 1.0;

	return dGammad0_2B(x, y) * M2_NeutrinoPseudoMeson(0., x, y);
}
*/

//CC version possible
double DecayRate::LeptonVectorMeson(double m_lepton, double m_meson)
{
	_m_parent = _N.M();
	double x = std::pow(m_lepton / _m_parent, 2);
        double y = std::pow(m_meson / _m_parent, 2);

	//integrated over angular dep.
	return dGammad0_2B(x, y) * M2_LeptonVectorMeson(0., x, y);
}

/*
//integrated over angular dep.
double DecayRate::I_LeptonVectorMeson(double x, double y)
{
	//double cos0 = theta < 0 ? 0.0 : cos(theta);
	//double fc = theta < 0.0 ? 2.0 : 1.0;

	return dGammad0_2B(x, y) * M2_LeptonVectorMeson(0., x, y);
}
*/

double DecayRate::NeutrinoVectorMeson(double m_neut, double m_meson)
{
	_m_parent = _N.M();
	double x = std::pow(m_neut / _m_parent, 2);
        double y = std::pow(m_meson / _m_parent, 2);

	//integrated over angular dep.
	return dGammad0_2B(x, y) * M2_NeutrinoVectorMeson(0., x, y);
}

/*
double DecayRate::I_NeutrinoVectorMeson(double x, double y)
{
	//double cos0 = theta < 0 ? 0.0 : cos(theta);
	//double fc = theta < 0.0 ? 2.0 : 1.0;

	return dGammad0_2B(x, y) * M2_NeutrinoVectorMeson(0., x, y);
}
*/

std::pair<double, double> DecayRate::NeutrinoLeptonAA(double m_neut, double m_lepton)
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
std::pair<double, double> DecayRate::NeutrinoLeptonAB(double m_neut,
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

double DecayRate::NeutrinoLeptonLepton(double x, double y, double z, double gL, double gR)
{
	_vars = {x, y, z, gL, gR};

	auto func = std::bind(&DecayRate::NeutrinoLeptonLepton_s, this, std::placeholders::_1);
	return dGammad2_3B() * Integration::Boole<1, double>(func); 
}

// integrand
double DecayRate::NeutrinoLeptonLepton_s(double s)
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

void DecayRate::Reset() {
	_table.clear();
}
