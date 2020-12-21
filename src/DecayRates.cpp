#include "DecayRates.h"

DecayRates::DecayRates(const Neutrino &N) : Amplitude(N)
{
}

double DecayRates::MassThreshold(Channel::Name chan)
{
	if (_channel != chan) {
		LoadMass(chan);
		_channel = chan;
	}

	return std::accumulate(_masspdg.begin(), _masspdg.end(), 0.,
			[](double sum, const std::pair<double, int> mp) {
				return sum + mp.first; });
}

bool DecayRates::IsAllowed(Channel::Name chan)
{
	return (_N.M() >= MassThreshold(chan));
}

//return the decay width (Gamma)
//
double DecayRates::Gamma(Channel::Name chan, std::array<double, 3> mix) {
	return Gamma(chan, mix[0], mix[1], mix[2])
}

double DecayRates::Gamma(Channel::Name chan, double ue, double um, double ut)
{
	if (!IsAllowed(chan))
		return 0.;

	using GammaF = double (DecayRates::*)(double, double, double);
	GammaF gf;

	switch(chan)
	{
		case Channel::_ALL:
			gf = &DecayRates::Total;
			break;
		case Channel::_nnn:
			gf = &DecayRates::nnn;
			break;
		case Channel::_nGAMMA:
			gf = &DecayRates::nGAMMA;
			break;
		case Channel::_nEE:
			gf = &DecayRates::nEE;
			break;
		case Channel::_nEM:
			gf = &DecayRates::nEM;
			break;
		case Channel::_nMM:
			gf = &DecayRates::nMM;
			break;
		case Channel::_nET:
			gf = &DecayRates::nET;
			break;
		case Channel::_nMT:
			gf = &DecayRates::nMT;
			break;
		case Channel::_nPI0:
			gf = &DecayRates::nPI0;
			break;
		case Channel::_EPI:
			gf = &DecayRates::EPI;
			break;
		case Channel::_MPI:
			gf = &DecayRates::MPI;
			break;
		case Channel::_TPI:
			gf = &DecayRates::TPI;
			break;
		case Channel::_EKA:
			gf = &DecayRates::EKA;
			break;
		case Channel::_MKA:
			gf = &DecayRates::MKA;
			break;
		case Channel::_EKAx:
			gf = &DecayRates::EKAx;
			break;
		case Channel::_MKAx:
			gf = &DecayRates::MKAx;
			break;
		case Channel::_nRHO0:
			gf = &DecayRates::nRHO0;
			break;
		case Channel::_ERHO:
			gf = &DecayRates::ERHO;
			break;
		case Channel::_MRHO:
			gf = &DecayRates::MRHO;
			break;
		case Channel::_nETA:
			gf = &DecayRates::nETA;
			break;
		case Channel::_nETAi:
			gf = &DecayRates::nETAi;
			break;
		case Channel::_nOMEGA:
			gf = &DecayRates::nOMEGA;
			break;
		case Channel::_nPHI:
			gf = &DecayRates::nPHI;
			break;
		case Channel::_ECHARM:
			gf = &DecayRates::ECHARM;
			break;
		case Channel::_ExpALL:
			gf = &DecayRates::ExpALL;
			break;
		default:
			throw std::invalid_argument("Channel " + toString(chan) + " unknown");
	}

	return (this->*gf)(ue, um, ut);
}

//Return Gamma_tot - Gamma of interest
//
double DecayRates::Other(Channel::Name chan, const std::array<double, 3> &mix) {
	return Other(chan, mix[0], mix[1], mix[2]);
}

double DecayRates::Other(Channel::Name chan, ue, um, ut)
{
	return Gamma(Channel::_ALL, ue, um, ut) - Gamma(chan, ue, um, ut);
}

//Return the branching ration
//
double DecayRates::Branch(Channel::Name chan, const std::array<double, 3> &mix) {
	if (chan == Channels::_ALL)
		return 1.;	// should save some time..
	return Branch(chan, ue, um, ut);
}

double DecayRates::Branch(Channel::Name chan, double ue, double um, double ut)
{
	if (Gamma(chan, ue, um, ut) <= 0.)
		return 0.;
	else
		return Gamma(chan, ue, um, ut)/Gamma(Channel::_ALL, ue, um, ut);
}

//total decay width
double DecayRates::Total(double ue, double um, double ut)
{
	return std::accumulate(std::begin(Channel::Decays), std::end(Channel::Decays), 0.,
			[=](double sum, Channel::Name chan) { return sum + Gamma(chan, ue, um, ut); });
	// check capture here
}

//special here
double DecayRates::ExpALL(double ue, double um, double ut)
{
	return (nEE(ue, um, ut) + nEM(ue, um, ut) + nMM(ue, um, ut) +
		EPI(ue, um, ut) + MPI(ue, um, ut) + nPI0(ue, um, ut) );
}

//individual decay channels
//all mixing factors are factorised out
double DecayRates::nGAMMA(double ue, double um, double ut)
{
	if (!fdecay.count(Channel::_nGAMMA))
		fdecay[Channel::_nGAMMA] = Const::GF2 * std::pow(_N.M(), 5)
					* (27.0 * Const::Aem) / (192.0 * 64.0 * Const::pi4);

	return (1.0 + !GetFermion()) * fdecay[Channel::_nGAMMA] * (ue*ue + um*um + ut*ut);
}

double DecayRates::nnn(double ue, double um, double ut)
{
	if (!fdecay.count(Channel::_nnn))
		fdecay[Channel::_nnn] = Const::GF2 * std::pow(_N.M(), 5) / (192.0 * Const::pi3);

	return (1.0 + !GetFermion()) * fdecay[Channel::_nnn] * (ue*ue + um*um + ut*ut);
}

//M_Sterile > 2 M_Electron (always)
double DecayRates::nEE(double ue, double um, double ut)
{
	if (!fdecay.count(Channel::_nEE)) {
		auto gamma = NeutrinoLeptonAA(Const::MNeutrino, Const::MElectron);
		fdecay[Channel::_nEE] = gamma.first;
		fdecay[Channel::_nEE_o] = gamma.second,
	}

	return fdecay[Channel::_nEE] * ue*ue + fdecay[Channel::_nEE_o] * (um*um + ut*ut);
}

//M_Sterile > Const::MMuon + M_Electron
double DecayRates::nEM(double ue, double um, double ut)	//Antiparticle is Elec
{
	if (!fdecay.count(Channel::_nEM)) {
		auto gamma NeutrinoLeptonAB(Const::MNeutrino, Const::MElectron, Const::MMuon);
		fdecay[Channel::_nEM] = gamma.first;
		fdecay[Channel::_nEM_o] = gamma.second;
	}

	return fdecay[Channel::_nEM] * ue*ue + fdecay[Channel::_nEM_o] * um*um;
}

//M_Sterile > 2 Const::MMuon
double DecayRates::nMM(double ue, double um, double ut)
{
	if (!fdecay.count(Channel::_nMM)) {
		auto gamma = NeutrinoLeptonAA(Const::MNeutrino, Const::MMuon);
		fdecay[Channel::_nMM] = gamma.first;
		fdecay[Channel::_nMM_o] = gamma.second;
	}

	return fdecay[Channel::_nMM] * um*um + fdecay[Channel::_nMM_o] * (ue*ue + ut*ut);
}

//M_Sterile > Const::MTau + M_Electron
double DecayRates::nET(double ue, double um, double ut)	//Antiparticle is Elec
{
	if (!fdecay.count(Channel::_nET)) {
		auto gamma = NeutrinoLeptonAB(Const::MNeutrino, Const::MElectron, Const::MTau);
		fdecay[Channel::_nET] = gamma.first;
		fdecay[Channel::_nET_o] = gamma.second;
	}

	return fdecay[Channel::_nET] * ue*ue + fdecay[Channel::_nET_o] * ut*ut;
}

//M_Sterile > Const::MTau + Const::MMuon
double DecayRates::nMT(double ue, double um, double ut)	//Antiparticle is Muon
{
	if (!fdecay.count(Channel::_nMT)) {
		auto gamma = NeutrinoLeptonAB(Const::MNeutrino, Const::MMuon, Const::MTau);
		fdecay[Channel::_nMT] = gamma.first;
		fdecay[Channel::_nMT_o] = gamma.second;
	}

	return fdecay[Channel::_nMT] * um*um + fdecay[Channel::_nMT_o] * ut*ut;
}

//M_Sterile > Const::MPion0
double DecayRates::nPI0(double ue, double um, double ut)
{
	if (!fdecay.count(Channel::_nPI0))
		fdecay[Channel::_nPI0] = std::pow(Const::DPion, 2)
				* NeutrinoPseudoMeson(Const::MNeutrino, Const::MPion0);

	return fdecay[Channel::_nPI0] * (ue*ue + um*um + ut*ut);
}

//M_Sterile > Const::MPion
double DecayRates::EPI(double ue, double um, double ut)
{
	if (!fdecay.count(Channel::_EPI))
		fdecay[Channel::_EPI] = std::pow(Const::U_ud * Const::DPion, 2)
				* LeptonPseudoMeson(Const::MElectron, Const::MPion);
	}
	
	return fdecay[Channel::_EPI] * ue*ue;
}

//M_Sterile > Const::MPion + Const::MMuon
double DecayRates::MPI(double ue, double um, double ut)
{
	if (!fdecay.count(Channel::_MPI))
		fdecay[Channel::_MPI] = std::pow(Const::U_ud * Const::DPion, 2)
				* LeptonPseudoMeson(Const::MMuon, Const::MPion);
	}
	
	return fdecay[Channel::_MPI] * um*um;
}

//M_Sterile > Const::MTau + Const::MPion
double DecayRates::TPI(double ue, double um, double ut)
{
	if (!fdecay.count(Channel::_TPI))
		fdecay[Channel::_TPI] = std::pow(Const::U_ud * Const::DPion, 2)
				* LeptonPseudoMeson(Const::MTau, Const::MPion);
	}
	
	return fdecay[Channel::_TPI] * ut*ut;
}

//M_Sterile > Const::MKaon + M_Electron
double DecayRates::EKA(double ue, double um, double ut)
{
	if (!fdecay.count(Channel::_EKA))
		fdecay[Channel::_EKA] = std::pow(Const::U_us * Const::DKaon, 2)
				* LeptonPseudoMeson(Const::MElectron, Const::MKaon);

	return fdecay[Channel::_EKA] * ue*ue;
}

//M_Sterile > Const::MKaon + Const::MMuon
double DecayRates::MKA(double ue, double um, double ut)
{
	if (!fdecay.count(Channel::_MKA))
		fdecay[Channel::_MKA] = std::pow(Const::U_us * Const::DKaon, 2)
				* LeptonPseudoMeson(Const::MMuon, Const::MKaon);

	return fdecay[Channel::_MKA] * um*um;
}

//M_Sterile > Const::MRho
double DecayRates::nRHO0(double ue, double um, double ut)
{
	if (!fdecay.count(Channel::_nRHO))
		fdecay[Channel::_nRHO] = std::pow(Const::DRho * Const::VRho, 2)
				* NeutrinoVectorMeson(Const::MNeutrino, Const::MRho0);

			

	return fdecay[Channel::_nRHO] * (ue*ue + um*um + ut*ut);
}

//M_Sterile > Const::MRho + M_Electron 
double DecayRates::ERHO(double ue, double um, double ut)
{
	if (!fdecay.count(Channel::_ERHO))
		fdecay[Channel::_ERHO] = std::pow(Const::U_ud * Const::DRho, 2)
				* LeptonVectorMeson(Const::MElectron, Const::MRho);

	return fdecay[Channel::_ERHO] * ue*ue;
}

//M_Sterile > Const::MRho + Const::MMuon 
double DecayRates::MRHO(double ue, double um, double ut)
{
	if (!fdecay.count(Channel::_MRHO))
		fdecay[Channel::_MRHO] = std::pow(Const::U_ud * Const::DRho, 2)
				* LeptonVectorMeson(Const::MMuon, Const::MRho);

	return fdecay[Channel::_MRHO] * um*um;
}

//M_Sterile > Const::MKaon* + M_Electron 
double DecayRates::EKAx(double ue, double um, double ut)
{
	if (fdecay.count(Channel::_EKAx))
		fdecay[Channel::_EKAx] = std::pow(Const::U_us * Const::DKaonx, 2)
				* LeptonVectorMeson(Const::MElectron, Const::MKaonx);

	return fdecay[Channel::_EKAx] * ue*ue;
}

//M_Sterile > Const::MKaon* + Const::MMuon 
double DecayRates::MKAx(double ue, double um, double ut)
{
	if (fdecay.count(Channel::_MKAx))
		fdecay[Channel::_MKAx] = std::pow(Const::U_us * Const::DKaonx, 2)
				* LeptonVectorMeson(Const::MMuon, Const::MKaonx);

	return fdecay[Channel::_MKAx] * um*um;
}

//M_Sterile > Const::MEta
double DecayRates::nETA(double ue, double um, double ut)
{
	if (!fdecay.count(Channel::_nETA))
		fdecay[Channel::_nETA] = std::pow(Const::DEta, 2)
				* NeutrinoPseudoMeson(Const::MNeutrino, Const::MEta);

	return fdecay[Channel::_nETA] * (ue*ue + um*um + ut*ut);
}

//M_Sterile > Const::MEta'
double DecayRates::nETAi(double ue, double um, double ut)
{
	if (!fdecay.count(Channel::_nETAi))
		fdecay[Channel::_nETAi] = std::pow(Const::DEta, 2)
				* NeutrinoPseudoMeson(Const::MNeutrino, Const::MEtai);

	return fdecay[Channel::_nETAi] * (ue*ue + um*um + ut*ut);
}

//M_Sterile > Const::MOmega 
double DecayRates::nOMEGA(double ue, double um, double ut)
{
	if (!fdecay.count(Channel::_nOMEGA))
		fdecay[Channel::_nOMEGA] = std::pow(Const::DOmega * Const::VOmega, 2)
				* NeutrinoVectorMeson(Const::MNeutrino, Const::MOmega);

	return fdecay[Channel::_nOMEGA] * (ue*ue + um*um + ut*ut);
}

//M_Sterile > Const::MPhi 
double DecayRates::nPHI(double ue, double um, double ut)
{
	if (!fdecay.count(Channel::_nPHI))
		fdecay[Channel::_nPHI] = std::pow(Const::DOmega * Const::VPhi, 2)
				* NeutrinoVectorMeson(Const::MNeutrino, Const::MPhi);

	return fdecay[Channel::_nPHI] * (ue*ue + um*um + ut*ut);
}

//M_Sterile > Const::MD + M_Electron
double DecayRates::ECHARM(double ue, double um, double ut)
{
	if (!fdecay.count(Channel::_ECHARM))
		fdecay[Channel::_ECHARM] = std::pow(Const::U_cd * Const::DCharm, 2)
				* LeptonPseudoMeson(Const::MElectron, Const::MD);

	return fdecay[Channel::_ECHARM] * ue*ue;
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

	double M2 = M2_LeptonPseudoMeson(0, x, y);
	return dGammad0_2B(M2, x, y);
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

	double M2 = M2_NeutrinoPseudoMeson(0, x, y);
	return dGammad0_2B(M2, x, y);
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

	double M2 = M2_LeptonVectorMeson(0, x, y);
	return dGammad0_2B(M2, x, y);
}

double DecayRates::NeutrinoVectorMeson(double m_neut, double m_meson)
{
	_m_parent = _N.M()
	return  I_NeutrinoVectorMeson(std::pow(m_neut / _m_parent, 2),
				      std::pow(m_meson / _m_parent, 2));
}

//integrated over angular dep.
double DecayRates::I_NeutrinoVectorMeson(double x, double y)
{
	//double cos0 = theta < 0 ? 0.0 : cos(theta);
	//double fc = theta < 0.0 ? 2.0 : 1.0;

	double M2 = M2_NeutrinoVectorMeson(0, x, y);
	return dGammad0_2B(M2, x, y);
}

std::pair<double, double> DecayRates::NeutrinoLeptonAA(double m_neut, double m_lepton)
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

	return std::make_pair(NeutrinoLeptonLepton(dMn2, dML2, dML2, gL_CC, gR_CC),
                              NeutrinoLeptonLepton(dMn2, dML2, dML2, gL_NC, gR_NC));
}

//CC version also available
std::pair<double, double. DecayRates::NeutrinoLeptonAB(double m_neut,
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
	double M2 = I_NeutrinoLeptonLepton(x, y, z, gL, gR);
	return dGammad2_3B(M2);
}

double DecayRates::I_NeutrinoLeptonLepton(double x, double y, double z, double gL, double gR)//, double theta)
{
	_vars = {x, y, z, gL, gR};
	//F_var.push_back(theta);	//3

	return BooleIntegration(&DecayRates::I_NeutrinoLeptonLepton_s); 
}

double DecayRates::I_NeutrinoLeptonLepton_s(double s)
{
	double x  = _vars[0];
	double y  = _vars[1];
	double z  = _vars[2];
	double gL = _vars[3];
	double gR = _vars[4];

	double s_ = s, t_ = s, u_ = s;
	double fcu = Limit(u_, y, z, x);
	double fct = Limit(t_, z, x, y);
	double fcs = Limit(s_, x, y, z);

	return     gL * gL * fcs * M2_WW(s_, 0, x, y, z) +
	           gR * gR * fct * M2_WW(t_, 0, x, y, z) +
	       2 * gL * gR * fcu * M2_WZ(u_, 0, x, y, z);
}

double Reset() {
	fdecay.clear();
}
