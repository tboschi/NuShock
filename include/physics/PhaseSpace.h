/*
 * PhaseSpace generator for simulation of sterile neutrino decays (it is not needed for production)
 * It also returns correct 
 *
 * Author: Tommaso Boschi
 */

#ifndef PHASESPACE_H
#define PHASESPACE_H

// base class
#include "physics/Amplitude.h"

//ROOT include
#include "TGenPhaseSpace.h"

#include <random>
#include "physics/Channels.h"
#include "physics/Particle.h"
#include "physics/Neutrino.h"

#include "tools/Optimization.h"
#include "tools/RNG.h"

class PhaseSpace : public Amplitude
{
	public:

		PhaseSpace(Neutrino N = Neutrino());
		~PhaseSpace();

		std::pair<std::vector<Particle>, double> Generate(Channel::Name chan, const Particle &part,
								const Mixing &mix = Mixing());
		std::pair<std::vector<Particle>, double> Generate(Channel::Name chan, const TLorentzVector &frame,
								const Mixing &mix = Mixing());

	private:
		bool SetDecay(Channel::Name chan);
		double Gamma(Channel::Name chan, const Mixing &mix = Mixing());

		double nEE(const Mixing &mix = Mixing());
		double nMM(const Mixing &mix = Mixing());
		std::pair<double, double> NeutrinoLeptonAA(Channel::Name chan, Channel::Name chan_o, double m_neut, double m_lepton);

		double nEM(const Mixing &mix = Mixing());
		double nET(const Mixing &mix = Mixing());
		double nMT(const Mixing &mix = Mixing());
		std::pair<double, double> NeutrinoLeptonAB(Channel::Name chan, Channel::Name chan_o, double m_neut, double m_leptonA, double m_leptonB);

		double NeutrinoLeptonLepton_max(double x, double y, double z, double gL, double gR);
		double F_NeutrinoLeptonLepton_max(double p[]);
		double NeutrinoLeptonLepton(double s, double t, double cos0, double cos1, double x, double y, double z, double gL, double gR);



		double nPI0(const Mixing &mix = Mixing());
		double nETA(const Mixing &mix = Mixing());
		double nETAi(const Mixing &mix = Mixing());
		double NeutrinoPseudoMeson(Channel::Name chan, double m_neut, double m_meson);

		double EPI(const Mixing &mix = Mixing());
		double MPI(const Mixing &mix = Mixing());
		double TPI(const Mixing &mix = Mixing());
		double EKA(const Mixing &mix = Mixing());
		double MKA(const Mixing &mix = Mixing());
		double ECHARM(const Mixing &mix = Mixing());
		double LeptonPseudoMeson(Channel::Name chan, double m_lepton, double m_meson);

		double ToPseudoMeson_max(double x, double y);
		double ToPseudoMeson_cos0_max(double cos0);
		double ToPseudoMeson(double cos0, double x, double y);



		double nRHO0(const Mixing &mix = Mixing());
		double nOMEGA(const Mixing &mix = Mixing());
		double nPHI(const Mixing &mix = Mixing());
		double NeutrinoVectorMeson(Channel::Name chan, double m_neut, double m_meson);

		double ERHO(const Mixing &mix = Mixing());
		double MRHO(const Mixing &mix = Mixing());
		double EKAx(const Mixing &mix = Mixing());
		double MKAx(const Mixing &mix = Mixing());
		double LeptonVectorMeson(Channel::Name chan, double M_Neut, double M_Meson);

		double ToVectorMeson_max(double x, double y);
		double ToVectorMeson_cos0_max(double cos0);
		double ToVectorMeson(double cos0, double x, double y);



		double MuonE(const Mixing &mix = Mixing());
		double MuonM(const Mixing &mix = Mixing());
		double TauEE(const Mixing &mix = Mixing());
		double TauET(const Mixing &mix = Mixing());
		double TauMM(const Mixing &mix = Mixing());
		double TauMT(const Mixing &mix = Mixing());

		double LeptonNeutrino(Channel::Name chan, double m_lepton0, double m_lepton, double m_neut = 0.);
		double LeptonNeutrino_max(double x, double y, double z);
		double LeptonNeutrino_u_max(double u);

		double AntileptonNeutrino(Channel::Name chan, double m_lepton0, double m_lepton, double m_neut = 0.);
		double AntileptonNeutrino_max(double x, double y, double z);
		double AntileptonNeutrino_s_max(double s);



		double TauPI(const Mixing &mix = Mixing());
		//double Tau2PI(const Mixing &mix = Mixing());

		// ratio could be 1
		double LeptonMeson(Channel::Name chan, double m_lepton, double m_meson);
		double LeptonMeson_max(double m_lepton, double m_meson);


		double PionE(const Mixing &mix = Mixing());
		double PionM(const Mixing &mix = Mixing());
		double KaonE(const Mixing &mix = Mixing());
		double KaonM(const Mixing &mix = Mixing());
		double CharmE(const Mixing &mix = Mixing());
		double CharmM(const Mixing &mix = Mixing());
		double CharmT(const Mixing &mix = Mixing());

		// ratio could be 1
		double MesonTwo(Channel::Name chan, double m_meson, double m_lepton);
		double MesonTwo_max(double m_meson, double m_lepton);



		double Kaon0E(const Mixing &mix = Mixing());
		double Kaon0M(const Mixing &mix = Mixing());
		double KaonCE(const Mixing &mix = Mixing());
		double KaonCM(const Mixing &mix = Mixing());

		double MesonThree(Channel::Name chan, double m_meson0, double m_meson, double m_lepton,
				  double L_, double L0);
		double MesonThree_max(double x, double y, double z, double L_, double L0);
		double F_MesonThree_max(double p[]);



		//KINEMATICS

		double Kinematic_2B();
		std::array<double, 6> Kinematic_3B();

	private:
		TGenPhaseSpace *_genps;

		double (Amplitude::*_M2_F)(double, double, double);

		// reverse sign to find minimum
		double (PhaseSpace::*_to_maximize)(double []);
};

#endif
