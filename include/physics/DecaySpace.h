/*
 * PhaseSpace generator for simulation of sterile neutrino decays (it is not needed for production)
 * It also returns correct 
 *
 * Author: Tommaso Boschi
 */

#ifndef DECAYSPACE_H
#define DECAYSPACE_H

// base class
#include "physics/Amplitude.h"
#include "physics/PhaseSpace.h"

//ROOT include
#include "TGenPhaseSpace.h"

#include <unordered_map>
#include <random>

#include "physics/Decays.h"
#include "physics/Particle.h"
#include "physics/Neutrino.h"

#include "tools/Optimization.h"
#include "tools/RNG.h"

class DecaySpace : public PhaseSpace<Decay::Channel>
{
	public:
		DecaySpace(Neutrino N = Neutrino()) : PhaseSpace(std::move(N)) { };

		/*
		std::pair<std::vector<Particle>, double> Generate(C chan, const Particle &part,
								  const Mixing &mix = Mixing());
		*/

		std::pair<std::vector<Particle>, double> Generate(Decay::Channel chan,
						const TLorentzVector &frame, const Mixing &mix = Mixing());

	private:
		bool SetDecay(Decay::Channel chan);
		double Gamma(Decay::Channel chan, const Mixing &mix = Mixing());

		double nEE(const Mixing &mix = Mixing());
		double nMM(const Mixing &mix = Mixing());
		std::pair<double, double> NeutrinoLeptonAA(Decay::Channel chan, Decay::Channel chan_o, double m_neut, double m_lepton);
		double NeutrinoLeptonAA_max(double p[]);

		double nEM(const Mixing &mix = Mixing());
		double nET(const Mixing &mix = Mixing());
		double nMT(const Mixing &mix = Mixing());
		std::pair<double, double> NeutrinoLeptonAB(Decay::Channel chan, Decay::Channel chan_o, double m_neut, double m_leptonA, double m_leptonB);
		double NeutrinoLeptonAB_max(double p[]);

		//double NeutrinoLeptonLepton_max(double x, double y, double z, double gL, double gR);
		//double F_NeutrinoLeptonLepton_max(double p[]);
		//double NeutrinoLeptonLepton(double s, double t, double cos0, double cos1, double x, double y, double z, double gL, double gR);


		double nPI0(const Mixing &mix = Mixing());
		double nETA(const Mixing &mix = Mixing());
		double nETAi(const Mixing &mix = Mixing());
		double NeutrinoPseudoMeson(Decay::Channel chan, double m_neut, double m_meson);
		double NeutrinoPseudoMeson_max(double cos0);

		double EPI(const Mixing &mix = Mixing());
		double MPI(const Mixing &mix = Mixing());
		double TPI(const Mixing &mix = Mixing());
		double EKA(const Mixing &mix = Mixing());
		double MKA(const Mixing &mix = Mixing());
		double EDs(const Mixing &mix = Mixing());
		double LeptonPseudoMeson(Decay::Channel chan, double m_lepton, double m_meson);
		double LeptonPseudoMeson_max(double cos0);

		//double ToPseudoMeson_max(double x, double y);
		//double ToPseudoMeson_cos0_max(double cos0);
		//double ToPseudoMeson(double cos0, double x, double y);



		double nRHO0(const Mixing &mix = Mixing());
		double nOMEGA(const Mixing &mix = Mixing());
		double nPHI(const Mixing &mix = Mixing());
		double NeutrinoVectorMeson(Decay::Channel chan, double m_neut, double m_meson);
		double NeutrinoVectorMeson_max(double cos0);

		double ERHO(const Mixing &mix = Mixing());
		double MRHO(const Mixing &mix = Mixing());
		double EKAx(const Mixing &mix = Mixing());
		double MKAx(const Mixing &mix = Mixing());
		double LeptonVectorMeson(Decay::Channel chan, double M_Neut, double M_Meson);
		double LeptonVectorMeson_max(double cos0);

		//double ToVectorMeson_max(double x, double y);
		//double ToVectorMeson_cos0_max(double cos0);
		//double ToVectorMeson(double cos0, double x, double y);

		void Reset();

		std::unordered_map<Decay::Channel, double> _table;
};

#endif
