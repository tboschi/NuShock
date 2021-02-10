/*
 * PhaseSpace generator for simulation of sterile neutrino decays (it is not needed for production)
 * It also returns correct 
 *
 * Author: Tommaso Boschi
 */

#ifndef PRODUCTIONSPACE_H
#define PRODUCTIONSPACE_H

// base class
#include "physics/Amplitude.h"
#include "physics/PhaseSpace.h"

//ROOT include
#include "TGenPhaseSpace.h"

#include <random>
#include "physics/Productions.h"
#include "physics/Particle.h"
#include "physics/Neutrino.h"

#include "tools/Optimization.h"
#include "tools/RNG.h"

class ProductionSpace :  public PhaseSpace<Production::Channel>
{
	public:
		ProductionSpace(Neutrino N = Neutrino()) : PhaseSpace(std::move(N)) { };

		/*
		std::pair<std::vector<Particle>, double> Generate(C chan, const Particle &part,
								  const Mixing &mix = Mixing());
		*/

		std::pair<std::vector<Particle>, double> Generate(Production::Channel chan,
						const TLorentzVector &frame, const Mixing &mix = Mixing());

	private:
		bool SetProduction(Production::Channel chan);
		double Gamma(Production::Channel chan, const Mixing &mix = Mixing());

		double MuonE(const Mixing &mix = Mixing());
		double MuonM(const Mixing &mix = Mixing());
		double TauEE(const Mixing &mix = Mixing());
		double TauET(const Mixing &mix = Mixing());
		double TauMM(const Mixing &mix = Mixing());
		double TauMT(const Mixing &mix = Mixing());

		double LeptonNeutrino(Production::Channel chan, double m_lepton0, double m_lepton, double m_neut = 0.);
		//double LeptonNeutrino_max(double x, double y, double z);
		double LeptonNeutrino_max(double u);

		double AntileptonNeutrino(Production::Channel chan, double m_lepton0, double m_lepton, double m_neut = 0.);
		//double AntileptonNeutrino_max(double x, double y, double z);
		double AntileptonNeutrino_max(double s);

		//double TauPI(const Mixing &mix = Mixing());
		//double Tau2PI(const Mixing &mix = Mixing());

		// ratio could be 1
		//double LeptonMeson(Production::Channel chan, double m_lepton, double m_meson);
		//double LeptonMeson_max(double m_lepton, double m_meson);


		//double PionE(const Mixing &mix = Mixing());
		//double PionM(const Mixing &mix = Mixing());
		//double KaonE(const Mixing &mix = Mixing());
		//double KaonM(const Mixing &mix = Mixing());
		//double CharmE(const Mixing &mix = Mixing());
		//double CharmM(const Mixing &mix = Mixing());
		//double CharmT(const Mixing &mix = Mixing());

		// ratio could be 1
		//double MesonTwo(Production::Channel chan, double m_meson, double m_lepton);
		//double MesonTwo_max(double m_meson, double m_lepton);

		double Kaon0E(const Mixing &mix = Mixing());
		double Kaon0M(const Mixing &mix = Mixing());
		double KaonCE(const Mixing &mix = Mixing());
		double KaonCM(const Mixing &mix = Mixing());

		double MesonThree(Production::Channel chan, double m_meson0, double m_meson, double m_lepton,
				  double L_, double L0);
		//double MesonThree_max(double x, double y, double z, double L_, double L0);
		double MesonThree_max(double p[]);

		void Reset();

		std::unordered_map<Production::Channel, double> _table;
};

#endif
