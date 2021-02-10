/*
 * Decay rate simplifier for a three body decay of muon and kaon(0)
 * Author: Tommaso Boschi
 */

#ifndef PRODUCTIONRATE_H
#define PRODUCTIONRATE_H

// base class
#include "physics/Amplitude.h"
#include "physics/Productions.h"

#include <unordered_map>
#include "tools/Integration.h"

class ProductionRate : public Amplitude
{
	public:
		ProductionRate(Neutrino N = Neutrino());

		double Gamma(Production::Channel chan, const Mixing &mix = Mixing());
		double Scale(Production::Channel chan, const Mixing &mix = Mixing());

	private:
		double Total(const Mixing &mix = Mixing());
		double MuonE(const Mixing &mix = Mixing());
                double MuonM(const Mixing &mix = Mixing());
                double TauEE(const Mixing &mix = Mixing());
                double TauET(const Mixing &mix = Mixing());
                double TauMM(const Mixing &mix = Mixing());
                double TauMT(const Mixing &mix = Mixing());
		double TauPI(const Mixing &mix = Mixing());
		double Tau2PI(const Mixing &mix = Mixing());
                double PionE(const Mixing &mix = Mixing());
                double PionM(const Mixing &mix = Mixing());
                double KaonE(const Mixing &mix = Mixing());
                double KaonM(const Mixing &mix = Mixing());
                double CharmE(const Mixing &mix = Mixing());
                double CharmM(const Mixing &mix = Mixing());
                double CharmT(const Mixing &mix = Mixing());
                double Kaon0E(const Mixing &mix = Mixing());
                double Kaon0M(const Mixing &mix = Mixing());
		double KaonCE(const Mixing &mix = Mixing());
                double KaonCM(const Mixing &mix = Mixing());

		double AntileptonNeutrinoDecay(double m_leptonA, double m_leptonB, double M_Neut);
		//double I_AntileptonNeutrino(double x, double y, double z);
		double AntileptonNeutrino_s(double s);

		double LeptonNeutrinoDecay(double m_leptonA, double m_leptonB, double M_Neutrino);
		//double I_LeptonNeutrino(double x, double y, double z);
		double LeptonNeutrino_u(double u);

		double LeptonTwoDecay(double m_lepton, double m_meson);
		//double I_LeptonTwo(double x, double y);

		double LeptonThreeDecay(double m_lepton, double m_meson, double m_meson0);
		//double I_LeptonThree(double x, double y, double z);
		double LeptonThree_s(double s);

		double MesonTwoDecay(double m_meson, double m_lepton);
		//double I_MesonTwo(double x, double y);

		double MesonThreeDecay(double m_meson0, double m_meson, double m_lepton, double L_, double L0);
		//double I_MesonThree(double x, double y, double z, double L_, double L0);
		double MesonThree_s_t(double s, double t);

		void Reset();

	private:
		static ProductionRate _sm;
		bool kSM;

		std::unordered_map<Production::Channel, double> _table;
};

#endif
