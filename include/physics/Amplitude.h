/*
 * Decay rate calculator for a sterile neutrino below 2GeV in a minimal model
 * Author: Tommaso Boschi
 */

#ifndef AMPLITUDE_H
#define AMPLITUDE_H

#include <iostream>
#include <fstream>
#include <vector>

#include "physics/Const.h"
#include "physics/Mixings.h"
#include "physics/Decays.h"
#include "physics/Productions.h"
#include "physics/Neutrino.h"

#include <numeric>

class Amplitude
{
	public:
		Amplitude(Neutrino N = Neutrino());
		void SetNeutrino(const Neutrino &N);
		Neutrino GetNeutrino() const;

		bool IsAllowed(Decay::Channel chan);
		static bool IsAllowed(const Neutrino &N, Decay::Channel chan) {
			return IsAllowed(N.M(), chan);
		}
		static bool IsAllowed(double mass, Decay::Channel chan) {
			return (mass >= MassThreshold(chan));
		}

		static double MassThreshold(Decay::Channel chan) {
			auto mass = Decay::Masses(chan);
			return std::accumulate(mass.begin(), mass.end(), 0.,
					[](double sum, double m) { return sum + m; });
		}

		bool IsAllowed(Production::Channel chan);
		static bool IsAllowed(const Neutrino &N, Production::Channel chan) {
			return IsAllowed(N.M(), chan);
		}
		static bool IsAllowed(double mass, Production::Channel chan) {
			return (mass <= MassThreshold(chan));
		}

		static double MassThreshold(Production::Channel chan) {
			auto mass = Production::Masses(chan);
			return std::accumulate(mass.begin()+1, mass.end(), mass.front(),
					[](double sum, double m) { return sum - m; });
		}

	protected:
		double Kallen(double X, double Y, double Z);

		double dGammad5_3B();
		double dGammad2_3B();
		double dGammad2_2B(double x, double y);
		double dGammad0_2B(double x, double y);

		double Limit(double &s, double x, double y, double z);
		double Limit(double &s, double &t, double x, double y, double z);

		double M2_LeptonPseudoMeson(double cos0, double x, double y);
		double M2_LeptonPseudoMeson(int hel, double cos0, double x, double y);
		double M2_NeutrinoPseudoMeson(double cos0, double x, double y);
		double M2_NeutrinoPseudoMeson(int hel, double cos0, double x, double y);
		double M2_LeptonVectorMeson(double cos0, double x, double y);
		double M2_LeptonVectorMeson(int hel, double cos0, double x, double y);
		double M2_NeutrinoVectorMeson(double cos0, double x, double y);
		double M2_NeutrinoVectorMeson(int hel, double cos0, double x, double y);
		double M2_WW(double s, double cos0, double x, double y, double z);
		double M2_WW(int hel, double s, double cos0, double x, double y, double z);
		double M2_WZ(double s, double t, double cos0s, double cos0t, double x, double y, double z);
		double M2_WZ(int hel, double s, double t, double cos0s, double cos0t, double x, double y, double z);
		double M2_WZ(double u, double cos0u, double x, double y, double z);
		double M2_WZ(int hel, double u, double cos0u, double x, double y, double z);

		double M2_LeptonNeutrino(double s, double x, double y, double z);
		double M2_LeptonNeutrino(int hel, double s, double x, double y, double z);
		double M2_AntileptonNeutrino(double s, double x, double y, double z);
		double M2_AntileptonNeutrino(int hel, double s, double x, double y, double z);
		double M2_LeptonTwo(double x, double y);
		double M2_LeptonTwo(int hel, double x, double y);
		double M2_LeptonThree(double x, double y, double z);
		double M2_LeptonThree(int hel, double x, double y, double z);
		double M2_MesonTwo(double x, double y);
		double M2_MesonTwo(int hel, double x, double y);
	public:
		double M2_MesonThree(double s, double t, double x, double y, double z, double L_, double L0);
		double M2_MesonThree(int hel, double s, double t, double x, double y, double z, double L_, double L0);

	protected:
		virtual void Reset() = 0;

		double _m_parent;
		Neutrino _N;

		std::vector<double> _vars;
};

#endif
