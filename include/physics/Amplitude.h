/*
 * Decay rate calculator for a sterile neutrino below 2GeV in a minimal model
 * Author: Tommaso Boschi
 */

#ifndef AMPLITUDE_H
#define AMPLITUDE_H

#include <iostream>
#include <fstream>
#include <vector>
#include <cstring>

#include "tools.h"

class Amplitude
{
	public:
		Amplitude();

		void LoadMass(Channel Name);

		double Kallen(double X, double Y, double Z);
		double SqrtKallen(double X, double Y, double Z);

		double dGammad5_3B(double M2);
		double dGammad2_3B(double M2);
		double dGammad2_2B(double M2, double x, double y);
		double dGammad0_2B(double M2, double x, double y);

		double Limit(double &s, double x, double y, double z);
		double Limit(double &s, double &t, double x, double y, double z);

		double M2_LeptonPseudoMeson(int hel, double cos0, double x, double y);
		double M2_LeptonPseudoMeson(double cos0, double x, double y);
		double M2_NeutrinoPseudoMeson(int hel, double cos0, double x, double y);
		double M2_NeutrinoPseudoMeson(double cos0, double x, double y);
		double M2_LeptonVectorMeson(int hel, double cos0, double x, double y);
		double M2_LeptonVectorMeson(double cos0, double x, double y);
		double M2_NeutrinoVectorMeson(int hel, double cos0, double x, double y);
		double M2_NeutrinoVectorMeson(double cos0, double x, double y);
		double M2_WW(int hel, double s, double cos0, double x, double y, double z);
		double M2_WW(double s, double cos0, double x, double y, double z);
		double M2_WZ(int hel, double s, double t, double cos0s, double cos0t, double x, double y, double z);
		double M2_WZ(double s, double t, double cos0s, double cos0t, double x, double y, double z);
		double M2_WZ(int hel, double u, double cos0u, double x, double y, double z);
		double M2_WZ(double u, double cos0u, double x, double y, double z);

		double M2_LeptonNeutrino(int hel, double s, double x, double y, double z);
		double M2_LeptonNeutrino(double s, double x, double y, double z);
		double M2_AntiLeptonNeutrino(int hel, double s, double x, double y, double z);
		double M2_AntiLeptonNeutrino(double s, double x, double y, double z);
		double M2_LeptonTwo(int hel, double x, double y);
		double M2_LeptonTwo(double x, double y);
		double M2_LeptonThree(int hel, double x, double y, double z);
		double M2_LeptonThree(double x, double y, double z);
		double M2_MesonTwo(int hel, double x, double y);
		double M2_MesonTwo(double x, double y);
		double M2_MesonThree(int hel, double s, double t, double x, double y, double z, double L_, double L0);
		double M2_MesonThree(double s, double t, double x, double y, double z, double L_, double L0);

		virtual void Reset() { ; }

		void SetNeutrino(const double &mass,
				 const std::array<double, 3> &mixings,
				 const size_t &opts);

	protected:
		template <typedef F> // using Boole's method
		double Integration(F func, size_t steps = 1000)
		{
			double h = 1./steps;
			double res = 0.;

			double left = 7. * this->func(0);
			for (size_t s = 0; s < steps; ++s) {
				res += left;

				res += 32. * this->func((s + .25) * h);
				res += 12. * this->func((s + 0.5) * h);
				res += 32. * this->func((s + .75) * h);

				left = 7. * this->func((s + 1.) * h);

				res += left;
			}	

			return res * h / 90.;
		}	

		Channel _channel;
		double _m_parent, _m_Nu;
		double _Ue2, _Um2, _Ut2;
		size_t _opts;

		std::vector<std::pair<double, int> > _masspdg;
		std::vector<double> _vars;

	private:
};

#endif
