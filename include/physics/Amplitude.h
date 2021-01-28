/*
 * Decay rate calculator for a sterile neutrino below 2GeV in a minimal model
 * Author: Tommaso Boschi
 */

#ifndef AMPLITUDE_H
#define AMPLITUDE_H

#include <iostream>
#include <fstream>
#include <vector>
#include <map>

#include "physics/Const.h"
#include "physics/Mixings.h"
#include "physics/Channels.h"
#include "physics/Neutrino.h"

#include <numeric>

class Amplitude
{
	public:
		Amplitude(Neutrino N = Neutrino());

		double MassThreshold(Channel::Name chan);
		void SetNeutrino(const Neutrino &N);

	protected:
		void LoadMass(Channel::Name chan);

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
		virtual void Reset() { _table.clear(); }

		double _m_parent;
		Channel::Name _channel;
		Neutrino _N;

		std::vector<std::pair<double, int> > _masspdg;
		std::vector<double> _vars;

		std::map<Channel::Name, double> _table;
};

#endif
