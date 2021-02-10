/*
 * Full decay width calculator for a sterile neutrino below 2GeV in a minimal model
 * Author: Tommaso Boschi
 */

#ifndef DECAYRATE_H
#define DECAYRATE_H

// base class
#include "physics/Amplitude.h"
#include "physics/Decays.h"

#include <unordered_map>
#include "tools/Integration.h"

class DecayRate : public Amplitude
{
	public:
		DecayRate(Neutrino N = Neutrino()) : Amplitude(std::move(N)) { };

		double Gamma(Decay::Channel chan, const Mixing &mix = Mixing());
		double Other(Decay::Channel chan, const Mixing &mix = Mixing());
		double Branch(Decay::Channel chan, const Mixing &mix = Mixing());

		double Total(const Mixing &mix = Mixing());
		double ExpALL(const Mixing &mix = Mixing());

	private:
		double nnn(const Mixing &mix = Mixing());
		double nGAMMA(const Mixing &mix = Mixing());
		double nEE(const Mixing &mix = Mixing());
		double nEM(const Mixing &mix = Mixing());
		double nMM(const Mixing &mix = Mixing());
		double nET(const Mixing &mix = Mixing());
		double nMT(const Mixing &mix = Mixing());
		double nPI0(const Mixing &mix = Mixing());
		double EPI(const Mixing &mix = Mixing());
		double MPI(const Mixing &mix = Mixing());
		double TPI(const Mixing &mix = Mixing());
		double EKA(const Mixing &mix = Mixing());
		double MKA(const Mixing &mix = Mixing());
		double EKAx(const Mixing &mix = Mixing());
		double MKAx(const Mixing &mix = Mixing());
		double nRHO0(const Mixing &mix = Mixing());
		double ERHO(const Mixing &mix = Mixing());
		double MRHO(const Mixing &mix = Mixing());
		double nOMEGA(const Mixing &mix = Mixing());
		double nETA(const Mixing &mix = Mixing());
		double nETAi(const Mixing &mix = Mixing());
		double nPHI(const Mixing &mix = Mixing());
		double EDs(const Mixing &mix = Mixing());

		double LeptonPseudoMeson(double m_lepton, double m_meson);
		//double I_LeptonPseudoMeson(double x, double y);

		double NeutrinoPseudoMeson(double m_lepton, double m_meson);
		//double I_NeutrinoPseudoMeson(double x, double y);

		double LeptonVectorMeson(double m_lepton, double m_meson);
		//double I_LeptonVectorMeson(double x, double y);

		double NeutrinoVectorMeson(double m_lepton, double m_meson);
		//double I_NeutrinoVectorMeson(double x, double y);

		std::pair<double, double> NeutrinoLeptonAA(double m_neut, double m_lepton);
		std::pair<double, double> NeutrinoLeptonAB(double m_neut, double m_leptonA, double m_leptonB);
		double NeutrinoLeptonLepton(double x, double y, double z, double gL, double gR);
		//double I_NeutrinoLeptonLepton(double x, double y, double z, double gL, double gR);//, double theta)
		double NeutrinoLeptonLepton_s(double s);

		void Reset();

	private:
		std::unordered_map<Decay::Channel, double> _table;
};

#endif
