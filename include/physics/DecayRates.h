/*
 * Full decay width calculator for a sterile neutrino below 2GeV in a minimal model
 * Author: Tommaso Boschi
 */

#ifndef DECAYRATES_H
#define DECAYRATES_H

#include <iostream>
#include <fstream>
#include <vector>
#include <cstring>

#include "tools.h"
#include "Amplitude.h"

class DecayRates : public Amplitude
{
	public:
		DecayRates(const Neutrino &N);

		bool IsAllowed(Channel::Name chan);
		double MassThreshold(Channel::Name chan);

		double Gamma(Channel::Name chan, const std::array<double, 3> &mix);
		double Gamma(Channel::Name chan, double ue = 0., double um = 0., double ut = 0.);
		double Other(Channel::Name chan, const std::array<double, 3> &mix);
		double Other(Channel::Name chan, double ue = 0., double um = 0., double ut = 0.);
		double Branch(Channel::Name chan, const std::array<double, 3> &mix);
		double Branch(Channel::Name chan, double ue = 0., double um = 0., double ut = 0.);

		double Total(double ue = 0., double um = 0., double ut = 0.);
		double ExpALL(double ue = 0., double um = 0., double ut = 0.);

		double nnn(double ue = 0., double um = 0., double ut = 0.);
		double nGAMMA(double ue = 0., double um = 0., double ut = 0.);
		double nEE(double ue = 0., double um = 0., double ut = 0.);
		double nEM(double ue = 0., double um = 0., double ut = 0.);
		double nMM(double ue = 0., double um = 0., double ut = 0.);
		double nET(double ue = 0., double um = 0., double ut = 0.);
		double nMT(double ue = 0., double um = 0., double ut = 0.);
		double nPI0(double ue = 0., double um = 0., double ut = 0.);
		double EPI(double ue = 0., double um = 0., double ut = 0.);
		double MPI(double ue = 0., double um = 0., double ut = 0.);
		double TPI(double ue = 0., double um = 0., double ut = 0.);
		double EKA(double ue = 0., double um = 0., double ut = 0.);
		double MKA(double ue = 0., double um = 0., double ut = 0.);
		double EKAx(double ue = 0., double um = 0., double ut = 0.);
		double MKAx(double ue = 0., double um = 0., double ut = 0.);
		double nRHO0(double ue = 0., double um = 0., double ut = 0.);
		double ERHO(double ue = 0., double um = 0., double ut = 0.);
		double MRHO(double ue = 0., double um = 0., double ut = 0.);
		double nOMEGA(double ue = 0., double um = 0., double ut = 0.);
		double nETA(double ue = 0., double um = 0., double ut = 0.);
		double nETAi(double ue = 0., double um = 0., double ut = 0.);
		double nPHI(double ue = 0., double um = 0., double ut = 0.);
		double ECHARM(double ue = 0., double um = 0., double ut = 0.);

		double LeptonPseudoMeson(double m_lepton, double m_meson);
		double I_LeptonPseudoMeson(double x, double y);

		double NeutrinoPseudoMeson(double m_lepton, double m_meson);
		double I_NeutrinoPseudoMeson(double x, double y);

		double LeptonVectorMeson(double m_lepton, double m_meson);
		double I_LeptonVectorMeson(double x, double y);

		double NeutrinoVectorMeson(double m_lepton, double m_meson);
		double I_NeutrinoVectorMeson(double x, double y);

		std::pair<double, double> NeutrinoLeptonAA(double m_neut, double m_lepton);
		std::pair<double, double> NeutrinoLeptonAB(double m_neut, double m_leptona, double m_leptonb);
		double NeutrinoLeptonLepton(double x, double y, double z, double gL, double gR);
		double I_NeutrinoLeptonLepton(double x, double y, double z, double gL, double gR);//, double theta)
		double I_NeutrinoLeptonLepton_s(double s);

		void Reset();

	private:
		std::map<Channel::Name, double> fdecay;
};

#endif
