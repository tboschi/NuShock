/*
 * Decay rate simplifier for a three body decay of muon and kaon(0)
 * Author: Tommaso Boschi
 */

#ifndef PRODUCTION_H
#define PRODUCTION_H

#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <cstring>

#include "Amplitude.h"

class Production : public Amplitude
{
	public:
		Production();

		bool IsAllowed(const Channel::Name &chan);
		double MassThreshold(const Channel::Name &chan);
		double Gamma(Channel::Name chan, const std::array<double, 3> &mix);
		double Gamma(Channel::Name chan, double ue = 0., double um = 0., double ut = 0.);
		double Scale(const Channel::Name &chan);

		double Total(double ue = 0., double um = 0., double ut = 0.);
		double MuonE(double ue = 0., double um = 0., double ut = 0.);
                double MuonM(double ue = 0., double um = 0., double ut = 0.);
                double TauEE(double ue = 0., double um = 0., double ut = 0.);
                double TauET(double ue = 0., double um = 0., double ut = 0.);
                double TauMM(double ue = 0., double um = 0., double ut = 0.);
                double TauMT(double ue = 0., double um = 0., double ut = 0.);
		double TauPI(double ue = 0., double um = 0., double ut = 0.);
		double Tau2PI(double ue = 0., double um = 0., double ut = 0.);
                double PionE(double ue = 0., double um = 0., double ut = 0.);
                double PionM(double ue = 0., double um = 0., double ut = 0.);
                double KaonE(double ue = 0., double um = 0., double ut = 0.);
                double KaonM(double ue = 0., double um = 0., double ut = 0.);
                double CharmE(double ue = 0., double um = 0., double ut = 0.);
                double CharmM(double ue = 0., double um = 0., double ut = 0.);
                double CharmT(double ue = 0., double um = 0., double ut = 0.);
                double Kaon0E(double ue = 0., double um = 0., double ut = 0.);
                double Kaon0M(double ue = 0., double um = 0., double ut = 0.);
		double KaonCE(double ue = 0., double um = 0., double ut = 0.);
                double KaonCM(double ue = 0., double um = 0., double ut = 0.);

		double AntiLeptonNeutrinoDecay(double M_LeptonA, double M_LeptonB, double M_Neut);
		double I_AntiLeptonNeutrino(double x, double y, double z);
		double I_AntiLeptonNeutrino_s(double s);

		double LeptonNeutrinoDecay(double M_LeptonA, double M_LeptonB, double M_Neutrino);
		double I_LeptonNeutrino(double x, double y, double z);
		double I_LeptonNeutrino_u(double u);

		double LeptonTwoDecay(double M_Lepton, double M_Meson);
		double I_LeptonTwo(double x, double y);

		double LeptonThreeDecay(double M_Lepton, double M_Meson, double M_Meson0);
		double I_LeptonThree(double x, double y, double z);
		double I_LeptonThree_s(double s);

		double MesonTwoDecay(double M_Meson, double M_Lepton);
		double I_MesonTwo(double x, double y);

		double MesonThreeDecay(double M_Meson0, double M_Meson1, double M_Lepton, double L_, double L0);
		double I_MesonThree(double x, double y, double z, double L_, double L0);
		double I_MesonThree_s(double s);
		double I_MesonThree_t(double t);

		void Reset();

	private:

		double fMuonE,
                       fMuonM,
                       fTauEE,
                       fTauET,
                       fTauMM,
                       fTauMT,
		       fTauPI,
		       fTau2PI,
                       fPionE,
                       fPionM,
                       fKaonE,
                       fKaonM,
                       fCharmE,
                       fCharmM,
                       fCharmT,
                       fKaon0E,
                       fKaon0M,
                       fKaonCE,
                       fKaonCM;
};

#endif
