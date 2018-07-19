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

#include "Physics/Amplitude.h"
#include "cuba.h"

class Production : public Amplitude
{
	public:
		Production();

		Amplitude::Channel FindChannel(std::string Name);
		std::string FindChannel(Amplitude::Channel Name);
		std::vector<Amplitude::Channel> ListChannels();

		bool IsAllowed(Channel Name);
		double Gamma(Channel Name, bool Unitary = false);
		double Scale(Channel Name);

		double Total();
		double MuonE();
                double MuonM();
                double TauEE();
                double TauET();
                double TauMM();
                double TauMT();
		double TauPI();
		double Tau2PI();
                double PionE();
                double PionM();
                double KaonE();
                double KaonM();
                double CharmE();
                double CharmM();
                double CharmT();
                double Kaon0E();
                double Kaon0M();
		double KaonCE();
                double KaonCM();

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
		double I_MesonThree_D(double *x);

		void Reset();
		void SetFunction(double (Production::*FF)(double));
		void SetFunction_D(double (Production::*FF)(double*));

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
