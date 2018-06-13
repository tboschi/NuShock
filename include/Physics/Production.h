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

#include "Tools.h"
#include "Physics/Amplitude.h"

class Production : public Amplitude
{
	public:
		Production();

		bool IsAllowed(Channel Name);
		double Gamma(Channel Name);
		double Scale(Channel Name);

		double Total();
		double MuonE();
                double MuonM();
                double TauEE();
                double TauET();
                double TauMM();
                double TauMT();
		double TauPion();
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

		double LeptonNeutrinoDecay(double M_LeptonA, double M_LeptonB, double M_Neutrino);
		double I_LeptonNeutrino(double x, double y, double z);
		double I_LeptonNeutrino_s(double s);
		double I_LeptonNeutrino_t(double t);
		double I_LeptonNeutrino(double X, double Y, double x, double y, double z);

		double LeptonAntineutrinoDecay(double M_LeptonA, double M_LeptonB, double M_Neut);
		double I_LeptonAntineutrino(double x, double y, double z);
		double I_LeptonAntineutrino_s(double s);
		double I_LeptonAntineutrino_t(double t);

		double LeptonMesonDecay(double M_Lepton, double M_Meson);
		double I_LeptonMeson(double x, double y);

		double MesonTwoDecay(double M_Meson, double M_Lepton);
		double I_MesonTwo(double x, double y);

		double MesonThreeDecay(double M_Meson0, double M_Meson1, double M_Lepton, double L_, double L0);
		double I_MesonThree(double x, double y, double z, double L_, double L0);
		double I_MesonThree_s(double s);
		double I_MesonThree_t(double t);

		void Reset();
		void SetFunction(double (Production::*FF)(double));

	private:

		double fMuonE,
                       fMuonM,
                       fTauEE,
                       fTauET,
                       fTauMM,
                       fTauMT,
		       fTauPion,
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
