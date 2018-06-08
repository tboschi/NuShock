/*
 * Decay rate simplifier for a three body decay of muon and kaon(0)
 * Author: Tommaso Boschi
 */

#ifndef THREEBODY_H
#define THREEBODY_H

#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <cstring>

#include "Tools/Tools.h"


class ProductionRates : protected DecayRates
{
	public:
		ProductionRates();

		bool IsAllowed(Channel Name)
		double dGamma(Channel Name);
		double Scale(Channel Name);

		double Total();
		double MuonE();
		double MuonE();
		double MuonE();
                double MuonM();
                double TauEE();
                double TauET();
                double TauMM();
                double TauMT();
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
		double d_LeptonNeutrino(double x, double y, double z, double theta);
		double d_LeptonNeutrino_s(double s);
		double d_LeptonNeutrino_t(double t);
		double d_LeptonNeutrino(double X, double Y, double x, double y, double z, double cos0);
		double LeptonAntineutrinoDecay(double M_LeptonA, double M_LeptonB, double M_Neutrino);
		double d_LeptonAntineutrino(double x, double y, double z, double theta);
		double d_LeptonAntineutrino_s(double s);
		double d_LeptonAntineutrino_s(double S, double x, double y, double z, double cos0);
		double LeptonMesonDecay(double M_Lepton, double M_Meson, double fDecay2);
		double d_LeptonMeson(double x, double y);
		double MesonTwoDecay(double M_Meson, double M_Lepton, double fDecay2);
		double d_MesonTwo(double x, double y);
		double MesonThreeDecay(double M_Meson0, double M_Meson1, double M_Lepton, double vCKM, double L_, double L0);
		double d_MesonThree(double x, double y, double z, double L_, double L0);
		double d_MesonThree_s(double s);
		double d_MesonThree_t(double t);
		double d_MesonThree(double Z, double Y, double x, double y, double cos0, double L_, double L0);

		bool IsChanged()

	private:

		double fMuonE,
                       fMuonM,
                       fTauEE,
                       fTauET,
                       fTauME,
                       fTauMT,
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
