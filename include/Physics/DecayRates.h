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

#include "Tools.h"
#include "Physics/Amplitude.h"

class DecayRates : public Amplitude
{
	public:
		DecayRates();

		Amplitude::Channel FindChannel(std::string Name);
		std::vector<std::string> ListChannel();

		bool IsAllowed(Channel Name);
		double Gamma(Channel Name);
		double Other(Channel Name);
		double Branch(Channel Name);

		double Total();
		double ExpALL();

		double nnn();
		double nGAMMA();
		double nEE();
		double nEM();
		double nME();
		double nMM();
		double nET();
		double nTE();
		double nMT();
		double nTM();
		double nPI0();
		double EPI();
		double MPI();
		double TPI();
		double EKA();
		double MKA();
		double EKAx();
		double MKAx();
		double nRHO0();
		double ERHO();
		double MRHO();
		double nOMEGA();
		double nETA();
		double nETAi();
		double nPHI();
		double ECHARM();

		double LeptonPseudoMeson(double M_Lepton, double M_Meson);
		double NeutrinoPseudoMeson(double M_Meson, double fDecay2);
		double LeptonVectorMeson(double M_Lepton, double M_Meson);
		double NeutrinoVectorMeson(double M_Meson, double fDecay2);
		double I_LeptonPseudoMeson(double x, double y);
		double I_NeutrinoPseudoMeson(double x, double y);
		double I_LeptonVectorMeson(double x, double y);
		double I_NeutrinoVectorMeson(double x, double y);

		double NeutrinoLeptonAA(double &fCC, double &fNC, double M_Neut, double M_Lepton);
		double NeutrinoLeptonAB(double M_Neut, double M_LeptonA, double M_LeptonB);
		double NeutrinoLeptonLepton(double x, double y, double z, double gL, double gR);
		double I_NeutrinoLeptonLepton(double x, double y, double z, double gL, double gR);//, double theta)
		double I_NeutrinoLeptonLepton_s(double s);

		//Set and Get
		//std::vector<std::string> ListChannels();
		//std::vector<ChannelName> ListNameKeys();

		void Reset();
		void SetFunction(double (DecayRates::*FF)(double));

	private:
		double fnnn,
                       fnGAMMA,
                       fnEE_e,	
                       fnEE_o,	
                       fnEM,	
                       fnME,	
                       fnMM_m,	
                       fnMM_o,	
                       fnET,	
                       fnTE,	
                       fnMT,	
                       fnTM,	
                       fnPI0,	
                       fEPI,	
                       fMPI,	
                       fTPI,	
                       fEKA,	
                       fMKA,	
                       fnRHO0,	
                       fERHO,	
                       fMRHO,	
                       fEKAx,	
                       fMKAx,	
                       fnETA,	
                       fnETAi,	
                       fnOMEGA,
                       fnPHI,	
                       fECHARM;
};

#endif
