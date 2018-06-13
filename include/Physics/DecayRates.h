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

class DecayRates : protected Amplitude
{
	public:
		DecayRates();

		bool IsAllowed(Channel Name)
		double Gamma(Channel Name);
		double Other(Channel Name);
		double Branch(Channel Name);

		double Total();
		double ExpALL();

		double nnn();
		double nGAMMA();
		double nEE();
		double nEMU();
		double nMUE();
		double nMUMU();
		double nET();
		double nTE();
		double nMUT();
		double nTMU();
		double nPI0();
		double EPI();
		double MUPI();
		double TPI();
		double EKA();
		double MUKA();
		double EKAx();
		double MUKAx();
		double nRHO0();
		double ERHO();
		double MURHO();
		double nOMEGA();
		double nETA();
		double nETAi();
		double nPHI();
		double ECHARM();

		double LeptonPseudoMeson(double M_Lepton, double M_Meson);
		double NeutrinoPseudoMeson(double M_Meson, double fDecay2);
		double I_PseudoMeson(double x, double y);
		double LeptonVectorMeson(double M_Lepton, double M_Meson);
		double NeutrinoVectorMeson(double M_Meson, double fDecay2);

		double NeutrinoLeptonAA(double &fCC, double &fNC, double M_Lepton);
		double NeutrinoLeptonAB(double M_LeptonA, double M_LeptonB);
		double I_WW(double x, double y, double z, double theta);
		double I_WW_s(double s);
		double I_WW_s(double s, double cos0, double x, double y, double z);
		double I_WZ(double x, double y, double z, double theta);
		double I_WZ_s(double s);
		double I_WZ_s(double s, double cos0, double x, double y, double z);

		//Set and Get
		//std::vector<std::string> ListChannels();
		//std::vector<ChannelName> ListNameKeys();

		void Reset();

	private:
		double fnnn,
                       fnGAMMA,
                       fnEE_e,	
                       fnEE_mt,	
                       fnEMU,	
                       fnMUE,	
                       fnMUMU_m,	
                       fnMUMU_et,	
                       fnET,	
                       fnTE,	
                       fnMUT,	
                       fnTMU,	
                       fnPI0,	
                       fEPI,	
                       fMUPI,	
                       fTPI,	
                       fEKA,	
                       fMUKA,	
                       fnRHO0,	
                       fERHO,	
                       fMURHO,	
                       fEKAx,	
                       fMUKAx,	
                       fnETA,	
                       fnETAi,	
                       fnOMEGA,
                       fnPHI,	
                       fECHARM;
};

#endif
