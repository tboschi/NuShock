/*
 * Full decay width calculator for a sterile neutrino below 2GeV in a minimal model
 * Author: Tommaso Boschi
 */

#ifndef FULLWIDTH_H
#define FULLWIDTH_H

#include <iostream>
#include <fstream>
#include <vector>
#include <cstring>

//ROOT include
#include "TLorentzVector.h"
#include "TGenPhaseSpace.h"

#include "Tools/Tools.h"
#include "Physics/DecayRates.h"


class FullWidth : protected DecayRates
{
	public:
		FullWidth(double MSterile = 0.0, double Ue = 0.0, double Um = 0.0, double Ut = 0.0);

		//Decay width with A, B, and K the enhancement factors
		double Gamma(std::string Channel, double B = 1.0);
		double Other(std::string Channel, double A = 1.0);
		double Branch(std::string Channel, double A = 1.0, double B = 1.0);
		int PhaseSpace(std::string Channel, double &Weight);
		bool IsAllowed(std::string Channel);

		void SetEnhancement(std::string Channel = "ALL", double K = 1.0);

		//M2
		double yLim(double &Min, double &Max, double Ex);	//y integration limits
		double xLim(double &Min, double &Max);
		double M2_W(double x, double y, double a, double b, double c);
		double M2_Z(double x, double y, double a, double b, double c);
		double M2_WZ(double x, double y, double a, double b, double c);
		double M2_nLL(double Ex, double Ey, double M_Lepton);
		double M2_nEE(double Ex, double Ey);
		double M2_nMUMU(double Ex, double Ey);
		double M2_nEMU(double Ex, double Ey);
		double ddGamma(double (Decay::*M2)(double, double), double Ex, double Ey);
		double ddGamma(std::string Channel, double Ex, double Ey);

		double Total(double A = 1.0);
		double ExpALL(double A = 1.0);

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

		double LeptonPseudoMeson(double M_Lepton, double M_Meson, double vCKM, double fDecay2);
		double NeutrinoPseudoMeson(double M_Meson, double fDecay2);
		double LeptonVectorMeson(double M_Lepton, double M_Meson, double vCKM, double fDecay2, double fVector);
		double NeutrinoVectorMeson(double M_Meson, double fDecay2, double fVector);
		double NeutrinoLeptonAA(double &fCC, double &fNC, double M_Lepton);
		double NeutrinoLeptonAB(double M_LeptonA, double M_LeptonB);

		//Set and Get
		std::vector<std::string> ListChannels();
		std::vector<ChannelName> ListNameKeys();

		TLorentzVector *GetNvec();
		TLorentzVector GetDecayProduct(int i, int &ID);

		bool IsChanged();

		double GetMass();
		double GetUe();
		double GetUm();
		double GetUt();

		void SetMass(double X);
		void SetUe(double X);
		void SetUm(double X);
		void SetUt(double X);

	private:
		int PdgCode[3];

		double M_Sterile, M_Sterile_prev;
		double Ue, Um, Ut;

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
