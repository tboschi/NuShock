/*
 * Decay rate calculator for a sterile neutrino below 500MeV in a minimal model
 * Other decay channels could be added
	 * Usage: create a generator Decay(f,f,f,f) defining sterile mass and the mixing elements.
	 *
 * Author: Tommaso Boschi
 */

#ifndef DECAYRATES_H
#define DECAYRATES_H

#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <cstring>

//ROOT include
#include "TLorentzVector.h"
#include "TGenPhaseSpace.h"

#include "Tools.h"
#include "ThreeBody.h"


class Decay
{
	public:
		Decay(double MSterile = 0.0, double Ue = 0.0, double Um = 0.0, double Ut = 0.0);	//Decay rates calculator

		void MapInit();
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
		double nKA0x();
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

		void SetNvec(TLorentzVector &X);
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
                       fnKA0x,	
                       fMUKAx,	
                       fnETA,	
                       fnETAi,	
                       fnOMEGA,
                       fnPHI,	
                       fECHARM;

		//Masses
		const double M_Neutrino;
		const double M_Photon;
		const double M_Electron;
		const double M_Muon;
		const double M_Pion;
		const double M_Pion0;
		const double M_Kaon;
		const double M_Kaon0;
		const double M_Eta;
		const double M_Rho;
		const double M_Rho0;
		const double M_Omega;
		const double M_Kaonx;
		const double M_Kaon0x;
		const double M_Etai;
		const double M_Phi;
		const double M_Tau;
		const double M_Charm;

		//Maps
		std::map<std::string, ChannelName> mapChannel;
		std::map<std::string, double> mapEnhance;

		//Generate PhaseSpace
		TGenPhaseSpace *Event;
		TLorentzVector *N_vec, *N_rest;
		ThreeBody *TheSpace;
};

#endif
