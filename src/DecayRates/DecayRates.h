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
#include "3Body.h"


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
		double nnn();
		double nGAMMA();
		double nEE();
		double nEMU();
		double nMUE();
		double nLeptonW(double m1, double m2);
		double nPI0();
		double EPI();
		double MUPI();
		double nMUMU();
		double EKA();
		double nKA0();

		//Set and Get
		std::vector<std::string> ListChannels();

		TLorentzVector *GetNvec();
		TLorentzVector *GetDecayProduct(int i);

		bool IsChanged(std::string Name);

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
		double M_Sterile;
		double U_e, U_m, U_t;
		double M_Sterile_prev;
		double U_e_prev, U_m_prev, U_t_prev;

		double fTotal;	//Total Gamma, to avoid recalculation
		double fnnn, fnGAMMA, fnEE, fnEMU, fnMUE, fnPI0, fEPI, fnMUMU, fMUPI, fEKA, fnKA0;

		//Masses
		const double M_Neutrino;
		const double M_Photon;
		const double M_Electron;
		const double M_Muon;
		const double M_Pion;
		const double M_Pion0;
		const double M_Kaon;
		const double M_Kaon0;

		//Maps
		std::map<std::string, ChannelName> mapChannel;
		std::map<std::string, double> mapEnhance;

		//Generate PhaseSpace
		TGenPhaseSpace *Event;
		TLorentzVector *N_vec;
		ThreeBody *TheSpace;
};

#endif
