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
#include "TH1.h"
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "TGenPhaseSpace.h"

//GENIE include
//#include "GHepParticle.h"
//#include "Constants.h"

#include "Tools.h"

enum ChannelName
{
	_undefined,
	_ALL,		//
	_nnn,		//3 body
	_nGAMMA,	//2 body
	_nEE,		//3 body
	_nEMU,		//3 body
	_nPI0,		//3 body
	_EPI,		//2 body
	_nMUMU,		//3 body
	_MUPI,		//2 body
	_EKA,		//2 body
	_nKA0		//3 body
};


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

		double Total();
		double nnn();
		double nGAMMA();
		double nEE();
		double nEMU();
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
};

#endif
