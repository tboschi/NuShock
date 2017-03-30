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

//GENIE include
//#include "GHepParticle.h"
//#include "Constants.h"

#include "Tools.h"

enum ChannelName
{
	_undefined,
	_ALL,
	_nnn,
	_nGAMMA,
	_nEE,
	_nEMU,
	_nPI0,
	_EPI,
	_MUPI,
	_nMUMU,
	_EKA,
	_nKA0
};


class Decay
{
	public:
		Decay(double MSterile, double Ue, double Um, double Ut);	//Decay rates calculator

		void MapInit();
		//Decay width with A, B, and K the enhancement factors
		double Gamma(std::string Channel, double B = 1.0);
		double Other(std::string Channel, double A = 1.0);
		double Branch(std::string Channel, double A = 1.0, double B = 1.0);
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
		double GetMSterile();
		double GetUe();
		double GetUm();
		double GetUt();

		void SetMSterile(double X);
		void SetUe(double X);
		void SetUm(double X);
		void SetUt(double X);

	private:
		double M_Sterile;
		double U_e, U_m, U_t;

		//Masses
		double M_Neutrino;
		double M_Electron;
		double M_Muon;
		double M_Pion;
		double M_Pion0;
		double M_Kaon;
		double M_Kaon0;

		//Maps
		std::map<std::string, ChannelName> mapChannel;
		std::map<std::string, double> mapEnhance;
};

#endif
