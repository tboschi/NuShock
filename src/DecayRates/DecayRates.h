/*
 * Decay rate calculator for a sterile neutrino below 500MeV in a minimal model
 * Other decay channels could be added
	 * Usage: create a generator Decay(f,f,f,f) defining sterile mass and the mixing elements.
	 *
 * Author: Tommaso Boschi
 */

#ifndef decayrates_H
#define decayrates_H

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
		void SetEnhancement(std::Channel = "ALL", double K = 1.0);

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

		double M_Neutrino = 0.0;
		double M_Electron = Tools::Const::fMElectron;
		double M_Muon = Tools::Const::fMMuon;
		double M_Pion = Tools::Const::fMPion;
		double M_Pion0 = Tools::Const::fMPion0;
		double M_Kaon = Tools::Const::fMKaon;
		double M_Kaon0 = Tools::Const::fMKaon0;

		std::map<std::string, ChannelName> mapChannel;
		std::map<std::string, double> mapEnhance;
};

#endif
