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

//Boost lib include
#include "boost/random.h"

//ROOT include
#include "TH1.h"
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"

//GENIE include
#include "GHepParticle.h"
#include "Constants.h"

#include "Tools.h"

static enum ChannelName
{
	ev_undefined,
	ev_ALL,
	ev_nnn,
	ev_nGAMMA,
	ev_nEE,
	ev_nEMU,
	ev_nPI0,
	ev_EPI,
	ev_MUPI,
	ev_nMUMU,
	ev_EKA,
	ev_nKA0
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
