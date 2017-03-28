/*
 * Event Generator for sterile neutrinos
 * Generates number of events given detector, flux and decay modes
 *
 * Author: Tommaso Boschi
 */

#ifndef eventgenerator_H
#define eventgenerator_H

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
#include "FluxDriver.h"
#include "Flux.h"
#include "DecayRates.h"

class EventGenerator
{
	public:
		EventGenerator(double M_Sterile, double U_e, double U_m, double U_t);	//
		~EventGenerator();

		double Probability();
		double Efficiency();
		double GetFlux();

		void MakeSterileFlux(double M_Sterile)
		Flux * SampleEnergy()

		void SetTotalFlux(double X);
	
	private:
		double M_Sterile, E_Sterile;
		double U_e, U_m, U_t;

		Flux *SterileFlux;
		Decay *TheGamma;

		TFile *SourceFile;

		TH1F *hTotalFlux, sTotalFlux;
		TH1F *hMuonPion, *sMuonPion;
		TH1F *hMuonKaon, *sMuonKaon;
		TH1F *hElectronPion, *sElectronPion;
		TH1F *hElectronKaon, *sElectronKaon;
		TH1F *hElectronKaon3, *sElectronKaon3;
		TH1F *hMuonKaonOther, *sMuonKaonOther;

		double M_Electron = Tools::Constants::fMElectron;
		double M_Muon = Tools::Constants::fMMuon;
		double M_Pion = Tools::Constants::fMPion;
		double M_Kaon = Tools::Constants::fMKaon;
}

#endif
