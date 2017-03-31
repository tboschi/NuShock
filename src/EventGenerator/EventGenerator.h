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
#include <sstream>
#include <string>
#include <vector>

//ROOT include
#include "TH1.h"
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TRandom3.h"

#include "Tools.h"
#include "FluxDriver.h"
#include "Flux.h"
#include "DecayRates.h"
#include "Detector.h"

class EventGenerator
{
	public:
		EventGenerator(std::string SMConfig, std::string FluxConfig, std::string DetectorConfig);
		~EventGenerator();

		double DrawEnergy();	//returns energy from distribution and also store components
		double SetEnergy(double X);
		void MakeFlux();
		double Probability(std::string Channel, double ESterile);
		double RandomDetectionEvent(std::string Channel);
		double NumberOfDetected(std::string Channel);
		double RandomEvent();
		std::string RandomChannel(double Energy);

		double GetMSterile();
		double GetUe();
		double GetUm();
		double GetUt();

		void SetMSterile(double X);
		void SetUe(double X);
		void SetUm(double X);
		void SetUt(double X);

	private:
		double M_Sterile, E_Sterile;
		double U_e, U_m, U_t;

		FluxDriver *TheFlux;
		Flux *SterileEnergy;
		Flux *StandardEnergy;
		Decay *TheGamma;
		Detector *TheBox;

		TFile *SourceFile;

		TH1F *hTotalFlux, sTotalFlux;
		TH1F *hMuonPion, *sMuonPion;
		TH1F *hMuonKaon, *sMuonKaon;
		TH1F *hElectronPion, *sElectronPion;
		TH1F *hElectronKaon, *sElectronKaon;
		TH1F *hElectronKaon3, *sElectronKaon3;
		TH1F *hMuonKaonOther, *sMuonKaonOther;

		const double M_Electron = Const::fMElectron;
		const double M_Muon = Const::fMMuon;
		const double M_Pion = Const::fMPion;
		const double M_Kaon = Const::fMKaon;

		TRandom3 *GenMT;
};

#endif
