/*
 * Event Generator for sterile neutrinos
 * Generates number of events given detector, wanted energy and decay modes
 * Developed as a super controller
 *
 * Author: Tommaso Boschi
 */

#ifndef eventgenerator_H
#define eventgenerator_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>
#include <vector>

//ROOT include
#include "TH1.h"
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TLorentzVector.h"

#include "Tools.h"
#include "FluxDriver.h"
#include "Flux.h"
#include "DecayRates.h"
#include "Detector.h"

class EventGenerator
{
	public:
		EventGenerator(std::string SMConfig, std::string DetectorConfig, std::string FluxConfig);
		//~EventGenerator();

		Detector* GetDetectorPtr();
		Decay* GetDecayPtr();
		FluxDriver* GetFluxDriverPtr();

		std::string RandomChannel();
		double Probability(std::string Channel = "ALL");
		bool Detectable(std::string Channel = "R");
		bool RandomDetectionEvent(std::string Channel = "R");
		int SimulateDecay(std::string Channel = "R");
		TLorentzVector *GetDecayProduct(int i);

		void MakeSterileFlux();
		void MakeStandardFlux();
		double SampleEnergy();

		double GetMass();
		double GetEnergy();
		double GetMomentum();
		double GetUe();
		double GetUm();
		double GetUt();

		void SetMass(double X);
		void SetEnergy(double X);
		void SetUe(double X);
		void SetUm(double X);
		void SetUt(double X);

	private:
		double M_Sterile, E_Sterile;
		double U_e, U_m, U_t;

		Decay *TheGamma;
		Detector *TheBox;
		FluxDriver *TheFlux;

		TFile *SourceFile;
		TRandom3 *GenMT;
};

#endif
