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

		//For a finer handle...it shouldn't be necessary
		Detector* GetDetectorPtr();
		Decay* GetDecayPtr();
		FluxDriver* GetFluxDriverPtr();

		//MC stuff
		double EventProbability();
		double EventEfficiency(double Efficiency = -1.0);
		double EventTotalNumber(double Efficiency = -1.0);
		//Random generators
		std::string RandomChannel();
		bool EventInDetector();
		bool EventDetectable();
		//Kinematics
		int EventKinematics();
		TLorentzVector *GetDecayProduct(int i);

		//Generate flux to be used as PDF
		void MakeSterileFlux(bool TotalPOT = false);
		void MakeStandardFlux(bool TotalPOT = false);
		double SampleEnergy();
		double FluxIntensity();

		//Get function
		bool IsChanged();
		std::string GetChannel();
		double GetMass();
		double GetEnergy();
		double GetMomentum();
		double GetUe();
		double GetUm();
		double GetUt();

		//Set function
		void SetChannel(std::string Ch = "R");
		void SetMass(double X);
		void SetEnergy(double X);
		void SetUe(double X);
		void SetUm(double X);
		void SetUt(double X);

	private:
		std::string sChannel;	//Channel is set globally in class

		double M_Sterile, E_Sterile;
		double U_e, U_m, U_t;
		double M_Sterile_prev, E_Sterile_prev;	//Previous declaration
		double U_e_prev, U_m_prev, U_t_prev;

		Decay *TheGamma;
		Detector *TheBox;
		FluxDriver *TheFlux;

		TFile *SourceFile;
		TRandom3 *GenMT;
};

#endif
