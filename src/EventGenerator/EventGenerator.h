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
#include <cmath>

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
#include "Particle.h"

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
		long double EventProbability();
		double EventEfficiency();
		void EventTotalNumber(bool Efficiency = false);
		//Random generators
		std::string RandomChannel();
		bool EventInDetector();
		bool EventDetectable();
		//Kinematics
		int EventKinematics();
		Particle *GetDecayProduct(int i, bool Smear = false);

		//Generate flux to be used as PDF
		void MakeSterileFlux(bool TotalPOT = true);
		void MakeInDetector(bool Efficiency = false);
		//void MakeStandardFlux(bool TotalPOT = true);
		double SampleEnergy(bool Set = true);
		double SampleInDetector(bool Set = true);
		long double FluxIntensity();
		long double BackgroundIntensity();	//Get the background intensity at given energy

		//Statistics
		void SmearVector(TLorentzVector* N, int Pdg);

		//Get function
		std::string GetChannel();
		double GetMass(int Pow = 1);
		double GetEnergy(int Pow = 1);
		double GetEnergyKin(int Pow = 1);
		double GetMomentum(int Pow = 1);
		double GetUe();
		double GetUm();
		double GetUt();

		long double GetSignal();
		long double GetBackground();
		long double GetReducedChi2();

		//Set function
		void SetChannel(std::string Ch = "R");
		void SetMass(double X);
		void SetEnergy(double X);
		void SetEnergyKin(double X);
		void SetUe(double X);
		void SetUm(double X);
		void SetUt(double X);

		void SetSignalNumber(long double X);
		void SetBackgroundNumber(long double X);
		void SetReducedChi2(long double X);

	private:
		std::string sChannel;	//Channel is set globally in class

		double M_Sterile, E_Sterile;
		double U_e, U_m, U_t;

		long double lSignal, lBackground, lRedChi2;

		Decay *TheGamma;
		Detector *TheBox;
		FluxDriver *TheFlux;

		TFile *SourceFile;
		TRandom3 *GenMT;

		TH1D *InDetector;

};

#endif
