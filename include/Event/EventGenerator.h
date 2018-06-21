/*
 * Event Generator for sterile neutrinos
 * Generates number of events given detector, wanted energy and decay modes
 * Developed as a super controller
 *
 * Author: Tommaso Boschi
 */

#ifndef EVENTGENERATOR_H
#define EVENTGENERATOR_H

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
#include "Flux.h"
#include "Physics.h"
#include "Detector.h"
//#include "Scattering/Nucleon.h"

class EventGenerator
{
	public:
		EventGenerator();
		EventGenerator(Neutrino *N, Driver *Driver, Detector *Box);
		~EventGenerator();

		//For a finer handle...it shouldn't be necessary

		void SetNeutrino(Neutrino* N);
		void SetDriver(Driver* Driver);
		void SetDetector(Detector* Box);
		Neutrino* GetNeutrino();
		Driver* GetDriver();
		Detector* GetDetector();

		double DecayProb();
		double ScatterProb(double Eh);
		double EventEfficiency();
		double DecayNumber(double Energy, bool Efficiency = false);
		double DecayNumber(bool Efficiency = false);
		double ScatterNumber(double Energy, bool Efficiency = false);

		//Random generators
		std::string RandomChannel();
		bool DecayInDetector();
		bool EventDetectable();

		//Kinematics
		bool IsAllowed();
		int EventKinematics();
		Particle *GetDecayProduct(int i, bool Smear = false);
		void Pi0Decay(Particle *Pi0, Particle *&PA, Particle *&PB, bool Smear = false);
		void GeneratePosition();

		//Generate flux to be used as PDF
		void MakeFlux(bool Mass, bool TotalPOT = true);
		void MakeSampler();
		double SampleEnergy(bool Set = true);

		double NeutIntensity(bool Mass);
		double AntiIntensity(bool Mass);
		double FluxIntensity(bool Mass);

		double ScaleXSec(double IntP, double IntN);
		double Variable(double dt);

		//Statistics
		//void SmearVector(TLorentzVector* N, int Pdg);

		void SetEnhancement(double X = 1.0);
		double GetEnhancement();
		double SetLambda(double X);
		void SetUserData(double X);
		double GetUserData();

	private:
		std::string sChannel;	//Channel is set globally in class

		double M_Sterile, E_Sterile;
		double M_Sterile_prev, E_Sterile_prev;
		double Ue, Um, Ut;
		double fTotalXSec;

		double fEnhance, fEnhance_prev;
		double fUserData, fUserData_prev;

		bool Sync;

		Neutrino *TheN;
		Driver *TheFlux;
		Detector *TheBox;
		//Nucleon *TheProton, *TheNeutron;
		//Nucleon *TheProton_a, *TheNeutron_a;

		TFile *SourceFile;
		TRandom3 *GenMT;

		TVector3* Position;

		TH1D *hSampler;

};

#endif
