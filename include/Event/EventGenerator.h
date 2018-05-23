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
#include "Flux.h"
#include "Decay/DecayRates.h"
#include "Detector/Detector.h"
#include "Scattering/Nucleon.h"

class EventGenerator
{
	public:
		EventGenerator();
		//~EventGenerator();

		//For a finer handle...it shouldn't be necessary
		Detector* GetDetectorPtr();
		Decay* GetDecayPtr();
		FluxDriver* GetFluxDriverPtr();

		double DecayProb();
		double ScatterProb(double Eh);
		double EventEfficiency();
		double DecayNumber(double Energy, bool Efficiency = false);
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

		//Get function
		std::string GetChannel();
		double GetMass(int Pow = 1);
		double GetEnergy(int Pow = 1);
		double GetEnergyKin(int Pow = 1);
		double GetMomentum(int Pow = 1);
		double GetUe(int Pow = 1);
		double GetUm(int Pow = 1);
		double GetUt(int Pow = 1);
		double GetRange(double &Start, double &End);
		int GetBinNumber();

		//Set function
		void SetChannel(std::string Ch = "R", bool Efficiency = false, char Couple = 'M');
		void SetMass(double X);
		void SetEnergy(double X);
		void SetEnergyKin(double X);
		void SetUe(double X, bool GvF = 0);
		void SetUm(double X, bool GvF = 0);
		void SetUt(double X, bool GvF = 0);
		void SyncUu(bool B = 1);

		bool IsChanged();

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

		Decay *TheGamma;
		Detector *TheBox;
		FluxDriver *TheFlux;
		Nucleon *TheProton, *TheNeutron;
		//Nucleon *TheProton_a, *TheNeutron_a;

		TFile *SourceFile;
		TRandom3 *GenMT;

		TVector3* Position;

		TH1D *hSampler;

};

#endif
