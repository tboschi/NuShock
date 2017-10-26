/*
 * Detector class
 * Collections of relevant properties of near detector
 *
 * Author: Tommaso Boschi
 */

#ifndef DETECTOR_H
#define DETECTOR_H

#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include <string>
#include <sstream>
#include <iomanip>

//ROOT include
#include "TRandom3.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
//#include "TTree.h"
//#include "TMath.h"

#include "Tools.h"
#include "Particle.h"
//#include "FluxDriver.h"
//#include "Flux.h"
//#include "DecayRates.h"

struct EnergyEfficiency
{
	double E;
	double f;
};

class Detector
{
	public:
		Detector(std::string ConfigFile, TRandom3 *Random);

		std::vector<std::string> ListKey();
		std::vector<std::string> ListChannel();
		double GetElement(std::string Key);
		//double Efficiency(std::string Channel, double Energy);
		double Efficiency(double Energy, double Mass);
		double SetEfficiency(std::string Channel, char Couple);

		void SignalSmearing(Particle *P);
		void TrackLength(Particle *P);
		double GammaDecay();
		double CriticalEnergy();
		double RadiationLength(bool Nuclear = false);
		double EnergyLoss(Particle *P, bool &Contained);
		double BetheLoss(Particle *P, double Target);

		bool IsDecayed(Particle *P, double dx);
		bool IsDetectable(Particle *P);	
		//bool IsDetectable(int Pdg, int Charge, double Ekin);	
		bool IsInside(Particle *P);
		bool IsContained(Particle *P);	//B2B = 1 can cross z wall, B2B = 0 can't cross z wall
		double GetXsize();
		double GetXstart();
		double GetXend();
		double GetYsize();
		double GetYstart();
		double GetYend();
		double GetZsize();
		double GetZstart();
		double GetZend();
	private:
		TRandom3 *RanGen;
		TFile *FuncFile;
		TH2D *hhFunc;
		TH1D *hTemp, *hEfficiency;
		bool Eff;

		std::map<std::string, double> mapDetector;
		std::map<std::string, std::string> mapEfficiencyE, mapEfficiencyM;
};

#endif
