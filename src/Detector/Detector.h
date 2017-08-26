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
		Detector(std::string ConfigFile);

		std::vector<std::string> ListKey();
		std::vector<std::string> ListChannel();
		double GetElement(std::string Key);
		//double Efficiency(std::string Channel, double Energy);
		double Efficiency(std::string Channel, double Energy, double Mass);
		double Background(std::string Channel, double Energy);
		double EnergySigma(std::string Channel, double Energy);

		double EnergySigma(Particle *P);
		void SignalSmearing(TRandom3 *RanGen, Particle *P);
		double TrackLength(TRandom3 *RanGen, Particle *P);
		double InteractionLength(int Target);
		double EnergyLoss(int Target, double Beta, double Mass);
		bool IsDetectable(Particle *P);	
		bool IsDetectable(int Pdg, int Charge, double Ekin);	
		bool IsInside(Particle *P);
		bool IsInside(TVector3 &P);
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

		TFile *EffFile;
		TH1D *hTemp, *hEfficiency;
		bool Eff;

		std::map<std::string, double> mapDetector;
		std::map<std::string, std::string> mapEfficiency;
};

#endif
