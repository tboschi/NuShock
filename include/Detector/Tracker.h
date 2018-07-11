/*
 * Detector class
 * Collections of relevant properties of near detector
 *
 * Author: Tommaso Boschi
 */

#ifndef TRACKER_H
#define TRACKER_H

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

#include "Tools.h"
#include "Detector/Detector.h"

class Tracker : public Detector
{
	public:
		Tracker(std::string ConfigFile);

		void TrackReconstruct(Particle *&P);
		void TrackVertex(Particle *&P);
		void TrackSmearing(Particle *&P);
		void TrackLength(Particle *&P);
		double GammaDecay();
		double CriticalEnergy();
		double RadiationLength(bool Nuclear = false);
		double EnergyLoss(Particle *P, bool &Contained);
		double BetheLoss(Particle *P, Material Target);
		double Bethe(Particle *P, double Density, double I, int Z, int A);

		bool IsDecayed(Particle *P, double dx);
		bool IsDetectable(Particle *P);	
		void Pi0Decay(Particle *Pi0, Particle *&PA, Particle *&PB);

		void Focus(Particle *P);

	private:
		//TRandom3 *GenMT;
		TFile *FuncFile;
		TH2D *hhFunc;
		TH1D *hTemp, *hEfficiency;
		bool Eff;
};

#endif
