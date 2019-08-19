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

#include "tools.h"
#include "src/detector/Detector.h"

class Detector;

class Tracker : public Detector
{
	public:
		Tracker(std::string ConfigFile, std::string mod = "");

		bool Reconstruct(Particle &P);
		void Vertex(Particle &P);
		void Smearing(Particle &P);

		void Length(Particle &P);

		double GammaDecay(const Particle &P);
		double CriticalEnergy(const Particle &P);
		double CriticalEnergy(Detector::Material Element);
		double RadiationLength(const Particle &P, bool Nuclear = false);
		double RadiationLength(Detector::Material Element, bool Nuclear = false);
		double EnergyLoss(const Particle &P, bool &Contained);
		double BetheLoss(const Particle &P, Material Target);
		double Bethe(const Particle &P, double Density, double I, int Z, int A);

		bool IsDecayed(const Particle &P, double dx);
		bool IsDetectable(const Particle &P, bool print = false);
		void Pi0Decay(Particle &Pi0, Particle &PA, Particle &PB);

		void Focus(Particle &P);

	private:
		//TRandom3 *GenMT;
		std::string module;
		TFile *FuncFile;
		TH2D *hhFunc;
		TH1D *hTemp, *hEfficiency;
		bool Eff;
};

#endif
