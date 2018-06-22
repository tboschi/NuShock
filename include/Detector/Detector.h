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

#include "Tools.h"
#include "Physics.h"

class Detector
{
	public:
		enum Coupling
		{
			E,
			M,
			T
		};

		enum Material
		{
			undefined, 

			LAr,
			LiquidArgon = LAr,

			GasAr,
			GasseousArgon = GasAr,

			Fe,
			Iron = Fe
		};

		Detector(std::string ConfigFile);

		std::vector<std::string> ListKey();
		std::vector<std::string> ListChannel();
		double Get(std::string Key);
		Detector::Material GetMaterial(std::string Key);
		Detector::Material FindMaterial(std::string Key);

		double Efficiency(double Energy, double Mass);
		void SetEfficiency(std::string Channel, Coupling U);

		double Xsize();
		double Xstart();
		double Xend();
		double Ysize();
		double Ystart();
		double Yend();
		double Zsize();
		double Zstart();
		double Zend();

		double DecayProb(Neutrino *N);
		double DecayProb(Particle *N, double Total, double Branch);

	private:
		TRandom3 *GenMT;
		TFile *FuncFile;
		TH2D *hhFunc;

		std::map<std::string, double> mDetector;
		std::map<std::string, std::string> mEfficiencyE, mEfficiencyM, mEfficiencyT;
		std::map<std::string, Material> mMaterial;
};

#endif
