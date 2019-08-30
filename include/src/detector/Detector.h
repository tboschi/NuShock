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

#include "tools.h"
#include "physics.h"

class Detector
{
	public:

		enum Material
		{
			undefined, 

			LAr,
			LiquidArgon = LAr,

			GasAr,
			GasseousArgon = GasAr,

			Fe,
			Iron = Fe,

			Pb,
			Lead = Pb,
		};

		Detector(std::string configName, std::string mod = "");

		std::vector<std::string> ListKey();
		std::vector<std::string> ListChannel();
		double Get(std::string Key);
		Detector::Material GetMaterial(std::string Key);
		Detector::Material FindMaterial(std::string Key);

		double Efficiency(const Neutrino &Nu);
		double Efficiency(double Energy, double Mass, std::string channel);
		void SetEfficiency(std::string key);

		double XsizeLAr();
		double XstartLAr();
		double XendLAr();
		double YsizeLAr();
		double YstartLAr();
		double YendLAr();
		double ZsizeLAr();
		double ZstartLAr();
		double ZendLAr();

		double XsizeFGT();
		double XstartFGT();
		double XendFGT();
		double YsizeFGT();
		double YstartFGT();
		double YendFGT();
		double YcentreFGT();
		double ZsizeFGT();
		double ZstartFGT();
		double ZendFGT();
		double ZcentreFGT();

		double Xstart();
		double Xend();
		double Xsize();
		double Ystart();
		double Yend();
		double Ysize();
		double Zsize();
		double Zstart();
		double Zend();

		double AreaLAr();
		double AreaFGT();
		double Area();
		double VolumeLAr();
		double VolumeFGT();
		double Volume();
		double RatioLAr();
		double RatioFGT();
		double Ratio();
		double WeightLAr();
		double WeightFGT();
		double Weight();

		double Radius();
		double AngularAcceptance();

		bool IsInsideLAr(const Particle &P);
		bool IsInsideFGT(const Particle &P);
		bool IsInside(const Particle &P);

		double DecayProb(Neutrino &N);
		double DecayProb(const Particle &P, double Total, double Branch);

	protected:
		TRandom3 *GenMT;

	private:
		std::string module;
		CardDealer dc;

		TFile *funcFile;
		std::map<std::string, TH2D*> mhFunc;
		bool effSet;

		std::map<std::string, double> mDetector;
		std::map<std::string, std::string> mEfficiencyD, mEfficiencyM;
		std::map<std::string, Material> mMaterial;
};

#endif
