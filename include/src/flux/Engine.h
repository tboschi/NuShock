/*
 * Flux engine, to handle more Driver at the same time
 *
 * Author: Tommaso Boschi
 */

#ifndef ENGINE_H
#define ENGINE_H

#include <iostream>
#include <string>
#include <fstream>
#include <cstring>
#include <sstream>

#include "tools.h"
#include "detector.h"
#include "physics.h"

#include "Flux.h"
#include "Driver.h"

class Engine
{
	public:

		enum Current
		{
			RHC = 0,
			FHC = 1,
			both = 2,
		};

		Engine(std::string fc);
		~Engine();
		void Reset();

		void BindNeutrino(std::string uuid, Neutrino &N, Current horn);
		void ReleaseNeutrino(std::string uuid);
		Neutrino& GetNeutrino(std::string uuid);

		void MakeFlux();
		void MakeFlux(Current horn);
		void MakeFlux(std::string uuid);

		void SampleEnergy(std::map<std::string, double> &mE,
				  std::map<std::string, double> &mI);
		void SampleEnergy(std::map<std::string, double> &mE,
				  std::map<std::string, double> &mI,
				  Current horn);
		double SampleEnergy(std::string uuid);

		double MakeSampler(Detector *box,
				   double ue = -1, double um = -1, double ut = -1);
		double MakeSampler(Detector *box, Current horn,
				   double ue = -1, double um = -1, double ut = -1);
		double MakeSampler(Detector *box, std::map<std::string, double> &mInt,
				   double ue = -1, double um = -1, double ut = -1);
		double MakeSampler(Detector *box, std::map<std::string, double> &mInt,
				   Current horn,
				   double ue = -1, double um = -1, double ut = -1);
		double MakeSampler(Detector *box, std::string uuid,
				   double ue = -1, double um = -1, double ut = -1);

		double DecayNumber(Detector *box);
		double DecayNumber(Detector *box, Current horn);
		double DecayNumber(Detector *box, std::string uuid);

		double Intensity(std::string uuid);
		double IntensitySample(std::string uuid);
		double IntensitySample(std::string uuid, double Energy);

		void ScaleToDetector(Detector *box);
		void ScaleBaseline(Detector *box);
		void ScalePOT(Detector *box);
		void ScaleArea(Detector *box);
		void ScaleBaseline(double Baseline);
		void ScalePOT(double POT);
		void ScaleArea(double Area);
		void Scale(double X);
		void Scale(double X, Current horn);
		void Scale(double X, std::string uuid);

		double BinNumber();
		double RangeWidth();
		double RangeWidth(double &start, double &end);

		std::string HornName(Current horn);

	private:
		std::string fluxConfig;

		std::map<std::string, Neutrino> mNeutrino;
		std::map<std::string, Neutrino>::iterator iN;

		std::map<std::string, Driver*> mDriver;
		std::map<std::string, Driver*>::iterator iD;

		std::map<std::string, TH1D*> sampleNu;
		std::map<std::string, TH1D*>::iterator iS;
};

#endif
