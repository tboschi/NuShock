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

#include "Tools.h"
#include "Detector.h"
#include "Physics.h"

#include "Flux/Flux.h"
#include "Flux/Driver.h"

class Engine
{
	public:

		enum Current
		{
			RHC = 0,
			FHC = 1,
		};

		Engine(std::string FluxConfig, unsigned int nFHC, unsigned int nRHC);
		~Engine();

		void BindNeutrino(Neutrino *N, Current Horn, unsigned int ID);

		void MakeFlux();
		void MakeFlux(Current Horn);
		void MakeFlux(Current Horn, unsigned int ID);
		void MakeSampler(Detector *Box);
		void MakeSampler(Detector *Box, Current Horn);
		void MakeSampler(Detector *Box, Current Horn, unsigned int ID);

		double Intensity(Current Horn, unsigned int ID);

		void ScaleBaseline(Detector *Box);
		void ScalePOT(Detector *Box);
		void ScaleArea(Detector *Box);
		void ScaleBaseline(double Baseline);
		void ScalePOT(double POT);
		void ScaleArea(double Area);
		void Scale(double X);
		void Scale(double X, Current Horn);
		void Scale(double X, Current Horn, unsigned int ID);

		double RangeWidth(double &Start, double &End);

		unsigned int vNeutrino(Current Horn);
		Neutrino* vNeutrino(Current Horn, unsigned int i);
		unsigned int vDriver(Current Horn);
		Driver* vDriver(Current Horn, unsigned int i);

	private:

		std::vector<Driver*> vDriverFHC, vDriverRHC;
		std::vector<Neutrino*> vNeutrinoFHC, vNeutrinoRHC;
		std::vector<TH1D*> vSampleNuFHC, vSampleNuRHC;
};

#endif
