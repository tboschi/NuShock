/*
 * Detector class
 * Collections of relevant properties of near detector
 *
 * Author: Tommaso Boschi
 */

#ifndef detector_H
#define detector_H

#include <iostream>
#include <fstream>

//Boost lib include
#include "boost/random.h"

//ROOT include
#include "TH1.h"
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"

//GENIE include
#include "GHepParticle.h"
#include "Constants.h"

#include "Tools.h"
#include "FluxDriver.h"
#include "Flux.h"
#include "DecayRates.h"

class Detector
{
	public:
		Detector(std::string ConfigFile);

		double GetBaseline();
		double GetLength();
		double GetWeight();
		double GetPOT();
		double GetEfficiency();

	private:
		std::map<std::string, double> mapDetector;
};
	
