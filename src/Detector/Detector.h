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
#include <map>
#include <vector>

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

static struct EnergyEfficiency
{
	double E;
	double f;
}

class Detector
{
	public:
		Detector(std::string ConfigFile);

		void ListKey();
		double GetElement(std::string Key);
		double Efficiency(std::string Channel, double Energy);

	private:
		std::map<std::string, double> mapDetector;
		std::map<std::string, std::vector<EnergyEfficiency> > mapEfficiency;
};
