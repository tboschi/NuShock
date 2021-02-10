#ifndef FLUX_H
#define FLUX_H

#include <iostream>
#include <string>
#include <unordered_map>
#include <memory>

#include "TH1.h"
#include "TFile.h"
#include "TKey.h"

namespace Flux
{
	enum class Parent
	{
		Pion,
		PPion,
		Kaon,
		Kaon0,
		Charm,
		Muon,
		TauE,
		TauM,
	};

	using Component = std::unordered_map<Parent, std::shared_ptr<TH1D> >;

	Parent fromString(std::string comp);
	std::string toString(Parent comp);

	// return map of histograms from ROOT file
	Component fromROOT(std::string file);

	std::vector<Parent> Parents();
};

#endif
