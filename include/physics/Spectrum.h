#ifndef SPECTRUM_H
#define SPECTRUM_H

#include <unordered_map>
#include <memory>

#include "physics/Mixings.h"
#include "physics/Flavours.h"
#include "physics/Neutrino.h"

#include "TH1D.h"

// aggreate of HNL flux and decay rate
// create by a Driver object
// use to combine production and decay rate
// with real flux information
struct Spectrum {
	using Distribution = std::unordered_map<Nu::Flavour, std::shared_ptr<TH1D> >;

	Neutrino N;
	Distribution distribution;

	// create flux
	std::shared_ptr<TH1D> operator()(const Mixing &mix = Mixing()) const;
	double Intensity(double energy, const Mixing &mix = Mixing()) const;
};

#endif
