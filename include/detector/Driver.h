/*
 * Class to drive neutrinos to the ND
 * Modifies published fluxes to extimate sterile neutrino flux
 *
 * Author: Tommaso Boschi
 */

#ifndef DRIVER_H
#define DRIVER_H

#include <iostream>
#include <string>
#include <fstream>
#include <string>
#include <unordered_map>

//ROOT include
#include "TH1.h"
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"

#include "tools/CardDealer.h"

#include "detector/Flux.h"

#include "physics/Const.h"
#include "physics/Flavours.h"
#include "physics/Mixings.h"

#include "physics/Productions.h"
#include "physics/ProductionRate.h"
#include "physics/Spectrum.h"

// combination of flux constructor
// and neutrino used to build spectrum
class Driver
{
	public:
		using Modifier = std::vector<std::array<double, 4> >;

		Driver(const std::string &card);
		void Init(const CardDealer &cd);
		//double Intensity(const Neutrino &N, const Mixing &mix);
		//double Intensity(double energy, const Mixing &mix);
		//std::shared_ptr<TH1D> Spectrum(const Neutrino &N, const Mixing &mix) const;
		//std::shared_ptr<TH1D> Spectrum(Nu::Flavour flv, const Mixing &mix) const;

		//bool MakeFlux(std::initializer_list<Neutrino> nus,
				//const Mixing &mix = Mixing());
		Spectrum MakeSpectrum(const Neutrino &N, const Mixing &mix = Mixing()) const;

		std::shared_ptr<TH1D> MakeComponent(ProductionRate &hnl, Nu::Flavour nu,
				double mass = 0.) const;

	private:
		std::shared_ptr<TH1D> MakeElectron(ProductionRate &heavy, const Flux::Component &fxNu) const;
		std::shared_ptr<TH1D> MakeMuon(ProductionRate &heavy, const Flux::Component &fxNu) const;
		std::shared_ptr<TH1D> MakeTau(ProductionRate &heavy, const Flux::Component &fxNu,
					double mass = 0.) const;
		double Stretch(std::shared_ptr<TH1D> hist, const Modifier &mod, double mass = 0.) const;


		std::unordered_map<Nu::Flavour, Flux::Component> _fxNu;
		std::unordered_map<Flux::Parent, Driver::Modifier> _modifiers;
		// full distribution

		// save infos
		double _mass;
		int _helix;
};

#endif
