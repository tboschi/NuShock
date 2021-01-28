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

#include "flux/Flux.h"

#include "physics/Const.h"
#include "physics/Flavours.h"
#include "physics/Mixings.h"
#include "physics/Channels.h"
#include "physics/Production.h"

class Driver
{
	public:
		using Modifier = std::vector<std::array<double, 4> >;

	private:
		std::unordered_map<Nu::Flavour, Flux::Component> _fxNu;
		std::unordered_map<Flux::Parent, Driver::Modifier> _modifiers;

		std::unordered_map<Nu::Flavour, std::shared_ptr<TH1D> > _distr;

	public:
		Driver(const std::string &card);
		void Init(const CardDealer &cd);
		double Intensity(const Neutrino &N, const Mixing &mix);
		double Intensity(double energy, const Mixing &mix);
		std::shared_ptr<TH1D> Spectrum(const Neutrino &N, const Mixing &mix);
		std::shared_ptr<TH1D> Spectrum(Nu::Flavour flv, const Mixing &mix);

		bool MakeFlux(std::initializer_list<Neutrino> nus,
				const Mixing &mix = Mixing());
		bool MakeFlux(const Neutrino &N, const Mixing &mix = Mixing());

		std::shared_ptr<TH1D> MakeComponent(Production &hnl, Nu::Flavour nu, double mass = 0.);
		std::shared_ptr<TH1D> MakeElectron(Production &heavy, const Flux::Component &fxNu);
		std::shared_ptr<TH1D> MakeMuon(Production &heavy, const Flux::Component &fxNu);
		std::shared_ptr<TH1D> MakeTau(Production &heavy, const Flux::Component &fxNu,
					double mass = 0.);
		double Stretch(std::shared_ptr<TH1D> hist, const Modifier &mod, double mass = 0.);
};

#endif
