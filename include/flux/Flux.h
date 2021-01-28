#ifndef FLUX_H
#define FLUX_H

#include <string>
#include <unordered_map>
#include <memory>

#include "TH1.h"
#include "TFile.h"
#include "TKey.h"

struct Flux
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
		undefined
	};

	using Component = std::unordered_map<Parent, std::shared_ptr<TH1D> >;

	static Parent fromString(std::string comp) {
		std::transform(comp.begin(), comp.end(), comp.begin(),
				[](unsigned char c) -> unsigned char { return std::tolower(c); });

		//if (comp == "total")
			//return Parent::Total;
		if (comp == "pion")
			return Parent::Pion;
		else if (comp == "ppion" || comp == "2pion")
			return Parent::PPion;
		else if (comp == "kaon")
			return Parent::Kaon;
		else if (comp == "kaon0")
			return Parent::Kaon0;
		else if (comp == "charm")
			return Parent::Charm;
		else if (comp == "muon")
			return Parent::Muon;
		else if (comp == "taue")
			return Parent::TauE;
		else if (comp == "taum")
			return Parent::TauM;
		else
			return Parent::undefined;
	}

	static std::string toString(const Parent &comp) {
		switch (comp) {
			//case Parent::Total:
				//return "Total";	
			case Parent::Pion:
				return "pion";
			case Parent::PPion:
				return "ppion";
			case Parent::Kaon:
				return "kaon";
			case Parent::Kaon0:
				return "kaon0";
			case Parent::Charm:
				return "charm";
			case Parent::Muon:
				return "muon";
			case Parent::TauE:
				return "taue";
			case Parent::TauM:
				return "taum";
			default:
				return "undefined";
		}
	}

	// return map of histograms from file
	static Component fromROOT(const std::string &file) {
		TFile infile(file.c_str(), "READ");
		if (infile.IsZombie())
			throw std::invalid_argument("Flux: file does not exist\n");

		TIter next(infile.GetListOfKeys());
		TKey* k;
		Component hists;
		while ((k = static_cast<TKey*>(next()))) {
			std::string name = k->GetName();
			std::shared_ptr<TH1D> hflux(static_cast<TH1D*>(infile.Get(name.c_str())));
			// remove leading 'h' from name 
			if (name.find_first_of("h") == 0)
				name.erase(0, 1);

			// if unknown or no object, skip
			if (fromString(name) == Parent::undefined || !hflux)
				continue;

			hflux->SetDirectory(0);
			hists[fromString(name)] = hflux;
		}

		if (!hists.size())
			throw std::logic_error("No histogram collected! Very bad\n");

		/*
		if (!fx.hists.count(Parent::Total)) {
			//make a copy of any other histogram
			fx.hists[Parent::Total] = std::shared_ptr<TH1D>(*hists.begin());
			fx.hists[Parent::Total]->Reset("ICES");
			for (const auto &ih : hists)
				if (ih != Parent::Total)
					hists[Parent::Total]->Add(ih.second.get());
		}
		*/

		return hists;
	}

	static std::vector<Flux::Parent> All() {
		return { Flux::Parent::Pion, Flux::Parent::PPion, Flux::Parent::Kaon,
			Flux::Parent::Kaon0, Flux::Parent::Charm, Flux::Parent::Muon,
			Flux::Parent::TauE, Flux::Parent::TauM };
	}


};

#endif
