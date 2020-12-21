/*
 * Flux class, container of various components as root histogram
 * 
 * Author: Tommaso Boschi
 */

#ifndef FLUX_H
#define FLUX_H

#include <iostream>
#include <string>
#include <memory>
#include <map>

#include "TObject.h"
#include "TFile.h"
#include "TH1D.h"

#include "tools.h"

class Flux
{
	public:
		enum Component
		{
			Total,
			Pion,
			PPion,
			Kaon,
			Kaon0,
			Charm,
			Muon,
			TauE,
			TauM,
			_undefined = -1;
		};

		Flux(std::string HistFile);
		Flux(const Flux & f);	//copy ctor
		~Flux();

		void CloneCopy(TH1D*& T, TObject* X);
		void CloneCopy(TH1D*& T, TH1D* X);

		TH1D* Get(Component Name) const;

		void Combine();
		void Add(Component Name);
		void Scale(double X);
		void Scale(Component Name, double X);
		bool Stretch(Component Name, double Sx, double Ex);

		double RangeStart();
		double RangeEnd();
		double BinNumber();
		double BinWidth();

	private:
		std::map<Component, TH1D*> hists;
};

#endif
