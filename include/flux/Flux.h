/*
 * Flux class, container of various components as root histogram
 * 
 * Author: Tommaso Boschi
 */

#ifndef FLUX_H
#define FLUX_H

#include <iostream>
#include <string>

#include "TObject.h"
#include "TFile.h"
#include "TH1D.h"

#include "tools.h"

class Flux
{
	public:
		enum Hist
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
		};

		Flux(std::string HistFile);
		Flux(const Flux & f);	//copy ctor
		~Flux();

		void CloneCopy(TH1D*& T, TObject* X);
		void CloneCopy(TH1D*& T, TH1D* X);

		TH1D* Get(Hist Name) const;

		void Add();
		void Add(Hist Name);
		void Scale(double X);
		void Scale(Hist Name, double X);
		bool Stretch(Hist Name, double Sx, double Ex);

		double RangeStart();
		double RangeEnd();
		double BinNumber();
		double BinWidth();

	private:
		TH1D *hTotal;
		TH1D *hPion;
		TH1D *hPPion;
		TH1D *hKaon;
		TH1D *hKaon0;
		TH1D *hCharm;
		TH1D *hMuon;
		TH1D *hTauE;
		TH1D *hTauM;
};

#endif
