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

class Flux
{
	public:
		Flux(std::string HistFile);
		Flux(const Flux & f);	//copy ctor
		~Flux();

		void CloneCopy(TH1D*& T, TObject* X);
		void CloneCopy(TH1D*& T, TH1D* X);

		TH1D* Get(std::string T) const;

	private:
		TH1D *hTotal;
		TH1D *hPion;
		TH1D *h2Pion;
		TH1D *hKaon;
		TH1D *hKaon0;
		TH1D *hCharm;
		TH1D *hMuon;
		TH1D *hTauE;
		TH1D *hTauM;
};

#endif
