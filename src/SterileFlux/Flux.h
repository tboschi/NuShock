/*
 * Flux class, container of various components as root histogram
 * 
 * Author: Tommaso Boschi
 */

#ifndef FLUX_H
#define FLUX_H

#include <iostream>

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

		TH1D* GetTotal() const;
		TH1D* GetPion() const;
		TH1D* Get2Pion() const;
		TH1D* GetKaon() const;
		TH1D* GetKaon0() const;
		TH1D* GetCharm() const;
		TH1D* GetMuon() const;
		TH1D* GetTauE() const;
		TH1D* GetTauM() const;

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
