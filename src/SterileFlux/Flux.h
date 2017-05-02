/*
 * Flux class, container of various components as root histogram
 * 
 * Author: Tommaso Boschi
 */

#ifndef FLUX_H
#define FLUX_H

#include <iostream>

#include "TFile.h"
#include "TH1D.h"

class Flux
{
	public:
		Flux(TH1D* Total,
		     TH1D* Pion, 
		     TH1D* Kaon,
		     TH1D* Kaon0,
		     TH1D* Muon);
		Flux(std::string HistFile);
		Flux(const Flux & f);	//copy ctor
		~Flux();

		void CloneTotal(TH1D* Hist);
		void ClonePion(TH1D* Hist);
		void CloneKaon(TH1D* Hist);
		void CloneKaon0(TH1D* Hist);
		void CloneMuon(TH1D* Hist);

		TH1D* GetTotal() const;
		TH1D* GetPion() const;
		TH1D* GetKaon() const;
		TH1D* GetKaon0() const;
		TH1D* GetMuon() const;

	private:
		TH1D *hTotal;
		TH1D *hPion;
		TH1D *hKaon;
		TH1D *hKaon0;
		TH1D *hMuon;
};

#endif
