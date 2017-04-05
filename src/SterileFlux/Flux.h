/*
 * Flux class, container of various components as root histogram
 * 
 * Author: Tommaso Boschi
 */

#ifndef FLUX_H
#define FLUX_H

#include <iostream>

#include "TH1D.H"

class Flux
{
	public:
		Flux(TH1D* Total,
		     TH1D* Pion, 
		     TH1D* Kaon,
		     TH1D* Kaon0,
		     TH1D* Muon);
		Flux(std::string HistFile);

		void CloneTotal(TH1D* Hist);
		void ClonePion(TH1D* Hist);
		void CloneKaon(TH1D* Hist);
		void CloneKaon0(TH1D* Hist);
		void CloneMuon(TH1D* Hist);

		TH1D* GetTotal();
		TH1D* GetPion();
		TH1D* GetKaon();
		TH1D* GetKaon0();
		TH1D* GetMuon();

	private:
		TH1D *hAll;
		TH1D *hPion;
		TH1D *hKaon;
		TH1D *hKaon0;
		TH1D *hMuon;
};

#endif
