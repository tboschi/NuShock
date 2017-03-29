/*
 * Flux creator for Monte Carlo using published fluxes
 *
 * Author: Tommaso Boschi
 */

#ifndef fluxdriver_H
#define fluxdriver_H

#include <iostream>
#include <fstream>

//Boost lib include
#include "boost/random.h"

//ROOT include
#include "TH1.h"
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"

//GENIE include
#include "GHepParticle.h"
#include "Constants.h"

#include "Tools.h"
#include "Flux.h"

class FluxDriver
{
	public:
		FluxDriver(std::string SourceName);	//Load from ROOT file
		~FluxDriver();

		TH1D *GetHist();
		void MakeSterileFlux(double M_Sterile);
		double SampleEnergy(Flux *StdFlux = 0, Flux *HeavyFlux = 0)

		void SetTotalFlux(double X);
	
	private:
		TFile *SourceFile;

		TH1F *hTotalFlux, sTotalFlux;
		TH1F *hMuonPion, *sMuonPion;
		TH1F *hMuonKaon, *sMuonKaon;
		TH1F *hElectronPion, *sElectronPion;
		TH1F *hElectronKaon, *sElectronKaon;
		TH1F *hElectronKaon3, *sElectronKaon3;
		TH1F *hMuonKaonOther, *sMuonKaonOther;

		double M_Electron = Tools::Constants::fMElectron;
		double M_Muon = Tools::Constants::fMMuon;
		double M_Pion = Tools::Constants::fMPion;
		double M_Kaon = Tools::Constants::fMKaon;
}

#endif
