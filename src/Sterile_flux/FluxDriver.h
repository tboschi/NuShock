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

#include "Flux.h"

class FluxDriver
{
	public:
		FluxDriver(std::string SourceName);	//Load from ROOT file

		double GetComponents;

		void SetTotalFlux(double X);

		double ShrockFactor(double M_Meson, double M_Lepton, double M_Sterile); 
		double ShrockRho(double X double Y);
		double ShrockFM(double X, double Y);
		double ShrockLambda(double X, double Y, double Z);
	
	private:
		Flux *SterileFlux;

		TFile *SourceFile;

		TH1F *hTotalFlux, sTotalFlux;
		TH1F *hMuonPion, *sMuonPion;
		TH1F *hMuonKaon, *sMuonKaon;
		TH1F *hElectronPion, *sElectronPion;
		TH1F *hElectronKaon, *sElectronKaon;
		TH1F *hElectronKaon3, *sElectronKaon3;
		TH1F *hMuonKaonOther, *sMuonKaonOther;

		double M_Pion;
		double M_Kaon;
		double M_Muon;
		double M_Electron;
}

#endif
