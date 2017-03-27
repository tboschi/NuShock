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
		TFile *SourceFile;

		TH1F *TotalFlux;
		TH1F *hMuonPion;
		TH1F *hMuonKaon;
		TH1F *hElectronPion;
		TH1F *hElectronKaon;
		TH1F *hElectronKaon3;
		TH1F *hMuonKaonOther;
}

#endif
