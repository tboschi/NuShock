/*
 * Flux creator for Monte Carlo using published fluxes
 *
 * Author: Tommaso Boschi
 */

#ifndef fluxdriver_H
#define fluxdriver_H

#include <iostream>
#include <fstream>

//ROOT include
#include "TH1.h"
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"

#include "Tools.h"
#include "Flux.h"

class FluxDriver
{
	public:
		FluxDriver(std::string ConfigFlux)
		~FluxDriver();

		TH1D *GetHist();
		void MakeSterileFlux(double M_Sterile, double U_e, double U_m, double U_t);
		double SampleEnergy();

		TH1D* GetSterile();
		TH1D* GetPion();
		TH1D* GetKaon();
		TH1D* GetKaon0();
		TH1D* GetMuon();
	
	private:
		TFile *SourceFile;

		Flux* fxNuMuon;
		Flux* fxNuMuonBar;
		Flux* fxNuElectron;
		Flux* fxNuElectronBar;

		TH1D* hTotalSterile;
		TH1D* hPionSterile;
		TH1D* hKaonSterile;
		TH1D* hKaon0Sterile;
		TH1D* hMuonSterile;

		double M_Sterile;
		double U_e, U_m, U_t;
		const double M_Electron = Const::fMElectron;
		const double M_Muon = Const::fMMuon;
		const double M_Pion = Const::fMPion;
		const double M_Kaon = Const::fMKaon;
};

#endif
