/*
 * Flux creator for Monte Carlo using published fluxes
 *
 * Author: Tommaso Boschi
 */

#ifndef fluxdriver_H
#define fluxdriver_H

#include <iostream>
#include <fstream>
#include <cstring>
#include <sstream>

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
		FluxDriver(std::string ConfigFlux);
		~FluxDriver();

		TH1D *GetHist();
		bool MakeSterileFlux(double M_Sterile, double U_e, double U_m, double U_t);
		void MakeMuonComponent(Flux &sxFlux, double M_Sterile, double U_e, double U_m, double U_t);
		void MakeElecComponent(Flux &sxFlux, double M_Sterile, double U_e, double U_m, double U_t);
		//void MakeStandardFlux();
		double SampleEnergy();

		double GetStartRange();
		double GetEndRange();
		double GetBinNumber();
		double GetIntensity(double Energy);
		void SetBaseline(double Baseline);
		void SetPOT(double POT);
		void SetArea(double Area);

		bool IsChanged(double M_Sterile, double U_e, double U_m, double U_t);

		TH1D* GetTotal();
		TH1D* GetPion();
		TH1D* GetKaon();
		TH1D* GetKaon0();
		TH1D* GetMuon();

		TH1D* GetTotalOriginal();
		TH1D* GetPionOriginal();
		TH1D* GetKaonOriginal();
		TH1D* GetKaon0Original();
		TH1D* GetMuonOriginal();
	
	private:
		double M_Sterile_prev, U_e_prev, U_m_prev, U_t_prev;

		TFile *SourceFile;
		TFile *KineFile;
		bool Kine;

		//Get fluxes from file
		Flux* fxNuMuon;
		Flux* fxNuMuonBar;
		Flux* fxNuElectron;
		Flux* fxNuElectronBar;

		//Kinematic factors
		TH1D *hTemp; 
		TH1D *hMuonMuon; 
		TH1D *hMuonElec;
		TH1D *hKaonMuon;
		TH1D *hKaonElec;
		TH1D *hKaon0Muon;
		TH1D *hKaon0Elec;

		//Output fluxes
		TH1D* hTotalSterile;
		TH1D* hPionSterile;
		TH1D* hKaonSterile;
		TH1D* hKaon0Sterile;
		TH1D* hMuonSterile;

		TH1D* hTotalStandard;
		TH1D* hPionStandard;
		TH1D* hKaonStandard;
		TH1D* hKaon0Standard;
		TH1D* hMuonStandard;

		const double M_Electron;
		const double M_Muon;
		const double M_Pion;
		const double M_Pion0;
		const double M_Kaon;
		const double M_Kaon0;
};

#endif
