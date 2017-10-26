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
		bool MakeFlux(double M_Sterile);
		void MakeMuonComponent(bool Neutrino, Flux &sxFlux, double M_Sterile);
		void MakeElecComponent(bool Neutrino, Flux &sxFlux, double M_Sterile);
		//void MakeStandardFlux();
		//double SampleEnergy();

		double GetRange();
		double GetRangeStart();
		double GetRangeEnd();
		int GetBinNumber();
		double GetIntensityNeut(double Energy, double Ue = 1.0, double Um = 1.0, double Ut = 1.0);
		double GetIntensityAnti(double Energy, double Ue = 1.0, double Um = 1.0, double Ut = 1.0);
		void SetBaseline(double Baseline);
		void SetPOT(double POT);
		void SetArea(double Area);

		bool IsChanged(double M_Sterile);

		/*
		TH1D* GetTotalMn();
		TH1D* GetPionMn();
		TH1D* GetKaonMn();
		TH1D* GetKaon0Mn();
		TH1D* GetMuonMn();
	
		TH1D* GetTotalEn();
		TH1D* GetPionEn();
		TH1D* GetKaonEn();
		TH1D* GetKaon0En();
		TH1D* GetMuonEn();

		TH1D* GetTotalMa();
		TH1D* GetPionMa();
		TH1D* GetKaonMa();
		TH1D* GetKaon0Ma();
		TH1D* GetMuonMa();
	
		TH1D* GetTotalEa();
		TH1D* GetPionEa();
		TH1D* GetKaonEa();
		TH1D* GetKaon0Ea();
		TH1D* GetMuonEa();
		*/

	private:
		int BinNumber;
		double RangeStart, RangeEnd;

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
		TH1D* hTotalMn, * hTotalMa;
		TH1D* hPionMn,  * hPionMa;
		TH1D* hKaonMn,  * hKaonMa;
		TH1D* hKaon0Mn, * hKaon0Ma;
		TH1D* hMuonMn,  * hMuonMa;

		TH1D* hTotalEn, * hTotalEa;
		TH1D* hPionEn,  * hPionEa;
		TH1D* hKaonEn,  * hKaonEa;
		TH1D* hKaon0En, * hKaon0Ea;
		TH1D* hMuonEn,  * hMuonEa;

		double M_Sterile_prev;

		const double M_Electron;
		const double M_Muon;
		const double M_Pion;
		const double M_Pion0;
		const double M_Kaon;
		const double M_Kaon0;
};

#endif
