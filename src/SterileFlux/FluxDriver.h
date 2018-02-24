/*
 * Flux creator for Monte Carlo using published fluxes
 *
 * Author: Tommaso Boschi
 */

#ifndef fluxdriver_H
#define fluxdriver_H

#include <iostream>
#include <string>
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

		void CloneCopy(TH1D*& T, TObject* X);
		bool MakeFlux(double M_Sterile);
		void MakeElecComponent(bool Neutrino, Flux &sxFlux, double M_Sterile);
		void MakeMuonComponent(bool Neutrino, Flux &sxFlux, double M_Sterile);
		void MakeTauComponent(bool Neutrino, Flux &sxFlux, double M_Sterile);
		//void MakeStandardFlux();
		//double SampleEnergy();

		double GetRange();
		double GetRangeStart();
		double GetRangeEnd();
		int GetBinNumber();
		double GetIntensity(double Energy, bool NvA, bool Uu);
		void SetBaseline(double Baseline);
		void SetPOT(double POT);
		void SetArea(double Area);

		bool IsChanged(double M_Sterile);

		void SetUe(double X);
		void SetUm(double X);
		void SetUt(double X);
		double GetUe();
		double GetUm();
		double GetUt();

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
		unsigned int BinNumber;
		double RangeStart, RangeEnd;
		double U_e, U_m, U_t;

		TFile *SourceFile;
		TFile *KineFile;
		bool Kine;

		//Get fluxes from file
		Flux* fxNuElectron;
		Flux* fxNuElectronBar;
		Flux* fxNuMuon;
		Flux* fxNuMuonBar;
		Flux* fxNuTau;
		Flux* fxNuTauBar;

		//Kinematic factors
		TH1D *hTemp; 
		TH1D *hMuonElec;
		TH1D *hMuonMuon; 
		TH1D *hKaonElec;
		TH1D *hKaonMuon;
		TH1D *hKaon0Elec;
		TH1D *hKaon0Muon;
		TH1D *hTauEElec;
                TH1D *hTauETau;
                TH1D *hTauMMuon;
                TH1D *hTauMTau;


		//Output fluxes
		TH1D* hTotalEn, * hTotalEa;
		TH1D* hPionEn,  * hPionEa;
		TH1D* hKaonEn,  * hKaonEa;
		TH1D* hKaon0En, * hKaon0Ea;
		TH1D* hMuonEn,  * hMuonEa;
		TH1D* hCharmEn, * hCharmEa;

		TH1D* hTotalMn, * hTotalMa;
		TH1D* hPionMn,  * hPionMa;
		TH1D* hKaonMn,  * hKaonMa;
		TH1D* hKaon0Mn, * hKaon0Ma;
		TH1D* hMuonMn,  * hMuonMa;
		TH1D* hCharmMn, * hCharmMa;

		TH1D* hTotalTn, * hTotalTa;
		TH1D* hPionTn,  * hPionTa;
		TH1D* h2PionTn, * h2PionTa;
		TH1D* hCharmTn, * hCharmTa;
		TH1D* hTauETn,  * hTauETa;
		TH1D* hTauMTn,  * hTauMTa;

		double M_Sterile_prev;

		const double M_Electron;
		const double M_Muon;
		const double M_Tau;
		const double M_Pion;
		const double M_Pion0;
		const double M_Kaon;
		const double M_Kaon0;
		const double M_Charm;
};

#endif
