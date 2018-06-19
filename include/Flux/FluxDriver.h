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
#include "Physics.h"
#include "Flux/Flux.h"

class FluxDriver
{
	public:
		FluxDriver(std::string ConfigFlux);
		~FluxDriver();

		void CloneCopy(TH1D*& T, TObject* X);
		bool MakeFlux(Neutrino *N);
		void MakeElecComponent(Flux *fxFlux, Neutrino *N);
		void MakeMuonComponent(Flux *fxFlux, Neutrino *N);
		void MakeTauComponent(Flux *fxFlux, Neutrino *N);
		//void MakeStandardFlux();
		//double SampleEnergy();

		double Range();
		double RangeBin();
		double Range(double &Start, double &End);
		double RangeBin(double &Start, double &End);
		double RangeStart();
		double RangeEnd();
		int BinNumber();

		double Intensity(Neutrino *N);
		double InterpolateIntensity(TH1D* Hist, double Energy);

		void SetBaseline(double Baseline);
		void SetPOT(double POT);
		void SetArea(double Area);
		void ScaleAll(double X);

		bool IsChanged(Neutrino* N);

		double Modify(double &xdir, double &ydir, double M_Sterile);

	private:
		TFile *SourceFile;
		TFile *KineFile;
		bool Kine, Mod;

		//Get fluxes from file
		Flux *fxNuElectron, *fxNuMuon, *fxNuTau;
		Flux *fxHeavyElectron, *fxHeavyMuon, *fxHeavyTau;

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

		double Mass_prev;
		int Helicity_prev, Particle_prev;

		const double M_Electron;
		const double M_Muon;
		const double M_Tau;
		const double M_Pion;
		const double M_Pion0;
		const double M_Kaon;
		const double M_Kaon0;
		const double M_Charm;

		std::vector<double> vMdir, vXdir, vYdir;
};

#endif
