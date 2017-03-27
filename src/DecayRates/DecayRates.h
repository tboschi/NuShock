/*
 * Flux creator for Monte Carlo using published fluxes
 *
 * Author: Tommaso Boschi
 */

#ifndef decayrates_H
#define decayrates_H

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

class Decay
{
	public:
		Decay();	//Decay rates calculator

		double Total(double M_Sterile, double U_e, double U_m, double U_t);
		double nEE(double M_Sterile, double U_e, double U_m, double U_t);
		double EPI(double M_Sterile, double U_e, double U_m, double U_t);
		double MUPI(double M_Sterile, double U_e, double U_m, double U_t);
		double nPI0(double M_Sterile, double U_e, double U_m, double U_t);
		double nGAMMA(double M_Sterile, double U_e, double U_m, double U_t);
		double nMUMU(double M_Sterile, double U_e, double U_m, double U_t);

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
