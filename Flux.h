/*
 * Flux class, container of various components
 * 
 * Author: Tommaso Boschi
 */

#ifndef flux_H
#define flux_H

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

class Flux
{
	public:
		Flux();
		Flux(double Energy);
		Flux(double MuonPion, double MuonKaon, double ElectronPion, double ElectronKaon, double ElectronKaon3, double MuonKaonOther);

		double GetMuonPion();
		double GetMuonKaon();
		double GetElectronPion();
		double GetElectronKaon();
		double GetElectronKaon3();
		double GetMuonKaonOther();
		double GetTotalFlux();

		void SetMuonPion(double X);
		void SetMuonKaon(double X);
		void SetElectronPion(double X);
		void SetElectronKaon(double X);
		void SetElectronKaon3(double X);
		void SetMuonKaonOther(double X);
		void SetTotalFlux(double X);

	private:
		double fMuonPion;
		double fMuonKaon;
		double fElectronPion;
		double fElectronKaon;
		double fElectronKaon3;
		double fMuonKaonOther;
}

#endif
