/*
 * Flux class, container of various components
 * 
 * Author: Tommaso Boschi
 */

#ifndef flux_H
#define flux_H

#include <iostream>
#include <fstream>

//ROOT include
#include "TH1.h"
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"

class Flux
{
	public:
		Flux(double Energy,
		     double MuonPion = 0,
		     double MuonKaon = 0,
		     double ElectronPion = 0,
		     double ElectronKaon = 0,
		     double ElectronKaon3 = 0,
		     double MuonKaonOther = 0);

		void SetAll(double Energy,
			    double MuonPion = 0,
			    double MuonKaon = 0,
			    double ElectronPion = 0,
			    double ElectronKaon = 0,
			    double ElectronKaon3 = 0,
			    double MuonKaonOther = 0);

		void SetEnergy(double X);
		void SetMuonPion(double X);
		void SetMuonKaon(double X);
		void SetElectronPion(double X);
		void SetElectronKaon(double X);
		void SetElectronKaon3(double X);
		void SetMuonKaonOther(double X);
		void SetTotalFlux(double X);

		double GetEnergy();
		double GetMuonPion();
		double GetMuonKaon();
		double GetElectronPion();
		double GetElectronKaon();
		double GetElectronKaon3();
		double GetMuonKaonOther();
		double GetTotalFlux();

	private:
		double fEnergy;
		double fMuonPion;
		double fMuonKaon;
		double fElectronPion;
		double fElectronKaon;
		double fElectronKaon3;
		double fMuonKaonOther;
};

#endif
