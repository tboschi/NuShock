/*
 * Class to drive neutrinos to the ND
 * Modifies published fluxes to extimate sterile neutrino flux
 *
 * Author: Tommaso Boschi
 */

#ifndef DRIVER_H
#define DRIVER_H

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

class Driver
{
	public:
		Driver(std::string ConfigFlux, bool FHC = 1);
		~Driver();

		void CloneCopy(TH1D*& T, TObject* X);
		bool MakeFlux(Neutrino *N);
		void MakeElecComponent(Flux *fxFlux, Neutrino *N);
		void MakeMuonComponent(Flux *fxFlux, Neutrino *N);
		void MakeTauComponent(Flux *fxFlux, Neutrino *N);

		double Intensity(Neutrino *N);
		double InterpolateIntensity(TH1D* Hist, double Energy);

		double Range();
		double RangeWidth();
		double Range(double &Start, double &End);
		double RangeWidth(double &Start, double &End);
		double RangeStart();
		double RangeEnd();
		int BinNumber();

		void Scale(double X);

		bool IsChanged(Neutrino* N);

		double Modify(double &xdir, double &ydir, double M_Sterile);

	private:
		TFile *SourceFile;
		TFile *KineFile;
		bool Kine, Mod;

		//Get fluxes from file
		Flux *fxNuElectron, *fxNuMuon, *fxNuTau;
		Flux *fxHeavyElectron, *fxHeavyMuon, *fxHeavyTau;

		double Mass_prev;
		int Helicity_prev;
		bool Particle_prev;

		std::vector<double> vMdir, vXdir, vYdir;
};

#endif
