/*
 * Efficiency analysis
 * Author: Tommaso Boschi
 */

#ifndef EFFICIENCY_H
#define EFFICIENCY_H

#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <cstring>
#include <sstream>

#include "TTree.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TRandom3.h"

class Efficiency
{
	public:
		Efficiency(std::string InFile);
		~Efficiency();
		void LoadSpectra(double mass = -1.0);
		double EventsLeft();
		double ReductionFactor();
		void LoadTree(TTree *mcData);
		void LoadCut(std::string cutFile);
		void SetCut(std::string BN, double *Address, double Lo, double Up);
		void SetSpecial(int cutNumber, double Lo, double Up);
		bool PassCut();
		bool SpecialCut();
		TH1D *GetAll();
		TH1D *GetCut();
		//TH2D* MakeFunction();

	private:
		TTree *Data;
		TFile *inFile;

		TH1D *hAll, *hCut; 
		std::map<double, TH1D*> mAll, mCut;

		std::map<std::string, double*> mRef;
		std::map<std::string, double*>::iterator im;
		std::map<std::string, double> mCutLo, mCutUp;
		std::map<int, double> mSpecialLo, mSpecialUp;
		std::map<int, double>::iterator is;

		double *Hist;
		double  True, W;

		double E_A;
		double P_A;
		double T_A;
		double TheA;
		double PhiA;
		double M_A;
		double In_A;
		double Out_A;
		double E_B;
		double P_B;
		double T_B;
		double TheB;
		double PhiB;
		double M_B;
		double In_B;
		double Out_B;
		double Angle;
		double E_0;
		double P_0;
		double T_0;
		double The0;
		double Phi0;
		double M_0;

		// List of branches
		TBranch        *b_iID;   //!
		TBranch        *b_fEnergyA;   //!
		TBranch        *b_fMomentA;   //!
		TBranch        *b_fTransvA;   //!
		TBranch        *b_fThetaA;   //!
		TBranch        *b_fPhiA;   //!
		TBranch        *b_fMassA;   //!
		TBranch        *b_fLengthA;   //!
		TBranch        *b_fLengthoA;   //!
		TBranch        *b_fEnergyB;   //!
		TBranch        *b_fMomentB;   //!
		TBranch        *b_fTransvB;   //!
		TBranch        *b_fThetaB;   //!
		TBranch        *b_fPhiB;   //!
		TBranch        *b_fMassB;   //!
		TBranch        *b_fLengthB;   //!
		TBranch        *b_fLengthoB;   //!
		TBranch        *b_fAngle;   //!
		TBranch        *b_fEnergy0;   //!
		TBranch        *b_fMoment0;   //!
		TBranch        *b_fTransv0;   //!
		TBranch        *b_fTh0ta0;   //!
		TBranch        *b_fPhi0;   //!
		TBranch        *b_fMass0;   //!
};

#endif
