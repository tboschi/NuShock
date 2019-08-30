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

#include "tools.h"

#include "TTree.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TRandom3.h"

class Efficiency
{
	public:
		enum Cut
		{
			_E_A,
			_P_A,
			_T_A,
			_TheA,
			_PhiA,
			_In_A,
			_Out_A,
			_E_B,
			_P_B,
			_T_B,
			_TheB,
			_PhiB,
			_In_B,
			_Out_B,
			_Angle,
			_E_0,
			_P_0,
			_T_0,
			_The0,
			_Phi0,
			_M_0,
			     
			_CosAB,
			_aCosAB,
			_CircAB,
			_atAB0,
			_T_AB,
			_T_AB0,
			_TTA,
			_EAB,
			_E0Ang
		};

		Efficiency(std::string mcFile = "");
		~Efficiency();
		void MapCuts();
		void LoadFile(std::string mcFile);
		std::vector<std::string> AvailableCuts(std::string name = "");
		void Reset();
		void LoadTree(TTree *mcData);
		void GetCutLimits(std::string name, double &cLo, double &cUp);
		void FindCut(std::string cutName, double &cLo, double &cUp, double CL, double mass = -1);
		void FindRange(TH1D *hist, int &sL, int &sR, int s0, double CL);
		int GetEntries();
		void GetEntry(int i);
		TH1D* LoadCutSpectrum(std::string name);
		void LoadCut(std::string cutFile);
		void SetCut(std::string name, double lower = 0, double upper = 0);
		void ApplyCut(double mass = -1.0);
		int EntriesLeft();
		double EventsLeft();
		double ReductionFactor();
		bool PassCut(std::string name = "");
		TH1D *GetAll();
		TH1D *GetCut();
		int FindFirstBin(TH1D* hist, double thr = 0, int start = -1, int end = -1);
		void MakeFunction();
		TH2D* CompleteFunction();

	private:
		bool funcSet;
		TTree *Data;
		TFile *inFile;

		TH1D *hAll, *hCut; 
		TH2D* hhFunc;
		//std::map<double, TH1D*> mAll, mCut;

		std::map<std::string, Cut> mCut;
		std::map<std::string, Cut>::iterator ic;
		std::map<std::string, double*> mRef;
		std::map<std::string, double*>::iterator im;
		std::map<std::string, double> mLower, mUpper;

		double *Hist;
		double  True, W;

		int PdgA;
		double E_A;
		double P_A;
		double T_A;
		double TheA;
		double PhiA;
		double In_A;
		double Out_A;
		int PdgB;
		double E_B;
		double P_B;
		double T_B;
		double TheB;
		double PhiB;
		double In_B;
		double Out_B;
		double Angle;
		double E_0;
		double P_0;
		double T_0;
		double The0;
		double Phi0;
		double M_0;
		//special
		double CosAB;
		double aCosAB;
		double CircAB;
		double atAB0;
		double T_AB;
		double T_AB0;
		double TTA;
		double EAB;
		double E0Ang;

		// List of branches
		TBranch        *b_fPdgA;   //!
		TBranch        *b_fEnergyA;   //!
		TBranch        *b_fMomentA;   //!
		TBranch        *b_fTransvA;   //!
		TBranch        *b_fThetaA;   //!
		TBranch        *b_fPhiA;   //!
		TBranch        *b_fLengthA;   //!
		TBranch        *b_fLengthoA;   //!
		TBranch        *b_fPdgB;   //!
		TBranch        *b_fEnergyB;   //!
		TBranch        *b_fMomentB;   //!
		TBranch        *b_fTransvB;   //!
		TBranch        *b_fThetaB;   //!
		TBranch        *b_fPhiB;   //!
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
