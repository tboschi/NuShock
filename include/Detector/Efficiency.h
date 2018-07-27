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
#include "TH1.h"
#include "TH2.h"
#include "TRandom3.h"

class Efficiency
{
	public:
		Efficiency(std::string InFile);

		void InitFunc();
		void LoopFile();
		void LoopTree();
		void InitTree();
		void LoadCut(std::string CutFile);
		void SetMap(std::string BN, double *Address, double Up, double Lo);
		void SetSpecial(int CutNumber, double Lo, double Up);
		void FillAll();
		void FillCut();
		bool PassCut();
		bool SpecialCut();
		void LoadFunction(double Mass);
		void MakeFunction();

		TH2D *GetFunction();
		TH1D *GetAll();
		TH1D *GetCut();

	private:
		TTree *Data, *Back;
		TFile *TreeFile, BkgFile;

		TH2D *hhFunc;
		TH1D *hCut, *hAll; 

		double dMass;
		std::map<std::string, double*> mRef;
		std::map<std::string, double*>::iterator im;
		std::map<std::string, double> mCutLo, mCutUp;
		std::map<int, double> mSpecialLo, mSpecialUp;
		std::map<int, double>::iterator is;
		std::vector<double> vMass;
		std::vector<std::string> vSim, vCut;

		double  True, W;

		double  E_A;
		double  P_A;
		double  T_A;
		double  TheA;
		double  PhiA;
		//Double_   M_A;	//not useful
		double  E_B;
		double  P_B;
		double  T_B;
		double  TheB;
		double  PhiB;
		//Double_   M_B;	//not useful
		double  Angle;
		double  E_0;
		double  P_0;
		double  T_0;
		double  The0;
		//Double_   Phi0;	//not useful
		double  M_0;
		
		TBranch  *b_fTrue;   //!
		TBranch  *b_fW;   //!
		TBranch  *b_fEnergyA;   //!
		TBranch  *b_fMomentA;   //!
		TBranch  *b_fTransvA;   //!
		TBranch  *b_fThetaA;   //!
		TBranch  *b_fPhiA;   //!
		//TBranch  *b_fMassA;   //not useful
		TBranch  *b_fEnergyB;   //!
		TBranch  *b_fMomentB;   //!
		TBranch  *b_fTransvB;   //!
		TBranch  *b_fThetaB;   //!
		TBranch  *b_fPhiB;   //!
		//TBranch  *b_fMassB;   //not useful
		TBranch  *b_fAngle;   //!
		TBranch  *b_fEnergy0;   //!
		TBranch  *b_fMoment0;   //!
		TBranch  *b_fTransv0;   //!
		TBranch  *b_fTheta0;   //!
		//TBranch  *b_fPhi0;   //not useful
		TBranch  *b_fMass0;   //!
};

#endif
