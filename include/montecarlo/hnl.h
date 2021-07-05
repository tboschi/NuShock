#ifndef HNL_H
#define HNL_H

#include <iostream>
#include "detector/Tracker.h"

#include "TFile.h"
#include "TTree.h"

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class hnl
{
	public: 
		double True;	// true neutrino energy
		double Vert[3]; // vertex position
		double W;	// event weight
		bool ChID;	// if event was correctly charge id'ed

		// particle A info
		int PdgA;
		int ChA;
		double E_A;
		double P_A;
		double T_A;
		double TheA;
		double PhiA;
		double M_A;
		double In_A;
		double Out_A;

		// particle B info
		int PdgB;
		int ChB;
		double E_B;
		double P_B;
		double T_B;
		double TheB;
		double PhiB;
		double M_B;
		double In_B;
		double Out_B;

		// reconstruncted A+B
		double Angle;
		double E_0;
		double P_0;
		double T_0;
		double The0;
		double Phi0;
		double M_0;

	private:
		TBranch * b_True;
		TBranch * b_Vert;
		TBranch * b_W;
		TBranch * b_ChID;

		// particle A info
		TBranch * b_PdgA;
		TBranch * b_ChA;
		TBranch * b_E_A;
		TBranch * b_P_A;
		TBranch * b_T_A;
		TBranch * b_TheA;
		TBranch * b_PhiA;
		TBranch * b_M_A;
		TBranch * b_In_A;
		TBranch * b_Out_A;

		// particle B info
		TBranch * b_PdgB;
		TBranch * b_ChB;
		TBranch * b_E_B;
		TBranch * b_P_B;
		TBranch * b_T_B;
		TBranch * b_TheB;
		TBranch * b_PhiB;
		TBranch * b_M_B;
		TBranch * b_In_B;
		TBranch * b_Out_B;

		// reconstruncted A+B
		TBranch * b_Angle;
		TBranch * b_E_0;
		TBranch * b_P_0;
		TBranch * b_T_0;
		TBranch * b_The0;
		TBranch * b_Phi0;
		TBranch * b_M_0;


	public: 
		hnl(std::string name = "hnl");	// write mode
		hnl(TTree *tree);		// read mode
		~hnl();

		size_t GetEntries();
		int GetEntry(size_t entry);
		void Fill();
		// special fill
		void Fill(double weight, bool chid,
			const Tracker::Event &e0,
			const Tracker::Event &e1,
			const Tracker::Event &e2);
		void Write();
		TTree* Chain();

		double Events();

	private:
		void Init(TTree *tree);
		void New(std::string name);

		TTree* fChain;
		bool own;
};

#endif
