//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Mar 20 21:53:51 2019 by ROOT version 5.34/36
// from TTree gst/GENIE Summary Event Tree
// found on file: numu.gst.root
//////////////////////////////////////////////////////////

#ifndef gst_h
#define gst_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <iostream>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class gst
{
	public :
		TTree         *fChain;   //!pointer to the analyzed TTree or TChain

		// Declaration of leaf types
		int           neu;
		bool          qel;
		bool          mec;
		bool          res;
		bool          dis;
		bool          coh;
		bool          cc;
		bool          nc;
		double        wght;
		double        Ev;
		double        pxv;
		double        pyv;
		double        pzv;
		double        En;
		double        pxn;
		double        pyn;
		double        pzn;
		double        El;
		double        pxl;
		double        pyl;
		double        pzl;
		int           nf;
		int           pdgf[47];   //[nf]
		double        Ef[47];   //[nf]
		double        pxf[47];   //[nf]
		double        pyf[47];   //[nf]
		double        pzf[47];   //[nf]
		double        pf[47];   //[nf]

		// List of branches
		TBranch        *b_neu;   //!
		TBranch        *b_qel;   //!
		TBranch        *b_mec;   //!

		TBranch        *b_res;   //!
		TBranch        *b_dis;   //!
		TBranch        *b_coh;   //!
		TBranch        *b_cc;   //!
		TBranch        *b_nc;   //!
		TBranch        *b_wght;   //!
		TBranch        *b_Ev;   //!
		TBranch        *b_pxv;   //!
		TBranch        *b_pyv;   //!
		TBranch        *b_pzv;   //!
		TBranch        *b_En;   //!
		TBranch        *b_pxn;   //!
		TBranch        *b_pyn;   //!
		TBranch        *b_pzn;   //!
		TBranch        *b_El;   //!
		TBranch        *b_pxl;   //!
		TBranch        *b_pyl;   //!
		TBranch        *b_pzl;   //!
		TBranch        *b_nf;   //!
		TBranch        *b_pdgf;   //!
		TBranch        *b_Ef;   //!
		TBranch        *b_pxf;   //!
		TBranch        *b_pyf;   //!
		TBranch        *b_pzf;   //!
		TBranch        *b_pf;   //!

		gst(TTree *tree=0);
		~gst();
		int    GetEntry(long entry);
		int    GetEntries();
		void     Init(TTree *tree);
};

#endif
