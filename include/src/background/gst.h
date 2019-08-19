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
		int           iev;
		int           neu;
		int           fspl;
		int           tgt;
		int           Z;
		int           A;
		int           hitnuc;
		int           hitqrk;
		int           resid;
		bool          sea;
		bool          qel;
		bool          mec;
		bool          res;
		bool          dis;
		bool          coh;
		bool          dfr;
		bool          imd;
		bool          imdanh;
		bool          singlek;
		bool          nuel;
		bool          em;
		bool          cc;
		bool          nc;
		bool          charm;
		int           neut_code;
		int           nuance_code;
		double        wght;
		double        xs;
		double        ys;
		double        ts;
		double        Q2s;
		double        Ws;
		double        x;
		double        y;
		double        t;
		double        Q2;
		double        W;
		double        EvRF;
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
		double        pl;
		double        cthl;
		int           nfp;
		int           nfn;
		int           nfpip;
		int           nfpim;
		int           nfpi0;
		int           nfkp;
		int           nfkm;
		int           nfk0;
		int           nfem;
		int           nfother;
		int           nip;
		int           nin;
		int           nipip;
		int           nipim;
		int           nipi0;
		int           nikp;
		int           nikm;
		int           nik0;
		int           niem;
		int           niother;
		int           ni;
		int           pdgi[14];   //[ni]
		int           resc[14];   //[ni]
		double        Ei[14];   //[ni]
		double        pxi[14];   //[ni]
		double        pyi[14];   //[ni]
		double        pzi[14];   //[ni]
		int           nf;
		int           pdgf[47];   //[nf]
		double        Ef[47];   //[nf]
		double        pxf[47];   //[nf]
		double        pyf[47];   //[nf]
		double        pzf[47];   //[nf]
		double        pf[47];   //[nf]
		double        cthf[47];   //[nf]
		double        vtxx;
		double        vtxy;
		double        vtxz;
		double        vtxt;
		double        sumKEf;
		double        calresp0;

		// List of branches
		TBranch        *b_iev;   //!
		TBranch        *b_neu;   //!
		TBranch        *b_fspl;   //!
		TBranch        *b_tgt;   //!
		TBranch        *b_Z;   //!
		TBranch        *b_A;   //!
		TBranch        *b_hitnuc;   //!
		TBranch        *b_hitqrk;   //!
		TBranch        *b_resid;   //!
		TBranch        *b_sea;   //!
		TBranch        *b_qel;   //!
		TBranch        *b_mec;   //!
		TBranch        *b_res;   //!
		TBranch        *b_dis;   //!
		TBranch        *b_coh;   //!
		TBranch        *b_dfr;   //!
		TBranch        *b_imd;   //!
		TBranch        *b_imdanh;   //!
		TBranch        *b_singlek;   //!
		TBranch        *b_nuel;   //!
		TBranch        *b_em;   //!
		TBranch        *b_cc;   //!
		TBranch        *b_nc;   //!
		TBranch        *b_charm;   //!
		TBranch        *b_neut_code;   //!
		TBranch        *b_nuance_code;   //!
		TBranch        *b_wght;   //!
		TBranch        *b_xs;   //!
		TBranch        *b_ys;   //!
		TBranch        *b_ts;   //!
		TBranch        *b_Q2s;   //!
		TBranch        *b_Ws;   //!
		TBranch        *b_x;   //!
		TBranch        *b_y;   //!
		TBranch        *b_t;   //!
		TBranch        *b_Q2;   //!
		TBranch        *b_W;   //!
		TBranch        *b_EvRF;   //!
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
		TBranch        *b_pl;   //!
		TBranch        *b_cthl;   //!
		TBranch        *b_nfp;   //!
		TBranch        *b_nfn;   //!
		TBranch        *b_nfpip;   //!
		TBranch        *b_nfpim;   //!
		TBranch        *b_nfpi0;   //!
		TBranch        *b_nfkp;   //!
		TBranch        *b_nfkm;   //!
		TBranch        *b_nfk0;   //!
		TBranch        *b_nfem;   //!
		TBranch        *b_nfother;   //!
		TBranch        *b_nip;   //!
		TBranch        *b_nin;   //!
		TBranch        *b_nipip;   //!
		TBranch        *b_nipim;   //!
		TBranch        *b_nipi0;   //!
		TBranch        *b_nikp;   //!
		TBranch        *b_nikm;   //!
		TBranch        *b_nik0;   //!
		TBranch        *b_niem;   //!
		TBranch        *b_niother;   //!
		TBranch        *b_ni;   //!
		TBranch        *b_pdgi;   //!
		TBranch        *b_resc;   //!
		TBranch        *b_Ei;   //!
		TBranch        *b_pxi;   //!
		TBranch        *b_pyi;   //!
		TBranch        *b_pzi;   //!
		TBranch        *b_nf;   //!
		TBranch        *b_pdgf;   //!
		TBranch        *b_Ef;   //!
		TBranch        *b_pxf;   //!
		TBranch        *b_pyf;   //!
		TBranch        *b_pzf;   //!
		TBranch        *b_pf;   //!
		TBranch        *b_cthf;   //!
		TBranch        *b_vtxx;   //!
		TBranch        *b_vtxy;   //!
		TBranch        *b_vtxz;   //!
		TBranch        *b_vtxt;   //!
		TBranch        *b_sumKEf;   //!
		TBranch        *b_calresp0;   //!

		gst(TTree *tree=0);
		~gst();
		int    GetEntry(long entry);
		int    GetEntries();
		void     Init(TTree *tree);
};

#endif
