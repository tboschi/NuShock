/*
 * Background analysis
 * Author: Tommaso Boschi
 */

#ifndef GenieBack_H
#define GenieBack_H

#include <iostream>
#include <vector>
#include <map>
#include <string>

#include "TLorentzVector.h"
#include "TVector3.h"
#include "TMath.h"
#include "TTree.h"
#include "TFile.h"
#include "TH1D.h"

#include "tools.h"
#include "detector.h"
#include "src/background/gst.h"

class GenieBack
{
	public:
		GenieBack(std::string backFile, std::string outFile,
				bool verb = false);
		~GenieBack();
		int Charge(int pdg);
		void InitInTree(std::string backFile);
		void InitOutTree(std::string outFile);
		void LoadTree(const std::vector<Particle> &particle,
			      const std::vector<Particle> &original);
		TTree *FindBackground(Tracker *theTrack,
				      const std::map<int, int> &process, int save = 1);
		bool MisIdentify(std::vector<Particle> &vPart, Tracker *theTrack);
		bool Identify(const std::vector<Particle> &particle,
			      const std::map<int, int> &process);
		bool IsHadron(const Particle &p);

	private:
		gst *genie;
		TTree *data;
		TFile *inBack, *outBack;
		TH1D* hist;

		bool chargeID, kVerbose;

		int checkPt;

		int ID, np, nr;

		double *energy,
		       *moment,
                       *transv,
                       *theta,
                       *phi,
                       *mass,
                       *lenOut,
                       *lenIn,
		       *r_energy,
		       *r_moment,
		       *r_transv,
		       *r_theta, 
		       *r_phi,   
		       *r_mass,  
		       *r_angle; 

		int PdgA, PdgB;
		double True;
		double E_A;
		double P_A;
		double T_A;
		double TheA;
		double PhiA;
		double M_A;
		double In_A;
		double Out_A;
		double e_A;
		double p_A;
		double t_A;
		double theA;
		double phiA;
		double E_B;
		double P_B;
		double T_B;
		double TheB;
		double PhiB;
		double M_B;
		double In_B;
		double Out_B;
		double e_B;
		double p_B;
		double t_B;
		double theB;
		double phiB;
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
