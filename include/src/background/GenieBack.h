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
		//
		double e_a;
		double p_a;
		double t_a;
		double thea;
		double phia;
		//
		double E_B;
		double P_B;
		double T_B;
		double TheB;
		double PhiB;
		double M_B;
		double In_B;
		double Out_B;
		//
		double e_b;
		double p_b;
		double t_b;
		double theb;
		double phib;
		//
		double Angle;
		//
		double angle;
		//
		double E_0;
		double P_0;
		double T_0;
		double The0;
		double Phi0;
		double M_0;
		//
		double e_0;
		double p_0;
		double t_0;
		double the0;
		double phi0;
		double m_0;
};

#endif
