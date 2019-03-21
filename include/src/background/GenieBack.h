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
#include "gst.h"

#include "tools.h"
#include "detector.h"

class GenieBack
{
	public:
		GenieBack(std::string backFile);
		~GenieBack();
		void InitInTree(std::string backFile);
		void InitOutTree();
		void LoadTree(const std::vector<Particle> &vPart);
		TTree *FindBackground(Tracker *theTrack, std::map<int, int> &process, int save);
		bool MisIdentify(std::vector<Particle> &vPart, Tracker *theTrack);
		bool Identify(const std::vector<Particle> &vPart, const std::map<int, int> mProc);

	private:
		gst *genie;
		TTree *data;
		TFile *inBack;

		int checkPt;

		int ID, np;

		double *energy,
		       *moment,
                       *transv,
                       *theta,
                       *phi,
                       *mass,
                       *lenOut,
                       *lenIn;
};

#endif
