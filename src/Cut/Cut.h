/*
 * Cut analysis
 * Author: Tommaso Boschi
 */

#ifndef CUT_H
#define CUT_H

#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <cstring>

#include <string>

#include "TLorentzVector.h"
#include "TVector3.h"
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TIterator.h"
#include "TH1.h"
#include "TGraph.h"

#include "EVGCore/EventRecord.h"
#include "GHEP/GHepParticle.h"
#include "Ntuple/NtpMCFormat.h"
#include "Ntuple/NtpMCTreeHeader.h"
#include "Ntuple/NtpMCEventRecord.h"
#include "Messenger/Messenger.h"
#include "PDG/PDGCodes.h"
#include "Utils/CmdLnArgParser.h"

class Cut
{
	public:

	private:
		//Masses
		const double M_Neutrino;
		const double M_Photon;
		const double M_Electron;
		const double M_Muon;
		const double M_Pion;
		const double M_Pion0;
		const double M_Kaon;
		const double M_Kaon0;

		//Maps
		std::map<std::string, ChannelName> mapCut;
		std::map<std::string, int> mapCount;
		std::vector<Particle*> vParticle;
		std::vector<Particle*>::iterator iP;

		TRandom3 *GenMT;

};

class Particle
{
	public:
		Particle(genie::gHepParticle *Candidate);
		Particle(int Pdg, TLorentzVector *P4);

		int Pdg();
		TLorentzVector GetP4();
		double E();
		double Px();
		double Py();
		double Pz();
		double P();
		double Theta();
		double Phi();

		void SetVector(TLorentzVector Vector);
		void SetPdg(int PdgCode);

	private:
		TLorentzVector P4;
		TVector3 Pos;
		int iPdg;
};

#endif
