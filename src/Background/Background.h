/*
 * Background analysis
 * Author: Tommaso Boschi
 */

#ifndef BACKGROUND_H
#define BACKGROUND_H

#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <cstring>

#include "TLorentzVector.h"
#include "TVector3.h"
#include "TRandom3.h"
#include "TMath.h"
#include "TTree.h"
#include "TFile.h"
#include "TIterator.h"

#include "EVGCore/EventRecord.h"
#include "GHEP/GHepParticle.h"
#include "Ntuple/NtpMCFormat.h"
#include "Ntuple/NtpMCTreeHeader.h"
#include "Ntuple/NtpMCEventRecord.h"
#include "Messenger/Messenger.h"
#include "PDG/PDGCodes.h"
#include "Utils/CmdLnArgParser.h"

#include "Particle.h"
#include "Detector.h"
#include "Tools.h"

class Background
{
	public:
		//Background(std::string EventDB, std::string DetectorConfig, std::string OutFile, std::string Channel);
		Background(std::string EventDB, std::string DetectorConfig, std::string Channel);
		~Background();
		void InitTree();
		void InitMap();
		TTree *GetTree();
		void LoadTree();

		std::string GetChannel();
		void Loop(unsigned int Save);
		Particle* CreateParticle(genie::GHepParticle *Hep, double PosX, double PosY, double PosZ);

		int Count(std::string PartName, int N = 1);
		void ListCount();
		bool CountParticles();

		bool Identify();
		void IdentifyHadron(Particle *iP);
		bool IdentifynEE();
		bool IdentifynEMU();
		bool IdentifynPI0();
		bool IdentifyEPI();
		bool IdentifynMUMU();
		bool IdentifyMUPI();

		double GammaDecay();
		void Pi0Decay(Particle *Pi0, Particle *&PA, Particle *&PB);

	private:
		int NEvt;
		TRandom3 *GenMT;
		TTree *Genie, *Data;
		TFile *InFile, *OutFile;
		genie::NtpMCEventRecord *gEvRec;
		Detector *TheBox;
		std::string TheChan;

		Particle *ParticleA, *ParticleB;

		//Maps
		std::map<std::string, ChannelName> mapChan;
		std::map<std::string, int> pCount;
		std::vector<Particle*> vParticle;
		std::vector<Particle*>::iterator iP;

		unsigned int ID, Global, Hadron;
		double EnergyA, MomentA, TransvA, ThetaA, PhiA, MassA;
		double EnergyB, MomentB, TransvB, ThetaB, PhiB, MassB;
		double Energy0, Moment0, Transv0, Theta0, Phi0, Mass0;

		//Masses
		const double M_Neutrino;
		const double M_Photon;
		const double M_Electron;
		const double M_Muon;
		const double M_Pion;
		const double M_Pion0;
		const double M_Kaon;
		const double M_Kaon0;
};

#endif
