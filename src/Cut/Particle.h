/*
 * Cut analysis
 * Author: Tommaso Boschi
 */

#ifndef PARTICLE_H
#define PARTICLE_H

#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <string>

#include "TLorentzVector.h"
#include "TVector3.h"

#include "EVGCore/EventRecord.h"
#include "GHEP/GHepParticle.h"
#include "Ntuple/NtpMCFormat.h"
#include "Ntuple/NtpMCTreeHeader.h"
#include "Ntuple/NtpMCEventRecord.h"
#include "Messenger/Messenger.h"
#include "PDG/PDGCodes.h"
#include "Utils/CmdLnArgParser.h"

class Particle
{
	public:
		Particle(genie::GHepParticle *Candidate, double PosX, double PosY, double PosZ);
		Particle(int PdgCode, TLorentzVector *Vector, double PosX, double PosY, double PosZ);

		int Pdg();
		TLorentzVector GetP4();
		TVector3 Position();
		double X();
		double Y();
		double Z();
		double M();
		double E();
		double Ekin();
		double Px();
		double Py();
		double Pz();
		double P();
		double Theta();
		double Phi();

		void SetP4(TLorentzVector &V);
		void SetPdg(int PdgCode);
		void SetE(double E);
		void SetMass(double M);
		void SetTheta(double Ang);
		void SetPhi(double Ang);
		void SetPosition(TVector3 &V);
		void SetPosition(double X, double Y, double Z);
		void SetX(double X);
		void SetY(double X);
		void SetZ(double X);

	private:
		TLorentzVector P4;
		TVector3 Pos;
		int iPdg;
};

#endif
