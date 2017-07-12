/*
 * Background analysis
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
		Particle(int PdgCode, double Charge, TLorentzVector &Vector, double PosX, double PosY, double PosZ);
		Particle(int PdgCode, double Charge, TLorentzVector &Vector, TVector3 &Position);
		Particle(const Particle &P);

		int Pdg() const;
		int Charge() const;
		TLorentzVector GetP4() const;
		TVector3 Position() const;
		double X() const;
		double Y() const;
		double Z() const;
		double M() const;
		double E() const;
		double Ekin() const;
		double Px() const;
		double Py() const;
		double Pz() const;
		double P() const;
		double Pt() const;
		double Theta() const;
		double Phi() const;

		void SetP4(TLorentzVector &V);
		void SetP4(double Px, double Py, double Pz, double E);
		void SetPdg(int PdgCode);
		void SetCharge(int Charge);
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
		int iPdg, iCharge;
};

#endif
