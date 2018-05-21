/*
 * Particle class, homemade
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
#include "TDatabasePDG.h"
#include "TParticlePDG.h"

class Particle
{
	public:
		Particle(int PdgCode, TLorentzVector &Vector, double PosX, double PosY, double PosZ);
		Particle(int PdgCode, TLorentzVector &Vector, TVector3 &Position);
		Particle(const Particle &P);

		int Pdg() const;
		double Charge() const;
		bool IsShower() const;
		double Tau() const;	//lifetime at rest
		double LabTau() const;		//lifetime in lab
		double LabSpace() const;		//space travelled in lab before decay
		TLorentzVector GetP4() const;
		TVector3 Direction() const;
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
		double TrackIn() const;
		double TrackOut() const;
		double TrackTot() const;

		void SetCharge(double Q);
		void SetP4(TLorentzVector &V);
		void SetP4(double Px, double Py, double Pz, double E);
		void SetPdg(int X);
		void SetCharge();
		void SetTau();
		void SetEnergy(double dE);
		void SetEnergyKin(double dE);
		void SetMomentum(double dP);
		void SetMass(double dM);
		void SetRho(double dR);
		void SetTheta(double Ang);
		void SetPhi(double Ang);
		void SetPosition(TVector3 &V);
		void SetPosition(double X, double Y, double Z);
		void SetX(double X);
		void SetY(double X);
		void SetZ(double X);
		void SetTrackIn(double X);
		void SetTrackOut(double X);
		void SetShower(bool X);

	private:
		TLorentzVector P4;
		TVector3 Pos;
		int iPdg, iCharge;
		double dTrackIn, dTrackOut, dTau;
		bool bShower;
};

#endif
