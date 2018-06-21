/*
 * Particle class, it is a wrapper of a TLorentzVector (fourmomentum) and TVector3 (position)
 * plus some extra details of the particle such as PDG ID code and charge, if it is howering etc..
 *
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

//#include "Tools.h"

class Particle
{
	public:
		Particle();
		Particle(int PdgCode, double Px, double Py, double Pz, double E, double X, double Y, double Z);
		Particle(int PdgCode, TLorentzVector *Vector, TVector3 *Position);
		Particle(int PdgCode, TLorentzVector &Vector, TVector3 &Position);
		Particle(const Particle &P);


		int Pdg() const;
		int Charge();
		double LifeTime();
		double LabLifeTime();
		double LabSpace();
		bool IsShower() const;

		TLorentzVector FourVector() const;
		double Mass();
		double Energy();
		double EnergyKin();
		double Momentum();
		double Transverse();
		double MomentumX();
		double MomentumY();
		double MomentumZ();
		double Theta();
		double Phi();
		double Beta();
		double Gamma();
		TVector3 Direction();

		TVector3 Position() const;
		double X();
		double Y();
		double Z();
		double TrackIn() const;
		double TrackOut() const;
		double TrackTot() const;
		
		void SetPdg(int X);
		void SetFourVector(TLorentzVector &V);
		void SetFourVector(double Px, double Py, double Pz, double E);
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
		TLorentzVector *ParticleVec;
		TVector3 *ParticlePos;
		int iPdg;
		double dTrackIn, dTrackOut;
		bool bShower;
};

#endif
