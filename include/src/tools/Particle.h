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
		Particle(int pdgCode = 0,  double Px = 0, double Py = 0, double Pz = 0,
			 double E = 0,     double X = 0,  double Y = 0,  double Z = 0);
		Particle(int pdgCode, TLorentzVector *Vector, TVector3 *Position = 0);
		Particle(int pdgCode, const TLorentzVector &Vector);
		Particle(int pdgCode, const TLorentzVector &Vector, const TVector3 &Position);
		Particle(const Particle &P);
		Particle& operator=(const Particle &rhs);
		~Particle();

		void Init(double Px = 0.0, double Py = 0.0, double Pz = 0.0, double E = 0.0,
			  double X = 0.0, double Y = 0.0, double Z = 0.0);

		int Pdg()  const;
		int Charge() const;
		double LifeTime() const;
		double LabLifeTime() const;
		double LabSpace() const;
		bool IsShower() const;

		TLorentzVector FourVector() const;
		double Mass() const;
		double Energy() const;
		double EnergyKin() const;
		double Momentum() const;
		double Transverse() const;
		double MomentumX() const;
		double MomentumY() const;
		double MomentumZ() const;
		double Theta() const;
		double Phi() const;
		double Beta() const;
		double Gamma() const;
		TVector3 Direction() const;

		TVector3 Position() const;
		double X() const;
		double Y() const;
		double Z() const;
		double Dist() const;
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
		TLorentzVector particleVec;
		TVector3 particlePos;
		int pdg;
		double trackIn, trackOut;
		bool kShower;
};

#endif
