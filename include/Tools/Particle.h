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
		Particle(int PdgCode, TLorentzVector *Vector);
		Particle(int PdgCode, TLorentzVector &Vector, TVector3 &Position);
		Particle(int PdgCode, TLorentzVector &Vector);
		Particle(const Particle &P);
		~Particle();
		void Init(double Px = 0.0, double Py = 0.0, double Pz = 0.0, double E = 0.0,
			  double X = 0.0, double Y = 0.0, double Z = 0.0);

		int Pdg() ;
		int Charge();
		double LifeTime();
		double LabLifeTime();
		double LabSpace();
		bool IsShower();

		//TLorentzVector FourVectorPtr();
		TLorentzVector FourVector() ;
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

		TVector3 Position();
		double X();
		double Y();
		double Z();
		double Dist();
		double TrackIn();
		double TrackBack();
		double TrackOut();
		double TrackTot();
		
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
		void SetTrackBack(double X);
		void SetTrackOut(double X);
		void SetShower(bool X);

	private:
		TLorentzVector ParticleVec;
		TVector3 ParticlePos;
		int iPdg;
		double dTrackIn, dTrackBack, dTrackOut;
		bool bShower;
};

#endif
