/*
 * Particle class, it is a wrapper of a TLorentzVector (fourmomentum) and TVector3 (position)
 * plus some extra details of the particle such as PDG ID code and charge, if it is howering etc..
 *
 * Author: Tommaso Boschi
 */

#ifndef PARTICLE_H
#define PARTICLE_H

#include <iostream>

#include "TLorentzVector.h"

class Particle : public TLorentzVector
{
	public:
		Particle(int pdgCode = 0, double E = 0, double Px = 0, double Py = 0, double Pz = 0);
		Particle(int pdgCode, const TLorentzVector &fv);

		int Pdg() const;
		int RealPdg() const;
		int Q() const;
		int RealQ() const;
		double LifeTime() const;
		double LabLifeTime() const;
		double LabSpace() const;

		double EKin() const;
		double M() const;

		void SetPdg(int X);
		void ChargeMisID();
		void ChargeID();
		void ChargeID(int X);
		void SetFourVector(const TLorentzVector &vv);
		void SetFourVector(double E, double Px, double Py, double Pz);
		void SetE(double dE);
		void SetEKin(double dE);
		void SetP(double dP);
		void SetM(double dM);
		void SetRho(double dR);

		// return charge associated to pdg code
		static int Q(int pdg) {
			int sign = (pdg > 0) - (pdg < 0);
			switch (std::abs(pdg))
			{
				case 1:
				case 3:
				case 5:
				case 7:		//down quarks
					return -1 * sign;
				case 2:
				case 4:
				case 6:
				case 8:		//up quarks
					return  2 * sign;
				case 12:
				case 14:
				case 16:
				case 18:
				case 21:
				case 22:
				case 23:
				case 111:
				case 221:
				case 331:
				case 223:
				case 333:
				case 130:
				case 310:
				case 311:
				case 421:
				case 2112:
				case 2114:
				case 3122:
				case 3212:	//chargeless
					return 0;
				case 11:
				case 13:
				case 15:
				case 17:	//negative particles
					return -3 * sign;
				default:	//positive particles
					return  3 * sign;
					break;
			}
		}

	protected:
		int _pdg, _charge;

	friend std::ostream & operator<<(std::ostream &os, const Particle &p);
};

#endif
