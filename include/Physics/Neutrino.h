/*
 * Neutrino class
 * Author: Tommaso Boschi
 */

#ifndef NEUTRINO_H
#define NEUTRINO_H

#include <iostream>
#include <fstream>
#include <vector>
#include <cstring>

//ROOT include
#include "TLorentzVector.h"
#include "TGenPhaseSpace.h"

#include "Tools/Tools.h"
#include "Physics/ThreeBody.h"
#include "Physics/DecayRates"

enum Neut
{
	Left         = 0;
	Right        = 1;
	Polarised    = 0;
	Unpolarised  = 2;
	Dirac	     = 0;
	Majorana     = 4;
	Particle     = 0;
	Antiparticle = 8;
};

class Neutrino:
{
	public:
		enum Neut
		{
			Particle     = 0;
			Antiparticle = 1;
			Left         = 2;
			Right        = 4;
		};

		Neutrino(double Mass, unsigned int Options = 2);
		~Neutrino();

		double Decay(Channel name);
		double Branch(Channel name);

		void SetMass(double Mass = 0.0);
		void SetMixings(double Ue, double Um, double Ut);
		void SetEnergy(double Energy);

		double GetMass();
		double *GetMixings();


	private:
		DecayRates *TheDecay;

		double fMass;
		double *fMixings;
		bool bParticle;
		bool bFermion;
		int  iHelicity;

};

#endif
