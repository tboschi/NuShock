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

#include "Tools.h"
#include "Physics/DecayRates.h"
#include "Physics/Production.h"
#include "Physics/PhaseSpace.h"

class Neutrino:
{
	public:
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

		Neutrino(double Mass, unsigned int Options = 2);
		~Neutrino();

		double DecayWidth(Channel name);
		double DecayBranch(Channel name);
		double ProductionWidth(Channel name);
		double ProductionScale(Channel name);
		double BranchWidth(Channel name);
		std::vector<TLorentzVector> PhaseSpace(Channel Name)

		void SetMass(double Mass = 0.0);
		void SetMixings(double Ue, double Um, double Ut);
		void SetEnergy(double Energy);
		void SetEnergyKin(double Energy)
		void SetHelicity(unsigned int Options)		//Left for particle is -1
		void SetFermion(unsigned int Options)
		void SetParticle(unsigned int Options)
		
		double GetMass()
		double* GetMixings()
		double GetEnergy()
		double GetEnergyKin()
		int GetHelicity()

		bool IsDirac()
		bool IsMajorana()
		bool IsParticle()
		bool IsAntiparticle()

	private:
		DecayRates *TheDecayRates;
		Production *TheProduction;
		PhaseSpace *ThePhaseSpace;

		double fMass;
		double *fMixings;
		bool bParticle;
		bool bFermion;
		int  iHelicity;

};

#endif
