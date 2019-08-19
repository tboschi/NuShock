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

#include "tools.h"
#include "DecayRates.h"
#include "Production.h"
#include "PhaseSpace.h"

class Particle;
class Neutrino : public Particle
{
	public:
		enum Neut
		{
			Left         = 0,
			Right        = 1,
			Polarised    = 0,
			Unpolarised  = 2,
			Dirac	     = 0,
			Majorana     = 4,
			//Fermion      = 0,
			Antiparticle = 8
		};

		Neutrino(double Mass = 0, unsigned int Options = 2);
		Neutrino(const Neutrino &N);
		~Neutrino();
		Neutrino & operator=(const Neutrino & N);

		void SetParent(Amplitude *Object);
		void SetUnitary(Amplitude *Object);

		double DecayThreshold();
		double DecayThreshold(std::string Name);
		double DecayThreshold(Amplitude::Channel Name);
		bool IsDecayAllowed();
		bool IsDecayAllowed(std::string Name);
		bool IsDecayAllowed(Amplitude::Channel Name);
		void DecayChannels(std::vector<std::string> &vChan);
		double DecayTotal();
		double DecayWidth();
		double DecayWidth(std::string Name);
		double DecayWidth(Amplitude::Channel name);
		double DecayBranch();
		double DecayBranch(std::string Name);
		double DecayBranch(Amplitude::Channel name);

		double ProductionThreshold();
		double ProductionThreshold(std::string Name);
		double ProductionThreshold(Amplitude::Channel Name);
		bool IsProductionAllowed();
		bool IsProductionAllowed(std::string Name);
		bool IsProductionAllowed(Amplitude::Channel Name);
		void ProductionChannels(std::vector<std::string> &vChan);
		double ProductionWidth();
		double ProductionWidth(std::string Name);
		double ProductionWidth(Amplitude::Channel name);
		double ProductionScale();
		double ProductionScale(std::string Name);
		double ProductionScale(Amplitude::Channel name);
		
		std::vector<Particle> DecayPS();
		std::vector<Particle> DecayPS(Amplitude::Channel Name);
		std::vector<Particle> ProductionPS(TLorentzVector &Vec);
		std::vector<Particle> ProductionPS(Amplitude::Channel Name, TLorentzVector &Vec);
		
		void SetDecayChannel(std::string Name);
		void SetProductionChannel(std::string Name);
		Amplitude::Channel DecayChannel();
		Amplitude::Channel ProductionChannel();
		std::string DecayChannelName();
		std::string ProductionChannelName();

		//void SetMass(double Mass = 0.0);
		void SetMixings(double *Mixings);
		void SetMixings(double Ue, double Um, double Ut);
		//void SetEnergy(double Energy);
		//void SetEnergyKin(double Energy);
		void SetHelicity(unsigned int Options);
		void SetFermion(unsigned int Options);
		void SetParticle(unsigned int Options);
		
		//double Mass();
		double* Mixings();
		double Ue(int E = 1.0) const;
		double Um(int E = 1.0) const;
		double Ut(int E = 1.0) const;
		//double Energy();
		//double EnergyKin();
		int Helicity();
		bool GetFermion();
		bool GetParticle();

		bool IsDirac();
		bool IsMajorana();
		bool IsParticle();
		bool IsAntiparticle();

	private:
		DecayRates *theDecayRates;
		Production *theProduction, *theProdLightN;
		PhaseSpace *thePhaseSpace;

		Amplitude::Channel chDecay, chProduction;

		//double fMass;
		//double fEnergy;

		double *fMixings;
		bool bParticle;
		bool bFermion;
		int  iHel;
};

#endif
