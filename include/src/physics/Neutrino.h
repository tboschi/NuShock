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

		double DecayThreshold(std::string name = "");
		double DecayThreshold(Amplitude::Channel name);
		bool IsDecayAllowed(std::string name = "");
		bool IsDecayAllowed(Amplitude::Channel name);
		void DecayChannels(std::vector<std::string> &vChan);
		double DecayTotal();
		double DecayWidth(std::string name = "");
		double DecayWidth(Amplitude::Channel name);
		double DecayBranch(std::string name = "");
		double DecayBranch(Amplitude::Channel name);

		double ProductionThreshold(std::string name = "");
		double ProductionThreshold(Amplitude::Channel name);
		bool IsProductionAllowed(std::string name = "");
		bool IsProductionAllowed(Amplitude::Channel name);
		void ProductionChannels(std::vector<std::string> &vChan);
		double ProductionWidth(std::string name = "");
		double ProductionWidth(Amplitude::Channel name);
		double ProductionScale(std::string name = "");
		double ProductionScale(Amplitude::Channel name);
		
		std::vector<Particle> DecayPS(std::string name = "");
		std::vector<Particle> DecayPS(Amplitude::Channel name);
		std::vector<Particle> ProductionPS(TLorentzVector &vec,
						   std::string name = "");
		std::vector<Particle> ProductionPS(TLorentzVector &vec,
						   Amplitude::Channel name);
		
		void SetDecayChannel(std::string name);
		void SetProductionChannel(std::string name);
		Amplitude::Channel DecayChannel() const;
		Amplitude::Channel ProductionChannel() const;
		std::string DecayChannelName() const;
		std::string ProductionChannelName() const;

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
