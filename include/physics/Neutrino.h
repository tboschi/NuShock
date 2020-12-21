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
		enum Options
		{
			left         = 0,
			right        = 1,
			polarised    = 0,
			unpolarised  = 2,
			dirac	     = 0,
			majorana     = 4,
			particle     = 0,
			antiparticle = 8
		};

		Neutrino(double mass = 0, unsigned int opts = 2);
		Neutrino(const Neutrino &N);

		void SetParent(Amplitude *Object);

		double DecayThreshold();
		double DecayThreshold(const std::string &name);
		double DecayThreshold(Channel::Name chan = Channel::_undefined);
		bool IsDecayAllowed();
		bool IsDecayAllowed(const std::string &name);
		bool IsDecayAllowed(Channel::Name chan = Channel::_undefined);
		void DecayChannels(std::vector<std::string> &vChan);

		double DecayTotal();
		double DecayWidth();
		double DecayWidth(const std::string &name);
		double DecayWidth(Channel::Name chan = Channel::_undefined);
		double DecayBranch();
		double DecayBranch(const std::string &name);
		double DecayBranch(Channel::Name chan = Channel::_undefined);

		double ProductionThreshold();
		double ProductionThreshold(const std::string &name);
		double ProductionThreshold(Channel::Name chan = Channel::_undefined);
		bool IsProductionAllowed();
		bool IsProductionAllowed(const std::string &name);
		bool IsProductionAllowed(Channel::Name chan = Channel::_undefined);
		void ProductionChannels(std::vector<std::string> &vChan);
		double ProductionWidth();
		double ProductionWidth(const std::string &name);
		double ProductionWidth(Channel::Name chan = Channel::_undefined);
		double ProductionScale();
		double ProductionScale(const std::string &name);
		double ProductionScale(Channel::Name chan = Channel::_undefined);
		
		std::vector<Particle> DecayPS();
		std::vector<Particle> DecayPS(const std::string &name);
		std::vector<Particle> DecayPS(Channel::Name chan = Channel::_undefined);
		std::vector<Particle> ProductionPS(const TLorentzVector &vec);
		std::vector<Particle> ProductionPS(const TLorentzVector &vec,
						   const std::string &name);
		std::vector<Particle> ProductionPS(const TLorentzVector &vec,
						   Channel::Name chan = Channel::_undefined);
		
		void SetDecayChannel(const std::string &name);
		void SetProductionChannel(const std::string &name);
		Channel::Name DecayChannel() const;
		Channel::Name ProductionChannel() const;
		std::string DecayChannelName() const;
		std::string ProductionChannelName() const;

		//void SetMass(double Mass = 0.0);
		void SetMixings(double *Mixings);
		void SetMixings(double Ue, double Um, double Ut);

		void SetOptions(size_t opts);
		void AddOptions(size_t opts);
		
		//double Mass();
		double* Mixings();
		double Ue(int E = 1.0) const;
		double Um(int E = 1.0) const;
		double Ut(int E = 1.0) const;
		//double Energy();
		//double EnergyKin();

		static int Helicity(const size_t &opts) const {
			if (opts & Neutrino::unpolarised)
				return 0;
			return 2 * int(opts | Neutrino::right) - 1;
		}

		static bool IsDirac(const size_t &opts) const {
			return !(opts & Neutrino::majorana);
		}

		static bool IsMajorana(const size_t &opts) const {
			return !IsDirac(opts);
		}

		// always true if majorana
		static bool IsParticle(const size_t &opts) const {
			return IsMajorana(opts) || !(opts & Neutrino::antiparticle);
		}

		// always true if majorana
		static bool IsAntiparticle(const size_t &opts) const {
			return IsMajorana(opts) || !IsParticle(opts);
		}

	private:
		DecayRates *theDecayRates;
		Production *theProduction, *theProdLightN;
		PhaseSpace *thePhaseSpace;

		Channel::Name _decay, _production;

		//double fMass;
		//double fEnergy;

		std::array<double, 3> _mixings;
};

#endif
