/*
 * PhaseSpace generator for simulation of sterile neutrino decays (it is not needed for production)
 * It also returns correct 
 *
 * Author: Tommaso Boschi
 */

#ifndef PHASESPACE_H
#define PHASESPACE_H

// base class
#include "physics/Amplitude.h"

//ROOT include
#include "TGenPhaseSpace.h"

#include <random>
#include "physics/Particle.h"
#include "physics/Neutrino.h"

#include "tools/Optimization.h"
#include "tools/RNG.h"

template <typename C>
class PhaseSpace : public Amplitude
{
	protected:

		PhaseSpace(Neutrino N = Neutrino()) : Amplitude(std::move(N)),
						      _genps(new TGenPhaseSpace()) { };
		~PhaseSpace() {
			delete _genps;
		}

		std::pair<std::vector<Particle>, double> Generate(C chan, const Particle &part,
								  const Mixing &mix = Mixing()) {
			TLorentzVector frame = static_cast<TLorentzVector>(part);
			return Generate(chan, frame, mix);
		}

		// defined in derived class
		//virtual std::pair<std::vector<Particle>, double>
		//	Generate(C chan, const TLorentzVector &frame, const Mixing &mix = Mixing()) = 0;

		//bool SetDecay(C chan);
		//double Gamma(Channel::Name chan, const Mixing &mix = Mixing());


		//KINEMATICS

		double Kinematic2B() {
			if (_genps->GetNt() < 2)
				throw std::logic_error("Kinematic_2B error");

			TLorentzVector vec1 = *(_genps->GetDecay(1));
			return std::cos(vec1.Theta());
		}

		//the mandelstam variables refer to a specific vector
		//s -> p3 //t -> p2 //u -> p1
		std::array<double, 6> Kinematic3B()
		{
			if (_genps->GetNt() < 3)
				throw std::logic_error("Kinematic_3B error");

			TLorentzVector vec1 = *(_genps->GetDecay(0));
			TLorentzVector vec2 = *(_genps->GetDecay(1));
			TLorentzVector vec3 = *(_genps->GetDecay(2));

			TLorentzVector rest(0, 0, 0, _m_parent);
			TLorentzVector vec_u = rest - vec1;
			TLorentzVector vec_t = rest - vec2;
			TLorentzVector vec_s = rest - vec3;

			std::array<double, 6> ret;

			ret[0] = vec_u.M2() / std::pow(_m_parent, 2); // s
			ret[1] = vec_t.M2() / std::pow(_m_parent, 2); // t
			ret[2] = vec_s.M2() / std::pow(_m_parent, 2); // u

			ret[3] = std::cos(vec1.Theta());	// cos0u 
			ret[4] = std::cos(vec2.Theta());	// cos0t
			ret[5] = std::cos(vec3.Theta());	// cos0s

			return ret;
		}

		TGenPhaseSpace *_genps;
};

#endif
