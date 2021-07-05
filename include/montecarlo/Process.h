#ifndef PROCESS_H
#define PROCESS_H

#include <string>
#include <vector>
#include <stdexcept>

#include "physics/Const.h"
#include "physics/Decays.h"
#include "physics/Particle.h"

#include "detector/Tracker.h"

namespace Process
{
	enum class Match {
		no_id = 0,
		pdg_id,
		charge_id
	};

	Match Identify(Decay::Channel chan, std::initializer_list<Particle> ps);
	Match Identify(Decay::Channel chan, const std::vector<Tracker::Event> &evts);
	Match Identify(Decay::Channel chan, const std::vector<Particle> &parts);

	bool IsDetectable(const Particle &p);
	std::vector<int> Pdgs(Decay::Channel chan);
}

inline std::ostream & operator<<(std::ostream &os, const Process::Match &id) {
	switch (id) {
		case Process::Match::no_id:
			return os << "no id";
		case Process::Match::pdg_id:
			return os << "particle id";
		case Process::Match::charge_id:
			return os << "charge id";
		default:
			return os;
	}
}

#endif
