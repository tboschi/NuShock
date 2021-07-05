#include "montecarlo/Process.h"

namespace Process {

	// verify that particle contains the particle sought for, in mProc
	// return true (first) if process is identified regardless of charge
	Match Identify(Decay::Channel chan, std::initializer_list<Particle> ps) {
		std::vector<Particle> parts(std::move(ps));
		return Identify(chan, parts);
	}

	Match Identify(Decay::Channel chan, const std::vector<Tracker::Event> &evts) {
		std::vector<Particle> parts;
		parts.reserve(evts.size());
		std::transform(evts.begin(), evts.end(), std::back_inserter(parts),
				[](const Tracker::Event &evt) { return evt.first; });
		return Identify(chan, parts);
	}

	// verify that particle contains the particle sought for, in mProc
	// return true (first) if process is identified regardless of charge
	Match Identify(Decay::Channel chan, const std::vector<Particle> &parts)
	{
		if (chan == Decay::Channel::ExpALL) {
			auto dets = Decay::Detections();
			for (const auto &c : Decay::Detections()) {
				auto id = Identify(c, parts);
				if (id != Match::no_id)
					return id;
			}
			return Match::no_id;
		}

		std::vector<int> proc = Process::Pdgs(chan);
		if (parts.size() != proc.size())
			return Match::no_id;

		// find which one matches first particle
		int pdg0 = parts.front().Pdg();
		int q0 = parts.front().Q();
		auto p = std::find_if(proc.begin(), proc.end(),
				[&](int pdg) { return std::abs(pdg) == std::abs(pdg0); });

		if (p == proc.end())	// no match
			return Match::no_id;

		int sign = pdg0 * (*p) < 0 ? -1 : 1;	// check if opposite sign

		// if signs match, potential charge ID
		bool charge = Particle::Q(*p * sign) == q0;

		proc.erase(p);
		for (auto ie = parts.begin() + 1; ie != parts.end(); ++ie) {
			// if cant find equivalent
			int pdgi = ie->Pdg();
			int qi = ie->Q();

			if (charge) {
				// verify charge match
				p = std::find_if(proc.begin(), proc.end(),
					[&](int pdg) -> bool {
						return (pdg*sign == pdgi // match pdg and charge
						     && Particle::Q(pdg*sign) == qi); } );
				charge = (p != proc.end());
			}

			// continue with pdg-only identification
			p = std::find_if(proc.begin(), proc.end(),
				[&](int pdg) -> bool {
					return std::abs(pdg) == std::abs(pdgi); } );

			if (p == proc.end())
				return Match::no_id;
			proc.erase(p);
		}

		// all match!
		return charge ? Match::charge_id : Match::pdg_id;
	}

	// list of particles which can be detected by detector
	bool IsDetectable(const Particle &p)
	{
		switch (std::abs(p.Pdg()))
		{
			case 11:
			case 13:
			case 15:
			case 22:
			case 111:
			case 211:
			case 321:
				return true;
			default:
				return false;
		}
	}

	// return particles which are signature of the decay process
	// only detectable particles are permitted
	// use sign as if particle HNL decayed (not antiparticle)
	std::vector<int> Pdgs(Decay::Channel chan) {
		switch(chan) {
			//DECAYS
			//case Channel::ALL:
			case Decay::Channel::nnn:
				return {};
			case Decay::Channel::nGAMMA:
				return {22};
			case Decay::Channel::nEE:
				return {-11, 11};
			case Decay::Channel::nMM:
				return {-13, 13};
			case Decay::Channel::nEM:
				return {-11, 13};
			case Decay::Channel::nET:
				return {-11, 15};
			case Decay::Channel::nMT:
				return {-13, 15};
			case Decay::Channel::nPI0:
				return {22, 22};
			case Decay::Channel::EPI:
				return {11, 211};
			case Decay::Channel::MPI:
				return {13, 211};
			case Decay::Channel::TPI:
				return {15, 211};
			case Decay::Channel::EKA:
				return {11, 321};
			case Decay::Channel::MKA:
				return {13, 321};
			default:
				throw std::invalid_argument("Process channel "
							+ Decay::toString(chan) + " is unknwon");
		}
	}
}
