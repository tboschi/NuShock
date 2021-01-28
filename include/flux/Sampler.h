/*
 * Author: Tommaso Boschi
 */

#ifndef SAMPLER_H
#define SAMPLER_H

#include <iostream>
#include <string>
#include <fstream>
#include <cstring>
#include <sstream>

#include "detector/Detector.h"
#include "physics/Neutrino.h"
#include "physics/Mixings.h"
#include "physics/DecayRates.h"

#include "flux/Driver.h"

class Sampler
{
	public:

		Sampler(const Detector &box, const Driver &drive);
				//Neutrino nu = Neutrino(), Mixing mix = Mixing());
		bool Bind(Channel::Name chan, Mixing mix = Mixing(), Neutrino nu = Neutrino());
		std::shared_ptr<TH1D> MakeSampler(const Mixing &mix = Mixing());

	private:

		Detector _box;
		Driver _drive;

		Channel::Name _chan;
		Neutrino _N;
		DecayRates _rate;
};

#endif
