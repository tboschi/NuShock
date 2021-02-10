#ifndef SAMPLER_H
#define SAMPLER_H

#include <iostream>
#include <string>
#include <fstream>
#include <cstring>
#include <functional>

#include "detector/Detector.h"
#include "detector/Driver.h"

#include "physics/Neutrino.h"
#include "physics/Mixings.h"
#include "physics/DecayRate.h"


namespace Sampler
{
	std::function<std::shared_ptr<TH1D> (Mixing mix)> Build(const Detector &box,
				const Driver &drive, Decay::Channel chan, Mixing mix,
				std::vector<Neutrino> nus);
	std::shared_ptr<TH1D> Compute(const Detector &box, const Driver &drive,
				Decay::Channel chan, Mixing mix, std::vector<Neutrino> nus);
	std::shared_ptr<TH1D> Compute(const Detector &box, const Driver &drive,
				Decay::Channel chan, Mixing mix, Neutrino nu);
};

#endif
