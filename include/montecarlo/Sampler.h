/* Sampler namespace is used to calculate number of HNL events
 * N = flux * decay prob * efficiency
 *
 * It combines the above elements:
 * 	- HNL flux
 * 	- decay probability
 * 	- detection efficiency if required
 */

#ifndef SAMPLER_H
#define SAMPLER_H

#include <iostream>
#include <string>
#include <fstream>
#include <cstring>
#include <functional>

#include "detector/Detector.h"
#include "detector/Performance.h"
#include "detector/Driver.h"

#include "physics/Neutrino.h"
#include "physics/Mixings.h"
#include "physics/DecayRate.h"


namespace Sampler
{
	using SamplerFn = std::function<std::shared_ptr<TH1D> (Mixing mix)>;
	SamplerFn Build(const Detector &box,
			const Performance &eff, const Driver &drive,
			const std::vector<Decay::Channel> &chans, Mixing mix,
			std::vector<Neutrino> nus);
	SamplerFn Build(const Detector &box,
			const Driver &drive,
			const std::vector<Decay::Channel> &chans, Mixing mix,
			std::vector<Neutrino> nus);
	SamplerFn Build(const Detector &box,
			const Performance &eff, const Driver &drive,
			Decay::Channel chan, Mixing mix,
			std::vector<Neutrino> nus);
	SamplerFn Build(const Detector &box,
			const Driver &drive,
			Decay::Channel chan, Mixing mix,
			std::vector<Neutrino> nus);

	std::shared_ptr<TH1D> Compute(const Detector &box, const Performance &eff, const Driver &drive,
				std::vector<Decay::Channel> &chans, Mixing mix, std::vector<Neutrino> nus);
	std::shared_ptr<TH1D> Compute(const Detector &box, const Driver &drive,
				std::vector<Decay::Channel> &chans, Mixing mix, std::vector<Neutrino> nus);
	std::shared_ptr<TH1D> Compute(const Detector &box, const Performance &eff, const Driver &drive,
				Decay::Channel chan, Mixing mix, std::vector<Neutrino> nus);
	std::shared_ptr<TH1D> Compute(const Detector &box, const Driver &drive,
				Decay::Channel chan, Mixing mix, std::vector<Neutrino> nus);

	SamplerFn Factory(const Detector &box,
		std::vector<std::pair<DecayRate, Spectrum> > &&spec_rate,
		const Performance &eff, const std::vector<Decay::Channel> &chan);

	SamplerFn Factory(const Detector &box,
		std::vector<std::pair<DecayRate, Spectrum> > &&spec_rate,
		const std::vector<Decay::Channel> &chan);
};

#endif
