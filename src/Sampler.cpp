#include "montecarlo/Sampler.h"

// call sampler function only after making flux
namespace Sampler {
	// return a sampler generator
	SamplerFn Build(const Detector &box,
			const Performance &eff, const Driver &drive,
			const std::vector<Decay::Channel> &chans, Mixing mix,
			std::vector<Neutrino> nus)
	{
		// combine decay rates and spectra in a single object
		// create rates with correct neutrino type
		std::vector<std::pair<DecayRate, Spectrum> > spec_rate;
		spec_rate.reserve(nus.size());
		for (auto &N : nus) {
			auto spec = drive.MakeSpectrum(N, mix);
			if (!spec(mix)) // empty
				continue;
			spec_rate.emplace_back(N, spec);
		}

		// create sampler-making function
		// mutable because it has to remember decay rates
		return Factory(box, std::move(spec_rate), eff, chans);
	}

	SamplerFn Build(const Detector &box,
			const Driver &drive,
			const std::vector<Decay::Channel> &chans, Mixing mix,
			std::vector<Neutrino> nus)
	{
		// combine decay rates and spectra in a single object
		// create rates with correct neutrino type
		std::vector<std::pair<DecayRate, Spectrum> > spec_rate;
		spec_rate.reserve(nus.size());
		for (auto &N : nus) {
			auto spec = drive.MakeSpectrum(N, mix);
			if (!spec(mix)) // empty
				continue;
			spec_rate.emplace_back(N, spec);
		}

		// create sampler-making function
		// mutable because it has to remember decay rates
		return Factory(box, std::move(spec_rate), chans);
	}

	SamplerFn Build(const Detector &box,
			const Performance &eff, const Driver &drive,
			Decay::Channel chan, Mixing mix,
			std::vector<Neutrino> nus)
	{
		std::vector<Decay::Channel> chans = {chan};
		return Build(box, eff, drive, chans, mix, nus);
	}

	SamplerFn Build(const Detector &box,
			const Driver &drive,
			Decay::Channel chan, Mixing mix,
			std::vector<Neutrino> nus)
	{
		std::vector<Decay::Channel> chans = {chan};
		return Build(box, drive, chans, mix, nus);
	}

	// compute directly a single sampler
	std::shared_ptr<TH1D> Compute(const Detector &box, const Performance &eff, const Driver &drive,
				std::vector<Decay::Channel> chans, Mixing mix, std::vector<Neutrino> nus)
	{
		auto func = Build(box, eff, drive, chans, mix, nus);
		return func(mix);
	}

	// compute directly a single sampler
	std::shared_ptr<TH1D> Compute(const Detector &box, const Driver &drive,
				std::vector<Decay::Channel> chans, Mixing mix, std::vector<Neutrino> nus)
	{
		auto func = Build(box, drive, chans, mix, nus);
		return func(mix);
	}

	// compute directly a single sampler
	std::shared_ptr<TH1D> Compute(const Detector &box, const Performance &eff, const Driver &drive,
				Decay::Channel chan, Mixing mix, std::vector<Neutrino> nus)
	{
		std::vector<Decay::Channel> chans = {chan};
		auto func = Build(box, eff, drive, chans, mix, nus);
		return func(mix);
	}

	// compute directly a single sampler
	std::shared_ptr<TH1D> Compute(const Detector &box, const Driver &drive,
				Decay::Channel chan, Mixing mix, std::vector<Neutrino> nus)
	{
		std::vector<Decay::Channel> chans = {chan};
		auto func = Build(box, drive, chans, mix, nus);
		return func(mix);
	}

	SamplerFn Factory(const Detector &box,
			std::vector<std::pair<DecayRate, Spectrum> > &&spec_rate,
			const Performance &eff,
			const std::vector<Decay::Channel> &chans) {

		return [&, spec_rate, chans](Mixing mix) mutable {

			std::shared_ptr<TH1D> all = nullptr;

			// loop over spectra
			for (auto & sr : spec_rate) {
				DecayRate &r = sr.first;
				Spectrum &s = sr.second;
				// copy neutrino object
				Neutrino nu = r.GetNeutrino();

				// compute flux
				std::shared_ptr<TH1D> flux = s(mix);
				if (!flux) // skip if empty
					continue;

				std::shared_ptr<TH1D> profile = nullptr;
				double branch = 0.;
				for (auto chan : chans) {
					auto prof = eff.EfficiencyProfile(chan, nu);
					if (!prof) {
						branch += r.Branch(chan, mix);
						continue;
					}

					prof->Scale(r.Branch(chan, mix));
					if (!profile)
						profile = prof;
					else
						profile->Add(prof.get());
				}

				double scale = box.Scaling() * box.POTs();

				for (int bin = 1; bin <= flux->GetNbinsX(); ++bin) {
					double energy = flux->GetBinCenter(bin);
					nu.SetE(energy);
					if (nu.EKin() <= 0. || nu.Beta() > 1.) {
 						// unphysical neutrino
						flux->SetBinContent(bin, 0.);
						continue;
					}
					double tby = r.Total(mix) / nu.Beta() / nu.Gamma();
					double eff = (profile ? profile->Interpolate(energy) : branch);

					// update flux content
					double cont = flux->GetBinContent(bin) * flux->GetBinWidth(bin)
						    * box.Probability(tby) * eff * scale;
					flux->SetBinContent(bin, cont);
					/*
					std::cout << "Sampler: " << energy
						  << "\t" << r.Total(mix)
						  << "\t" << nu.Beta() * nu.Gamma()
						  << "\t" << box.Probability(tby)
						  << "\t" << scale
						  << "\t" << eff
						  << " = " << cont << "\n";
					*/
				}

				if (!all)
					all = flux;
				else
					all->Add(flux.get());
			}
			return all;
		};
	}

	// same as factory above, but without Performance object
	// bad coding practice, but quickest fix I thought about
	SamplerFn Factory(const Detector &box,
			std::vector<std::pair<DecayRate, Spectrum> > &&spec_rate,
			const std::vector<Decay::Channel> &chans) {

		return [&, spec_rate, chans](Mixing mix) mutable {

			std::shared_ptr<TH1D> all = nullptr;

			// loop over spectra
			for (auto & sr : spec_rate) {
				DecayRate &r = sr.first;
				Spectrum &s = sr.second;
				// copy neutrino object
				Neutrino nu = r.GetNeutrino();

				// compute flux
				std::shared_ptr<TH1D> flux = s(mix);
				if (!flux) // skip if empty
					continue;

				double branch = std::accumulate(chans.begin(), chans.end(), 0.,
						[&](double sum, const Decay::Channel &chan) {
							return sum + r.Branch(chan, mix);
						});
				double scale = box.Scaling() * box.POTs() * branch;

				for (int bin = 1; bin <= flux->GetNbinsX(); ++bin) {
					double energy = flux->GetBinCenter(bin);
					nu.SetE(energy);
					if (nu.EKin() <= 0. || nu.Beta() > 1.) {
 						// unphysical neutrino
						flux->SetBinContent(bin, 0.);
						continue;
					}
					double tby = r.Total(mix) / nu.Beta() / nu.Gamma();

					// update flux content
					double cont = flux->GetBinContent(bin) * flux->GetBinWidth(bin)
						    * box.Probability(tby) * scale;
					flux->SetBinContent(bin, cont);
					/*
					std::cout << "Sampler: " << energy
						  << "\t" << r.Total(mix)
						  << "\t" << r.Branch(chan, mix)
						  << "\t" << nu.Beta() * nu.Gamma()
						  << "\t" << box.Probability(tby)
						  << "\t" << eff.Efficiency(chan, nu)
						  << " = " << cont << "\n";
					*/
				}

				if (!all)
					all = flux;
				else
					all->Add(flux.get());
			}
			return all;
		};
	}
}
