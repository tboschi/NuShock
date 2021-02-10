#include "montecarlo/Sampler.h"

// call sampler function only after making flux
namespace Sampler {
	// return a sampler generator
	std::function<std::shared_ptr<TH1D> (Mixing mix)> Build(const Detector &box,
			const Driver &drive, Decay::Channel chan, Mixing mix,
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
		auto func = [&, spec_rate, chan](Mixing mix) mutable {

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

				for (int bin = 1; bin <= flux->GetNbinsX(); ++bin) {
					nu.SetE(flux->GetBinCenter(bin));
					if (nu.EKin() <= 0. || nu.Beta() > 1.) {
 						// unphysical neutrino
						flux->SetBinContent(bin, 0.);
						continue;
					}
					double tby = r.Total(mix) / nu.Beta() / nu.Gamma();

					// update flux content
					double cont = flux->GetBinContent(bin) * flux->GetBinWidth(bin);
					//std::cout << "spec " << flux->GetBinLowEdge(bin)
					//	  << "\t" << cont  << "\t" << r.Total(mix)
					//	  << "\t" << r.Branch(chan, mix)
					//	  << "\t" << nu.Beta() * nu.Gamma()
					//	  << "\t" << box.Probability(tby)
					//	  << "\t" << _box.Efficiency(chan, N.M(), N.E())
					//	  << "\t = " << flux->GetBinContent(bin) << "\n";
					flux->SetBinContent(bin, cont
							* box.Probability(tby));
							//* box.Efficiency(chan, is->N.M(), is->N.E())
				}

				// scale by box dimension and branching ratio
				flux->Scale(box.Scaling() * box.POTs() * r.Branch(chan, mix));
				if (!all)
					all = flux;
				else
					all->Add(flux.get());
			}
			return all;
		};

		return func;
	}

	// compute directly a single sampler
	std::shared_ptr<TH1D> Compute(const Detector &box, const Driver &drive,
				Decay::Channel chan, Mixing mix, std::vector<Neutrino> nus)
	{
		auto func = Build(box, drive, chan, mix, std::move(nus));
		return func(mix);
	}

	// compute directly a single sampler
	std::shared_ptr<TH1D> Compute(const Detector &box, const Driver &drive,
				Decay::Channel chan, Mixing mix, Neutrino nu)
	{
		std::vector<Neutrino> nus = {std::move(nu)};
		return Compute(box, drive, chan, mix, nus);
	}
}
