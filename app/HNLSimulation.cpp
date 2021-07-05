#include <iostream>
#include <fstream>
#include <vector>
#include <cstring>
#include <getopt.h>

#include "tools/CardDealer.h"
#include "tools/RootFinding.h"

#include "detector/Tracker.h"
#include "detector/Driver.h"

#include "physics/Const.h"
#include "physics/Neutrino.h"
#include "physics/Decays.h"
#include "physics/DecaySpace.h"

#include "montecarlo/hnl.h"
#include "montecarlo/Sampler.h"
#include "montecarlo/Process.h"

int main(int argc, char** argv)
{
	const struct option longopts[] = 
	{
		{"channel", 	required_argument,	0, 'c'},
		{"threshold", 	required_argument,	0, 't'},
		{"massdepend", 	required_argument,	0, 'q'},
		{"efficiency", 	no_argument,		0, 'W'},
		{"dirac", 	no_argument,		0, 'r'},
		{"majorana", 	no_argument,		0, 'j'},
		{"help", 	no_argument,	 	0, 'h'},
		{0,	0, 	0,	0},
	};

	int index; 
	int iarg = 0;
	opterr = 1;
	
	//Initialize variables
	std::string channel = "ALL";
	std::string background;
	
	double mass = 0.;
	std::string ferm;
	double ue = 0., um = 0., ut = 0.;

	while((iarg = getopt_long(argc,argv, "m:c:E:M:T:rj", longopts, &index)) != -1)
	{
		switch(iarg)
		{
			case 'm':
				mass = std::strtod(optarg, NULL);
				break;
			case 'c':
				channel.assign(optarg);
				break;
			case 'E':
				ue = std::strtod(optarg, NULL);
				break;
			case 'M':
				um = std::strtod(optarg, NULL);
				break;
			case 'T':
				ut = std::strtod(optarg, NULL);
				break;
			case 'r':
				ferm = "d";
				break;
			case 'j':
				ferm = "m";
				break;
			default:
				break;
		}
	}

	//constructing the neutrino, engine and detector classes

	//output file name
	//
	if (channel.empty()) {
		std::cerr << "HNLSimulation: a channel must be selected with -c <channel>\n";
		return 0;
	}

	if (!ue && !um && !ut) {
		std::cerr << "HNLSimulation: at least one mixing must be selected with -E -M or -T\n";
		return 0;
	}

	if (ferm.empty()) {
		std::cerr << "HNLSimulation: fermion type must be selected with --dirac or --majorana\n";
		return 0;
	}

	CardDealer cd(argv[optind]);

	std::string config;
	if (!cd.Get("detector_configuration", config)) {
		std::cerr << "HNLSimulation: There is no detector configuration\n";
		return 0;
	}
	Tracker box(config);

	if (!cd.Get("flux_configuration", config)) {
		std::cerr << "HNLSimulation: There is no flux configuration\n";
		return 0;
	}
	Driver drive(config);

	Mixing mix(ue, um, ut);
	Decay::Channel chan = Decay::fromString(channel);

	if (!DecayRate::IsAllowed(mass, chan)) {
		std::cerr << "HNLSimulation: mass " << mass << " is not allowed for " << channel << " decay\n";
		return 0;
	}
	// MC config

	size_t nEvents;
	if (!cd.Get("number_events", nEvents))
		nEvents = 1e5;


	std::vector<size_t> nus;
	if (ferm == "m")
		nus = { Neutrino::majorana | Neutrino::left,
		        Neutrino::majorana | Neutrino::right };
	else if (ferm == "d") 
		nus = { Neutrino::dirac | Neutrino::particle | Neutrino::left,
			Neutrino::dirac | Neutrino::particle | Neutrino::right,
			Neutrino::dirac | Neutrino::antiparticle | Neutrino::left,
			Neutrino::dirac | Neutrino::antiparticle | Neutrino::right } ;

	// store samples for left and right neutrinos
	std::unordered_map<size_t, std::shared_ptr<TH1D> > samples;
	std::unordered_map<size_t, double> weights;
	std::unordered_map<size_t, DecaySpace> decayps;
	for (size_t opt : nus) {
		Neutrino N(mass, opt);
		auto smp = Sampler::Compute(box, drive, chan, mix, {N});
		if (!smp)
			continue;
		samples[opt] = smp;
		weights[opt] = smp->Integral();
		decayps.emplace(opt, std::move(N));
	}

	// clear neutrino vector
	if (!samples.size()) {
		std::cerr << "HNLSimulation: HNL decay of mass " << std::to_string(mass)
			<< " into " << channel << " is not allowed\n";
		return 0;
	}


	// output creation
	
	std::string outname;
	if (!cd.Get("output", outname)) {
		std::cerr << "HNLSimulation: no output \".root\" file specified in card\n";
		return 0;
	}

	if (outname.find(".root") == std::string::npos)
		outname += ".root";

	outname.insert(outname.find(".root"), "_" + channel);

	if (ue > 0.)
		outname.insert(outname.find(".root"), "_E");
	if (um > 0.)
		outname.insert(outname.find(".root"), "_M");
	if (ut > 0.)
		outname.insert(outname.find(".root"), "_T");
	
	outname.insert(outname.find(".root"), "_" + ferm);

	std::stringstream ssm;
	ssm << "_" << std::setw(4) << std::fixed << std::setfill('0')
		<< std::setprecision(0) << mass * 1000;
	outname.insert(outname.find(".root"), ssm.str());

	TFile out(outname.c_str(), "RECREATE");
	std::shared_ptr<hnl> data(new hnl);	// 

	std::cout << "Starting HNL simulation of " << nEvents << " events of "
		  << mass << " GeV neutrinos\nweights:";
	for (auto &w : weights) {
		std::cout << "\t" << w.second << "\t";
		w.second /= nEvents;	// update weight
		// divide by 2 because of helicity average!
	}
	std::cout << "\ndecay in channel " << channel
		  << " with mixing " << mix << "\n"
		  << "Creating file " << outname << "\n";

	// use volume ratios to determine probability
	std::unordered_map<std::string, double> mod_ratios;
	for (std::string mod : box.Modules())
		mod_ratios.emplace(mod, box.Volume(mod) / box.Volume());

	for (size_t id = 0; id < nEvents; ++id) {
		//std::cout << " Event " << id << "\n";
		for (const auto & s : samples) {
			// random energy for this neutrino
			double energy = (s.second)->GetRandom();
			// generate same event in all detector modules
			for (const auto & mrat : mod_ratios) {

				//neutrino probe
				Particle nu(12, energy, 0, 0, std::sqrt(energy*energy - mass*mass));
				Tracker::Event nu_evt = box.GenerateEvent(mrat.first, std::move(nu));

				auto prods = decayps[s.first].Generate(chan, nu_evt.first, mix);
				if (prods.first.size() != Decay::Ns(chan))
					continue; 
				std::vector<Tracker::Event> events;
				events.reserve(prods.first.size());
				for (Particle &part : prods.first) {
					Tracker::Event part_evt(std::move(part), nu_evt.second); //create particle
					if (std::abs(part_evt.first.Pdg() == 111)) //special treatments 
					{					//almost 100% into 2photons
						auto decay_res = box.Pi0Decay(std::move(part_evt));
						Tracker::Event gA_evt = std::move(decay_res[0]);
						Tracker::Event gB_evt = std::move(decay_res[1]);


						if (box.Reconstruct(gA_evt))
							events.push_back(std::move(gA_evt));
						if (box.Reconstruct(gB_evt))
							events.push_back(std::move(gB_evt));
					}
					else if (box.Reconstruct(part_evt))
						events.push_back(std::move(part_evt));
				}

				if (events.size() >= 2) {
					//make some simple misidentification of events
					events = box.MisIdentify(std::move(events));

					auto id = Process::Identify(chan, events);
					if (id == Process::Match::no_id)
						continue;

					data->Fill(weights[s.first] * prods.second * mrat.second,
						   (id == Process::Match::charge_id),
						   nu_evt, events[0], events[1]);
				}
			}
		}
	}
	std::cout << "Successfully simulated " << data->GetEntries() << " events!\n";

	out.cd();
	data->Write();
	for (const auto & s : samples) {
		std::string name = "sample" + std::to_string(s.first);
		s.second->Scale(1, "WIDTH");	// it looks nicer
		s.second->Write(name.c_str(), TObject::kOverwrite);
	}

	return 0;
}
