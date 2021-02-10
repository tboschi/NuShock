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
	bool UeFlag = false, UmFlag = false, UtFlag = false;

	while((iarg = getopt_long(argc,argv, "m:c:EMTrj", longopts, &index)) != -1)
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
				UeFlag = true;
				break;
			case 'M':
				UmFlag = true;
				break;
			case 'T':
				UtFlag = true;
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
	if (channel.empty())
		throw std::logic_error("HNLSimulation: a channel must be selected with -c <channel>");

	if (!UeFlag && !UmFlag && !UtFlag)
		throw std::invalid_argument("HNLSimulation: at least one mixing must be selected with -E -M or -T");

	if (ferm.empty())
		throw std::logic_error("HNLSimulation: fermion type must be selected with --dirac or --majorana");

	CardDealer cd(argv[optind]);

	std::string config;
	if (!cd.Get("detector_configuration", config))
		throw std::logic_error("HNLSimulation: There is no detector configuration\n");
	Tracker box(config);

	if (!cd.Get("flux_configuration", config))
		throw std::logic_error("HNLSimulation: There is no flux configuration\n");
	Driver drive(config);

	Mixing mix(UeFlag, UmFlag, UtFlag);
	mix *= 1e-4;
	Decay::Channel chan = Decay::fromString(channel);

	if (!DecayRate::IsAllowed(mass, chan))
		throw std::invalid_argument("HNLSimulation: mass " + std::to_string(mass)
						+ "does not allow for decay");
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
		auto smp = Sampler::Compute(box, drive, chan, mix, N);
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
		return 1;
	}


	// output creation
	
	std::string outname;
	if (!cd.Get("output", outname))
		throw std::logic_error("HNLSimulation: no output \".root\" file specified in card");

	if (outname.find(".root") == std::string::npos)
		outname += ".root";

	outname.insert(outname.find(".root"), channel);

	if (UeFlag)
		outname.insert(outname.find(".root"), "_E");
	if (UmFlag)
		outname.insert(outname.find(".root"), "_M");
	if (UtFlag)
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
		for (const auto & s : samples) {
			// random energy for this neutrino
			double energy = s.second->GetRandom();
			// generate same event in all detector modules
			for (const auto & mrat : mod_ratios) {

				//neutrino probe
				Particle nu(12, energy, 0, 0, std::sqrt(energy*energy - mass*mass));
				Tracker::Event nu_evt = box.GenerateEvent(mrat.first,
								          std::move(nu));

				auto prods = decayps[s.first].Generate(chan, nu_evt.first);
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

				if (events.size() >= 2)
					//make some simple misidentification of events
					events = box.MisIdentify(std::move(events));
				// if still ok then fill
				if (events.size() >= 2)
					data->Fill(weights[s.first] * prods.second, // * mrat.second,
							nu_evt, events[0], events[1]);
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
