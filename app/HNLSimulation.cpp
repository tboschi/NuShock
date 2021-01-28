#include <iostream>
#include <fstream>
#include <vector>
#include <cstring>
#include <getopt.h>

#include "tools/CardDealer.h"
#include "tools/RootFinding.h"

#include "detector/Tracker.h"

#include "flux/Driver.h"
#include "flux/Sampler.h"

#include "physics/Const.h"
#include "physics/Neutrino.h"
#include "physics/PhaseSpace.h"

#include "montecarlo/hnl.h"


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
	size_t ferm = 0;
	std::string ferm_append;
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
				ferm = Neutrino::dirac;
				ferm_append = "_d";
				break;
			case 'j':
				ferm = Neutrino::majorana;
				ferm_append = "_m";
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

	if (ferm_append.empty())
		throw std::logic_error("HNLSimulation: fermion type must be selected with --dirac or --majorana");

	CardDealer cd(argv[optind]);

	std::string config;
	if (!cd.Get("detector_configuration", config))
		throw std::logic_error("There is no detector configuration\n");
	Tracker box(config);

	if (!cd.Get("flux_configuration", config))
		throw std::logic_error("There is no flux configuration\n");
	Driver drive(config);

	Mixing mix(UeFlag, UmFlag, UtFlag);
	mix *= 1e-4;
	Channel::Name chan = Channel::fromString(channel);

	std::vector<Neutrino> nus;
	if (ferm == Neutrino::majorana)
		nus = { Neutrino(mass, ferm | Neutrino::left),
			Neutrino(mass, ferm | Neutrino::right) };
	else // is dirac
		nus = { Neutrino(mass, ferm | Neutrino::particle | Neutrino::left),
			Neutrino(mass, ferm | Neutrino::particle | Neutrino::right),
			Neutrino(mass, ferm | Neutrino::antiparticle | Neutrino::left),
			Neutrino(mass, ferm | Neutrino::antiparticle | Neutrino::right) } ;

	Sampler smp(box, drive);
	// store samples for left and right neutrinos
	std::map<size_t, std::shared_ptr<TH1D> > samples;
	std::map<size_t, double> weights;
	std::map<size_t, PhaseSpace> ps;
	for (const Neutrino &N : nus)
		if (smp.Bind(chan, mix, N)) {
			samples[N.GetOptions()] = smp.MakeSampler(mix);
			weights[N.GetOptions()] = samples[N.GetOptions()]->Integral("WIDTH");
			ps.emplace(N.GetOptions(), std::move(N));
		}

	// reweighting
	if (ferm == Neutrino::majorana) {
		double sum = std::accumulate(weights.begin(), weights.end(), 0.,
				[](double sum, const std::pair<size_t, double> &w) {
					return sum + w.second;
					} );
		for (auto & w : weights)
			w.second /= sum;
	}
	else {
		for (bool o : {true, false}) {
			double sum = std::accumulate(weights.begin(), weights.end(), 0.,
					[&](double sum, const std::pair<size_t, double> &w) {
						if (Neutrino::IsParticle(w.first) == o)
							return sum + w.second;
						return sum;
						} );
			for (auto & w : weights)
				if (Neutrino::IsParticle(w.first) == o)
					w.second /= sum;
		}
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
	
	outname.insert(outname.find(".root"), ferm_append);

	std::stringstream ssm;
	ssm << "_" << std::setw(4) << std::fixed << std::setfill('0')
		<< std::setprecision(0) << mass * 1000;
	outname.insert(outname.find(".root"), ssm.str());

	// MC config

	size_t nEvents;
	if (!cd.Get("number_events", nEvents))
		nEvents = 1e5;

	TFile out(outname.c_str(), "RECREATE");
	std::shared_ptr<hnl> data(new hnl);	// 

	std::cout << "Starting HNL simulation of " << nEvents << " events of\n";
	for (const auto &n : nus)
		std::cout << "\t" << n << " (" << weights[n.GetOptions()] << ")\n";
	std::cout << "decay in channel " << channel
		  << " with mixing " << mix << "\n"
		  << "Creating file " << outname << "\n";

	for (size_t id = 0; id < nEvents; ++id) {
		for (const auto & s : samples) {
			double energy = s.second->GetRandom();
			for (const auto &mod : box.Modules()) {

				// use volume ratios to determine probability
				double rati = box.Volume(mod) / box.Volume();
				//neutrino probe
				Particle nu(12, energy, 0, 0, std::sqrt(energy*energy - mass*mass));
				Tracker::Event nu_evt = box.GenerateEvent(mod, std::move(nu));

				auto prods = ps[s.first].Generate(chan, nu_evt.first);
				std::vector<Tracker::Event> events;
				events.reserve(prods.first.size());
				for (Particle &part : prods.first) {
					Tracker::Event part_evt(std::move(part), nu_evt.second); //create particle
					if (std::abs(part_evt.first.Pdg() == 111))	//special treatments for pi0
					{						//almost 100% into 2photons
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
					if (events.size() >= 2)
						data->Fill(weights[s.first] * prods.second * rati,
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
