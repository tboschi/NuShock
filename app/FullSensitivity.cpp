#include <iostream>
#include <fstream>
#include <vector>
#include <cstring>
#include <getopt.h>

#include "tools/CardDealer.h"
#include "tools/RootFinding.h"
#include "tools/FeldmanCousins.h"

#include "detector/Detector.h"
#include "detector/Performance.h"
#include "detector/Driver.h"

#include "montecarlo/Sampler.h"

#include "physics/Const.h"
#include "physics/Neutrino.h"
#include "physics/Decays.h"
#include "physics/DecayRate.h"
#include "physics/ProductionRate.h"


void usage(char* argv0);
int main(int argc, char** argv)
{
	const struct option longopts[] = 
	{
		{"CL", 		required_argument,	0, 'C'},
		{"channels", 	required_argument,	0, 'c'},
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
	std::string channels = "";
	
	double CL = 0.90;
	std::string ferm;
	bool UeFlag = false, UmFlag = false, UtFlag = false;

	while((iarg = getopt_long(argc,argv, "C:c:w:EMTrjh", longopts, &index)) != -1)
	{
		switch(iarg)
		{
			case 'C':
				CL = std::strtod(optarg, NULL);
				break;
			case 'c':
				channels.assign(optarg);
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
			case 'h':
				usage(argv[0]);
				return 1;
			default:
				break;
		}
	}

	//constructing the neutrino, engine and detector classes

	if (channels.empty())
		throw std::logic_error("FullSens: at least one channel must be selected with -c <channel>");

	if (!UeFlag && !UmFlag && !UtFlag)
		throw std::invalid_argument("FullSens: at least one mixing must be selected with -E -M or -T");

	if (ferm.empty())
		throw std::logic_error("FullSens: fermion type must be selected with --dirac or --majorana");

	CardDealer cd(argv[optind]);

	std::string config;
	if (!cd.Get("detector_configuration", config))
		throw std::logic_error("There is no detector configuration\n");
	Detector box(config);

	// detector performance is described in input card
	if (!cd.Get("efficiency_cuts", config)) {
		std::cerr << "Sensitivity: no efficiency cut specified\n";
		config = "";
	}
	Performance eff(config);

	if (!cd.Get("flux_configuration", config))
		throw std::logic_error("There is no flux configuration\n");
	Driver drive(config);

	std::string outname;
	if (!cd.Get("output", outname))
		throw std::logic_error("FullSens: no output \".dat\" file specified in card");
	if (outname.find(".dat") == std::string::npos)
		outname += ".dat";

	std::replace(channels.begin(), channels.end(), ',', '_');
	outname.insert(outname.find(".dat"), "_" + channels);

	if (UeFlag)
		outname.insert(outname.find(".dat"), "_E");
	if (UmFlag)
		outname.insert(outname.find(".dat"), "_M");
	if (UtFlag)
		outname.insert(outname.find(".dat"), "_T");
	
	outname.insert(outname.find(".dat"), "_" + ferm);


	// from card
	size_t grid;	// from card
	if (!cd.Get("discretization", grid))
		grid = 500;
	std::vector<double> mass_range;
	if (!cd.Get("mass_range", mass_range))
		mass_range = {0.01, 2.00};

	std::vector<double> mix_range;
	if (!cd.Get("lmix_range", mix_range))
		mix_range = {-16., 0.};

	Mixing mix(UeFlag, UmFlag, UtFlag);
	// detection channels
	std::vector<Decay::Channel> chans;
	std::stringstream sschan(channels);
	std::string llchan;
	while (std::getline(sschan, llchan, '_'))
		chans.push_back(Decay::fromString(llchan));

	// store first result
	bool first = true;
	double m0, u0;

	std::vector<std::pair<double, double> > lo_limit, up_limit;
	lo_limit.reserve(grid);
	up_limit.reserve(grid);

	//Sampler sample(box, drive); //, Neutrino(mass, ferm), mix);
	std::ofstream out(outname.c_str());
	for (size_t g = 0; g < grid; ++g) {
		double mass = std::pow(mass_range[1] / mass_range[0],
				       g / double(grid)) * mass_range[0];

		// there is no channel for which decay is allowed
		if (std::none_of(chans.begin(), chans.end(),
				[&](Decay::Channel chan) { return DecayRate::IsAllowed(mass, chan); }))
			continue;

		std::cout << "mass " << mass << "\n";

		// store here all neutrinos needed for the analysis
		std::vector<Neutrino> nus;
		if (ferm == "m")
			nus = {Neutrino(mass, Neutrino::majorana)};
		else if (ferm == "d")
			nus = {Neutrino(mass, Neutrino::dirac | Neutrino::particle),
			       Neutrino(mass, Neutrino::dirac | Neutrino::antiparticle)};

		std::array<double, 2> extremes{{mix_range[1], mix_range[0]}};
		bool lo_lim = false, up_lim = false;
		for (const auto &chan : chans) {
			// create a function that generate a sample for given mixing
			auto sample = Sampler::Build(box, eff, drive, chan, mix, nus);

			double bkg = eff.Background(chan, (ferm == "d"), mass);
			double thr = FeldmanCousins::RejectBackground(bkg, CL, 0.01);

			// create function that computes number of events using sampler
			// lx is mixing in log scale
			auto func = [&](double lx) -> double {
				auto hist = sample(mix * std::pow(10., lx / 2.));
				return (hist ? hist->Integral() : 0.) - thr;
			};

			double half = RootFinding::CheckInterval(func, mix_range[0], mix_range[1], 1.e-4);
			double lim = RootFinding::BinarySearch(func, mix_range[0], half, 1.e-4);
			if (half > mix_range[0]) {
				lo_lim = true;
				if (lim < extremes[0])
					extremes[0] = lim;
				if (half < mix_range[1]) {
					up_lim = true;
					lim = RootFinding::BinarySearch(func, half, mix_range[1], 1.e-4);
					if (lim > extremes[1])
						extremes[1] = lim;
				}
			}
		}
		if (lo_lim) {
			lo_limit.emplace_back(mass, std::pow(10., extremes[0]));
			if (up_lim)
				up_limit.emplace_back(mass, std::pow(10., extremes[1]));
			else
				up_limit.emplace_back(mass, std::pow(10., mix_range[1]));
		}
	}

	for (auto i = lo_limit.begin(); i != lo_limit.end(); ++i)
		out << i->first << "\t" << i->second << "\n";
	for (auto i = up_limit.rbegin(); i != up_limit.rend(); ++i)
		out << i->first << "\t" << i->second << "\n";
	if (lo_limit.size()) {
		out << lo_limit[0].first << "\t" << lo_limit[0].second << "\n";
	}

	out.close();

	return 0;
}
	
void usage(char* argv0)
{
	std::cout << "Description" << std::endl;
	std::cout << "usage : " << std::endl;
	std::cout << argv0 << " [OPTIONS]" << std::endl;
	std::cout <<"\n  -d,  --detconfig [CONFIG]" << std::endl;
	std::cout << "\t\tDetector configuration file" << std::endl;
	std::cout <<"\n  -f,  --fluxconfig [CONFIG]" << std::endl;
	std::cout << "\t\tFlux configuration file" << std::endl;
	std::cout <<"\n  -o,  --output <path>" << std::endl;
	std::cout << "\t\tBase name of output file" << std::endl;
	std::cout <<"\n  -c,  --channels <chan>[,<chan2>,...]" << std::endl;
	std::cout << "\t\tDecay channels, multiple ones separated by ," << std::endl;
	std::cout <<"\n  -t,  --threshold [FLOAT]" << std::endl;
	std::cout << "\t\tThreshold intercept. If -q not set, then is mass independent. Defaul 2.44" << std::endl;
	std::cout <<"\n  -q,  --massdepend [FLOAT]" << std::endl;
	std::cout << "\t\tMass dependance of threshold. Needs -t flag. Defaul 0" << std::endl;
	std::cout <<"\n  -E,  -M,  -T" << std::endl;
	std::cout << "\t\tSelect which mixing element. Multiple selection allowed" << std::endl;
	std::cout <<"\n  -L,  -R,  -U" << std::endl;
	std::cout << "\t\tSelect neutrino polarisation." << std::endl;
	std::cout <<"\n  --dirac,  --majorana" << std::endl;
	std::cout << "\t\tSelect fermionic nature of neutrino. Default Majorana" << std::endl;
	std::cout <<"\n  -P, -A {incompatible with --majorana}" << std::endl;
	std::cout << "\t\tSelect either neutrino or antineutrino components for dirac. Default, both are used." << std::endl;
	std::cout <<"\n  -h,  --help" << std::endl;
	std::cout << "\t\tPrint this message and exit" << std::endl;
}
