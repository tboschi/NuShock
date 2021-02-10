#include <iostream>
#include <cstring>
#include <getopt.h>

#include "physics/Neutrino.h"
#include "physics/Mixings.h"
#include "physics/ProductionRate.h"

void usage(char* argv0)
{
	std::cout << "usage: " << argv0 << " <mass> <channel> [<options>]\n\n";
	std::cout << "Compute production scale for given channel of HNL with given mass\n";

	std::cout << "\nOptional parameters:\n"
		  << "  -r, --dirac\tuse a Dirac HNL [default]\n"
		  << "  -j, --majorana\tuse a Majorana HNL\n"
		  << "  -L, --left\tset Left helicity\n"
		  << "  -R, --right\tset Right helicity\n"
		  << "  -E, --ue\tset electron mixing value [default = 1]\n"
		  << "  -M, --um\tset muon mixing value [default = 1]\n"
		  << "  -T, --ut\tset tau mixing value [default = 1]\n"
		  << "  -q, --quiet\tonly print production scale value\n"
		  << "  -h, --help\tprint this message and exit\n\n";
}

int main(int argc, char** argv)
{

	const struct option longopts[] = 
	{
		{"dirac", 	no_argument,		0, 'r'},
		{"majorana", 	no_argument,		0, 'j'},
		{"left", 	no_argument,		0, 'L'},
		{"right", 	no_argument,		0, 'R'},
		{"channel", 	required_argument, 	0, 'c'},
		{"mass", 	required_argument, 	0, 'm'},
		{"ue", 		required_argument, 	0, 'E'},
		{"um", 		required_argument, 	0, 'M'},
		{"ut", 		required_argument, 	0, 'T'},
		{"help",	required_argument, 	0, 'h'},
		{0,	0, 	0,	0},
	};

	int index; 
	int iarg = 0;
	opterr = 1;
	
	size_t ferm = Neutrino::dirac, helix = Neutrino::left;
	Production::Channel chan;
	double mass = -1;
	bool verbose = true;
	double ue = 1., um = 1., ut = 1.;
	while((iarg = getopt_long(argc,argv, "m:c:E:M:T:LRUrjqh", longopts, &index)) != -1)	
	{
		switch(iarg)
		{
			case 'c':
				break;
			case 'm':
				mass = std::strtod(optarg, NULL);
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
				ferm = Neutrino::dirac;
				break;
			case 'j':
				ferm = Neutrino::majorana;
				break;
			case 'L':
				helix = Neutrino::left;
				break;
			case 'R':
				helix = Neutrino::right;
				break;
			case 'U':
				helix = Neutrino::unpolarized;
				break;
			case 'q':
				verbose = false;
				break;
			case 'h':
			default:
				usage(argv[0]);
				return 1;
		}
	}

	mass = std::strtod(argv[optind], NULL);
	chan = Production::fromString(std::string(argv[optind+1]));

	// unpolairsed?
	Mixing mix(ue, um, ut);
	Neutrino nu(mass, ferm | helix);

	if (verbose) {
		std::cout << "Neutrino " << nu << "\n";
		std::cout << "Mixing " << mix << "\n";
		std::cout << "Production scale for channel " << Production::toString(chan) << ": ";
	}

	if (!ProductionRate::IsAllowed(nu, chan)) {
		if (verbose)
			std::cout << "not allowed\n";
		std::cout << 0.0;
		return 1;
	}


	ProductionRate hnl(nu);
	std::cout << hnl.Scale(chan, mix) << "\n";

	return 0;
}
