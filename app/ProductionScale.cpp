#include <iostream>
#include <fstream>
#include <vector>
#include <cstring>
#include <getopt.h>

#include "tools.h"
#include "physics.h"

void Usage(char* argv0);
int main(int argc, char** argv)
{

	const struct option longopts[] = 
	{
		{"output", 	required_argument, 	0, 'o'},
		{"channel", 	required_argument, 	0, 'c'},
		{"help", 	no_argument,	 	0, 'h'},
		{0,	0, 	0,	0},
	};

	int index; 
	int iarg = 0;
	opterr = 1;
	
	std::string channel;
	std::ofstream outf;
	double ue = 0.0, um = 0.0, ut = 0.0;
	
	while((iarg = getopt_long(argc,argv, "c:o:E:M:T:h", longopts, &index)) != -1)	
	{
		switch(iarg)
		{
			case 'c':
				channel.assign(optarg);
				break;
			case 'o':
				outf.open(optarg);
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
			case 'h':
			default:
				Usage(argv[1]);
				return 1;
		}
	}

	//To have an output, depending on usage usage
	std::ostream &out = (outf.is_open()) ? outf : std::cout;

	/* Neutrino options available
	 * Dirac, Majorana, Antiparticle	for fermionic nature
	 * Left, Right, Unpolarised		for helicity
	 *
	 * Concate with bitwise OR operator '|'
	 *
	 */
	
	//Neutrino nu(0, Neutrino::Dirac | Neutrino::Unpolarised );
	Neutrino nu(0, Neutrino::Majorana | Neutrino::Unpolarised );

	nu.SetProductionChannel(channel);
	nu.SetMixings(sqrt(ue), sqrt(um), sqrt(ut));

	out << "#Mass\tScale\n";
	for (double m = 0.0; m < 2.0; m += 0.001)	//computing production scale is slower for 3body decays!
	{						//because of many integrals involved
		nu.SetMass(m);

		double ps = nu.ProductionScale();

		//printing 1e-30 instead of 0 makes it nicer for log scale plot
		out << m << "\t" << (ps == 0 ? 1e-30 : ps) << std::endl;
	}

	return 0;
}

void Usage(char* argv0)
{
	std::cout << "Compute production scales" << std::endl;
	std::cout << "Usage : " << std::endl;
	std::cout << "decayplot [OPTIONS]" << std::endl;
	std::cout <<"\n  -o,  --output" << std::endl;
	std::cout << "\t\tOutput file to save plot" << std::endl;
	std::cout <<"\n  -c,  --channel" << std::endl;
	std::cout << "\t\tDecay channel" << std::endl;
	std::cout <<"\n  -E, -M, -T" << std::endl;
	std::cout << "\t\tChoose which mixing U2, giving also a magnitude" << std::endl;
	std::cout <<"\n  -h,  --help" << std::endl;
	std::cout << "\t\tPrint this message and exit" << std::endl;
}
