#include <iostream>
#include <fstream>
#include <vector>
#include <cstring>
#include <getopt.h>

#include "physics/Production.h"

void Usage(char* argv0);
int main(int argc, char** argv)
{

	const struct option longopts[] = 
	{
		{"left", 	no_argument,		0, 'L'},
		{"right", 	no_argument,		0, 'R'},
		{"help", 	no_argument,	 	0, 'h'},
		{0,	0, 	0,	0},
	};

	int index; 
	int iarg = 0;
	opterr = 1;
	
	size_t ferm = Neutrino::majorana | Neutrino::antiparticle;
	while((iarg = getopt_long(argc,argv, "RLh", longopts, &index)) != -1)	
	{
		switch(iarg)
		{
			case 'R':
				ferm = ferm | Neutrino::right;
				break;
			case 'L':
				ferm = ferm | Neutrino::left;
				break;
			case 'h':
			default:
				Usage(argv[1]);
				return 1;
		}
	}


	Neutrino nu(0.1, Neutrino::dirac | Neutrino::unpolarized);

	Mixing mix(1., 1., 1.);

	Production hnl(nu);
	std::cout << "PionM " << hnl.Scale(Channel::PionM, mix) << "\n";
	std::cout << "gamma " << hnl.Gamma(Channel::PionM, mix) << "\n";
	std::cout << "KaonM " << hnl.Scale(Channel::KaonM, mix) << "\n";
	std::cout << "gamma " << hnl.Gamma(Channel::KaonM, mix) << "\n";
	std::cout << "KaonCM " << hnl.Scale(Channel::KaonCM, mix) << "\n";
	std::cout << "gamma " << hnl.Gamma(Channel::KaonCM, mix) << "\n";
	std::cout << "Kaon0M " << hnl.Scale(Channel::Kaon0M, mix) << "\n";
	std::cout << "gamma " << hnl.Gamma(Channel::Kaon0M, mix) << "\n";
	std::cout << "MuonM " << hnl.Scale(Channel::MuonM, mix) << "\n";
	std::cout << "gamma " << hnl.Gamma(Channel::MuonM, mix) << "\n";
	std::cout << "CharM " << hnl.Scale(Channel::CharmM, mix) << "\n";
	std::cout << "gamma " << hnl.Gamma(Channel::CharmM, mix) << "\n";

	//std::cout << " M2   " << hnl.M2_MesonThree(0.5, 0.5, 0.1, 0.2, 0.3, Const::KCL_, Const::KCL0) << "\n";
	return 1;

	for (const auto &chan : Channel::Productions()) {
		std::cout << "Scale for " << Channel::toString(chan)
			  << " is " << hnl.Scale(chan) << "\n";
	}

	return 0;
}

void Usage(char* argv0)
{
	std::cout << "Compute decay widths/branching ratios" << std::endl;
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
