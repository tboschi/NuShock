#include <iostream>
#include <fstream>
#include <vector>
#include <cstring>
#include <getopt.h>

#include "physics/DecayRates.h"

void Usage(char* argv0);
int main(int argc, char** argv)
{

	const struct option longopts[] = 
	{
		{"dirac", 	no_argument,		0, 'r'},
		{"majorana", 	no_argument,		0, 'j'},
		{"output", 	required_argument, 	0, 'o'},
		{"channel", 	required_argument, 	0, 'c'},
		{"help", 	no_argument,	 	0, 'h'},
		{0,	0, 	0,	0},
	};

	int index; 
	int iarg = 0;
	opterr = 1;
	
	size_t ferm;
	std::string ferm_append;
	while((iarg = getopt_long(argc,argv, "c:o:rjh", longopts, &index)) != -1)	
	{
		switch(iarg)
		{
			case 'r':
				ferm = Neutrino::dirac;
				ferm_append = "_d";
				break;
			case 'j':
				ferm = Neutrino::majorana;
				ferm_append = "_m";
				break;
			case 'h':
			default:
				Usage(argv[1]);
				return 1;
		}
	}


	// unpolairsed?
	Neutrino nu(0.5, ferm | Neutrino::left);
	std::cout << "for neutrino " << nu << "\n";

	DecayRates hnl(nu);

	for (const auto &chan : Channel::Decays()) {
		std::cout << "BR for " << Channel::toString(chan)
			  << " is " << hnl.Branch(chan) << "\n";
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
