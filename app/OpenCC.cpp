#include <iostream>
#include <fstream>
#include <cstring>
#include <vector>
#include <getopt.h>

#include "physics/OpenQQ.h"

void Usage(char* argv0);
int main(int argc, char** argv)
{
	const struct option longopts[] = 
	{
		{"", 	required_argument, 	0, 'b'},
		{"confidence", 	required_argument, 	0, 'C'},
		{"help", 	no_argument,	 	0, 'h'},
		{0,	0, 	0,	0},
	};

	int index; 
	int iarg = 0;
	opterr = 1;

	std::ofstream outfile;

	while((iarg = getopt_long(argc,argv, "o:h", longopts, &index)) != -1)
	{
		switch(iarg)
		{
			case 'o':
				outfile.open(optarg);
				break;
			case 'h':
				Usage(argv[0]);
				return 1;
			default:
				break;
		}
	}

	std::ostream &out = (outfile.is_open()) ? outfile : std::cout;

	OpenQQ *xsec = new OpenCC(argv[optind]);

	// scan over center of mass energy in GeV
	for (double cme = 40; cmE < 450; cmE += 40) {
		xsec->SetCMEnergy(cme);
		double err, chi2;
		out << "Center of mass energy " << cme << " GeV:\t"
		    << xsec->Calculate(err, chi2) << " ± " << std::sqrt(chi2) << std::endl;
	}
	
	if (outfile.is_open())
		outfile.close();

	return 0;
}

void Usage(char* argv0)
{
	std::cout << "usage: " << argv0 << " <card> [-o <output>]\n\n"
		  << "Computes the open charm cross sections on target specified in\n"
		  << "the <card> file passed. Optionally saves output to a file <output>\n\n";
}
