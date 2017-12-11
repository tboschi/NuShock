#include <iostream>
#include <fstream>
#include <vector>
#include <cstring>
#include <algorithm>
#include <getopt.h>

#include "Tools.h"
#include "Nucleon.h"
#include "TLorentzVector.h"
#include "cuba.h"

void Usage(char* argv0);
int main(int argc, char** argv)
{
	const struct option longopts[] = 
	{
		{"neutrino", 	no_argument,	 	0, 'n'},
		{"antineutrino",no_argument,	 	0, 'a'},
		{"neutron", 	no_argument,	 	0, 'N'},
		{"proton", 	no_argument,	 	0, 'P'},
		{"mass", 	required_argument, 	0, 'M'},
		{"energy", 	required_argument, 	0, 'E'},
		{"output", 	required_argument, 	0, 'o'},
		{"vegas", 	no_argument,	 	0, 'V'},
		{"help", 	no_argument,	 	0, 'h'},
		{0,	0, 	0,	0},
	};

	int index; 
	int iarg = 0;
	opterr = 1;
	
	std::ofstream OutFile;
	double Mass = 0.0;	//sterile neutrino mass (GeV)
	double Energy = 20;
	bool Nucl = false;	//T = proton, F = neutron
	bool Neut = true;	//T = neutrino, F = antineutrino
	bool Vegas = false;	//MC evaluation instead of numeric sum
	
//Initialize variables
	
	while((iarg = getopt_long(argc,argv, "naNPE:M:o:Vh", longopts, &index)) != -1)	
	{
		switch(iarg)
		{
			case 'n':
				Neut = true;
				break;
			case 'a':
				Neut = false;
				break;
			case 'N':
				Nucl = false;
				break;
			case 'P':
				Nucl = true;
				break;
			case 'E':
				Energy = strtod(optarg, NULL);
				break;
			case 'M':
				Mass = strtod(optarg, NULL);
				break;
			case 'o':
				OutFile.open(optarg);
				break;
			case 'V':
				Vegas = true;
				break;
			case 'h':
				Usage(argv[0]);
				return 1;
			default:
				return 1;
		}
	}

	//To have multiple output, handled by usage
	std::ostream &Out = (OutFile.is_open()) ? OutFile : std::cout;

	Nucleon *NN = new Nucleon(Neut, Nucl); 

	//TLorentzVector Ns(0, 0, Mass, Mass);
	//Nucl->SetSterile(Ns);

	TLorentzVector Nu(0, 0, Energy, Energy);
	NN->SetProbe(Nu);
	NN->SetPS(Mass);
	for (unsigned int i = 0; i < 1000; ++i)
	{
		NN->GeneratePS();
		Out << NN->Q2() << "\t" << NN->dSigmadQ2() << std::endl;
	}

	return 0;
}

void Usage(char* argv0)
{
	std::cout << "Cross section for neutrino - free nucleon scattering" << std::endl;
	std::cout << "Usage : " << std::endl;
	std::cout << argv0 << " [OPTIONS]" << std::endl;
	std::cout <<"\n  -N,  --neutron" << std::endl;
	std::cout << "\t\tSet the target as a neutron" << std::endl;
	std::cout <<"\n  -P,  --proton" << std::endl;
	std::cout << "\t\tSet the target as a proton" << std::endl;
	std::cout <<"\n  -M,  --mass" << std::endl;
	std::cout << "\t\tSet heavy neutrino mass (GeV)" << std::endl;
	std::cout <<"\n  -U,  --coupling" << std::endl;
	std::cout << "\t\tSet coupling U_a4" << std::endl;
	std::cout <<"\n  -E,  --energy" << std::endl;
	std::cout << "\t\tSet neutrino probe energy" << std::endl;
	std::cout <<"\n  -o,  --output" << std::endl;
	std::cout << "\t\tOutput file to save" << std::endl;
	std::cout <<"\n  -h,  --help" << std::endl;
	std::cout << "\t\tPrint this message and exit" << std::endl;
}
