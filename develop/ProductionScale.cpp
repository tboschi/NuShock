#include <iostream>
#include <fstream>
#include <vector>
#include <cstring>
#include <getopt.h>

#include "Tools.h"
#include "Physics.h"

int main(int argc, char** argv)
{

	const struct option longopts[] = 
	{
		{"output", 	required_argument, 	0, 'o'},
		{"help", 	no_argument,	 	0, 'h'},
		{0,	0, 	0,	0},
	};

	int index; 
	int iarg = 0;
	opterr = 1;
	
	std::string Channel;
	std::ofstream OutFile;
	
//Initialize variables
	
	while((iarg = getopt_long(argc,argv, "c:o:h", longopts, &index)) != -1)	
	{
		switch(iarg)
		{
			case 'c':
				Channel.assign(optarg);
				break;
			case 'o':
				OutFile.open(optarg);
				break;
			case 'h':
				std::cout << "Compute Shrock factors" << std::endl;
				std::cout << "Usage : " << std::endl;
				std::cout << "decayplot [OPTIONS]" << std::endl;
				std::cout <<"\n  -o,  --output" << std::endl;
				std::cout << "\t\tOutput file to save plot" << std::endl;
				std::cout <<"\n  -h,  --help" << std::endl;
				std::cout << "\t\tPrint this message and exit" << std::endl;
				return 1;
			default:
				return 1;
		}
	}



	//To have multiple output, handled by usage
	std::ostream &Out = (OutFile.is_open()) ? OutFile : std::cout;

	Neutrino *NuL = new Neutrino(0, Neutrino::Dirac | Neutrino::Left );
	Neutrino *NuR = new Neutrino(0, Neutrino::Dirac | Neutrino::Right);
	Neutrino *Nu0 = new Neutrino(0, Neutrino::Dirac | Neutrino::Unpolarised);

	NuL->SetMixings(1, 1, 1);
	NuR->SetMixings(1, 1, 1);
	Nu0->SetMixings(1, 1, 1);

	NuL->SetProductionChannel(Channel);
	NuR->SetProductionChannel(Channel);
	Nu0->SetProductionChannel(Channel);

	//Out << "#MS\tElPi\tElKa\tElCh\tMuPi\tMuKa\tMuCh" << std::endl;
	for (double t = 0.0; t < 0.35; t += 0.001)
	{
		NuL->SetMass(t);
		NuR->SetMass(t);
		Nu0->SetMass(t);

		double pL = NuL->ProductionScale();
		double pR = NuR->ProductionScale();
		double p0 = Nu0->ProductionScale();
		//std::cout << "Mass " << t << std::endl;
		Out << t;
		Out << "\t" << (pL == 0.0 ? 1e-30 : pL);
		Out << "\t" << (pR == 0.0 ? 1e-30 : pR);
		Out << "\t" << (p0 == 0.0 ? 1e-30 : p0);
		Out << std::endl;
	}

	//Main body

	//Garbage collection

	return 0;
}
