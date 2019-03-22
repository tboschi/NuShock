#include <iostream>
#include <fstream>
#include <cstring>
#include <sstream>
#include <vector>
#include <cmath>
#include <getopt.h>
#include <cstdlib>

#include "Tools.h"
#include "DecayRates.h"

double Time(double Mass, bool EM);
void Usage(char* argv0);
int main(int argc, char** argv)
{
	const struct option longopts[] = 
	{
		{"electron", 	no_argument,		0, 'E'},
		{"muon", 	no_argument,		0, 'M'},
		{"output", 	required_argument,	0, 'o'},
		{"help", 	no_argument,	 	0, 'h'},
		{0,	0, 	0,	0},
	};

	int index; 
	int iarg = 0;
	opterr = 1;
	
	//Initialize variables
	std::ofstream OutFile;
	bool EM = false;
	
	while((iarg = getopt_long(argc,argv, "EMo:h", longopts, &index)) != -1)
	{
		switch(iarg)
		{
			case 'E':
				EM = true;
				break;
			case 'M':
				EM = false;
				break;
			case 'o':
				OutFile.open(optarg);
				break;
			case 'h':
				Usage(argv[1]);
				break;
			default:
				return 1;
				break;
		}
	}

	std::ostream &Out = (OutFile.is_open()) ? OutFile : std::cout;

	Decay * TheGamma = new Decay();
	unsigned int Grid = 500;
	double Mass, Uu;

	for (double logMass = -2.0; logMass < -0.3; logMass += 1.7/Grid)	//increase mass log
	{
		Mass = pow(10.0, logMass);
		TheGamma->SetMass(Mass);

		double CosmoTime = Time(Mass, EM);

		for (double logUu2 = -15.0; logUu2 < -0.0; logUu2 += 15.0/Grid)	//increase Uu logarithmically
		{
			Uu = pow(10.0, 0.5*logUu2);
			if (EM)
				TheGamma->SetUe(Uu);
			else
				TheGamma->SetUm(Uu);
			double LifeTime = Const::fhBar / TheGamma->Total();

			if (LifeTime < CosmoTime)
			{
				Out << Mass << "\t" << Uu*Uu << "\t" << LifeTime << "\t" << CosmoTime << std::endl;
				break;
			}
		}
	}

	return 0;
}
	
double Time(double Mass, bool EM)
{
	/*
	if (Mass < Const::fMPion0)
	{
		if (EM)		//Electronic mixing
			return 1699 * pow(Mass*1000, -2.652) + 0.0544;
		else		//Muonic mixing
			return 1.287 * pow(Mass*1000, -1.828) + 0.04179;
	}
	else return 0.1;
	*/
	return 1;
}

void Usage(char* argv0)
{
	std::cout << "Description" << std::endl;
	std::cout << "Usage : " << std::endl;
	std::cout << argv0 << " [OPTIONS]" << std::endl;
	std::cout <<"\n  -i,  --input" << std::endl;
	std::cout << "\t\tInput file, tabulated as Mass\tUu\tnEvt" << std::endl;
	std::cout <<"\n  -o,  --output" << std::endl;
	std::cout << "\t\tOutput file" << std::endl;
	std::cout <<"\n  -t,  --threshold" << std::endl;
	std::cout << "\t\tEvent threshold for signal. As mass dependent threshold, t is the y-intercept" << std::endl;
	std::cout <<"\n  -m,  --massdep" << std::endl;
	std::cout << "\t\tIf specified, a mass dependance is considered for threshold, [thr] = t + m * [Mass]" << std::endl;
	std::cout <<"\n  -h,  --help" << std::endl;
	std::cout << "\t\tPrint this message and exit" << std::endl;
}
