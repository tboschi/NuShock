#include <iostream>
#include <fstream>
#include <vector>
#include <cstring>
#include <getopt.h>

#include "Tools.h"
#include "EventGenerator.h"
#include "3Body.h"

int main(int argc, char** argv)
{
	const struct option longopts[] = 
	{
		{"output", 	required_argument, 	0, 'o'},
		{"sterile", 	required_argument, 	0, 's'},
		{"ue4", 	required_argument, 	0, 'e'},
		{"um4", 	required_argument, 	0, 'm'},
		{"ut4", 	required_argument, 	0, 't'},
		{"help", 	no_argument,		0, 'h'},
		{0,	0, 	0,	0},
	};

	int index; 
	int iarg = 0;
	opterr = 1;
	
	std::ofstream OutFile;
	double M_Sterile = 0.300;	//350 MeV
	double U_e = 1.0/sqrt(3.0);	//All mixing enabled to be maximal
	double U_m = 1.0/sqrt(3.0);
	double U_t = 1.0/sqrt(3.0);
	
	while((iarg = getopt_long(argc,argv, "o:s:e:m:t:h", longopts, &index)) != -1)	
	{
		switch(iarg)
		{
			case 'o':
				OutFile.open(optarg);
				if (!OutFile.is_open())
				{
					std::cout << "Wrong input!" << std::endl;
					return 1;
				}
				break;
			case 'e':
				U_e = strtod(optarg,NULL);
				break;
			case 'm':
				U_m = strtod(optarg,NULL);
				break;
			case 't':
				U_t = strtod(optarg,NULL);
				break;
			case 'h':
				std::cout << "Compute decay spectrum from 0 MeV to 500 MeV" << std::endl;
				std::cout << "Usage : " << std::endl;
				std::cout << "decayplot [OPTIONS]" << std::endl;
				std::cout <<"\n  -o,  --output" << std::endl;
				std::cout << "\t\tOutput file to save plot; if not specified the standard output is used instead" << std::endl;
				std::cout <<"\n  -e,  --ue4" << std::endl;
				std::cout << "\t\tSpecifiy electron-heavy mixing, default 1/sqrt(3)" << std::endl;
				std::cout <<"\n  -m,  --um4" << std::endl;
				std::cout << "\t\tSpecifiy muon-heavy mixing, default 1/sqrt(3)" << std::endl;
				std::cout <<"\n  -t,  --ut4" << std::endl;
				std::cout << "\t\tSpecifiy tau-heavy mixing, default 1/sqrt(3)" << std::endl;
				std::cout <<"\n  -h,  --help" << std::endl;
				std::cout << "\t\tPrint this message and exit" << std::endl;
				return 1;
			default:
				break;
		}
	}

	std::ostream &Out = (OutFile.is_open()) ? OutFile : std::cout;

	Decay * SuperGamma = new Decay(M_Sterile, U_e, U_m, U_t);

	std::vector<std::string> vChannel(SuperGamma->ListChannels());
	
	//ThreeBody * SuperGamma = new ThreeBody("ALL", M_Sterile, U_e, U_m, U_t);
/*
	double xa, xb, dx;
	dx = SuperGamma->xLim(xa, xb);
	double MaxGamma = SuperGamma->MaxGamma();
	for (double x = xa; x < xb; x += dx/100)
	{
		SuperGamma->SetX(x);
		double ya, yb, dy;
		dy = SuperGamma->yLim(ya, yb);
		for (double y = ya; y < yb; y += dy/100)
		{
			SuperGamma->SetY(y);
			double GG = SuperGamma->ddGamma();
			Out << x << "\t" << y << "\t" << GG/MaxGamma << std::endl; 
		}
	}
*/
	
	Out << "#Mass\t";
	for (int i = 1; i < vChannel.size(); ++i)
		Out << vChannel.at(i) << "\t";
	Out << std::endl;

	for (M_Sterile = 0; M_Sterile < 0.500; M_Sterile += 0.001)
	{
		SuperGamma->SetMass(M_Sterile);

		Out << M_Sterile << "\t";
		for (int i = 1; i < vChannel.size(); ++i)
			Out << SuperGamma->Branch(vChannel.at(i)) << "\t";
		Out << std::endl;
	}
	

	return 0;
}
