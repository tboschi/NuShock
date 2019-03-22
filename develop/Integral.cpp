#include <iostream>
#include <fstream>
#include <vector>
#include <cstring>
#include <getopt.h>

#include "Tools.h"

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
	
	std::ofstream OutFile;
	
//Initialize variables
	
	while((iarg = getopt_long(argc,argv, "o:h", longopts, &index)) != -1)	
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
				return 1;
		}
	}



	double M_Electron = Const::fMElectron;	//350 MeV
	//To have multiple output, handled by usage
	std::ostream &Out = (OutFile.is_open()) ? OutFile : std::cout;

	std::cout << "ratio " << M_Electron/0.35 << std::endl;
	std::cout << "0.350 -> " << Kine::ShrockLambda(0.5,0.0,M_Electron/0.35) << std::endl;
	std::cout << "0.350 -> " << Kine::ShrockLambda(1.0,0.5,M_Electron/0.35) << std::endl;
	std::cout << "0.350 -> " << Kine::I1_f(0.5,0,M_Electron/0.35,M_Electron/0.35) << std::endl;
	std::cout << "0.350 -> " << Kine::I1_xyz(0,M_Electron/0.35,M_Electron/0.35) << std::endl;
	Out << "#MS\tIntegral\n";
	for (double t = 0; t < 0.5; t += 0.5/1000)
	{
		Out << t << "\t";
		Out << Kine::I1_xyz(0, M_Electron/t, M_Electron/t) << std::endl;
	}
	//Create object from classes
	//e.g.	Decay * SuperGamma = new Decay(M_Sterile, U_e, U_m, U_t);

	//Main body

	//Garbage collection

	return 0;
}
