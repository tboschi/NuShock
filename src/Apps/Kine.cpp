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



	double Elec = Const::fMElectron;
	double Muon = Const::fMMuon;
	double Pion = Const::fMPion;
	double Pion0 = Const::fMPion0;
	double Kaon = Const::fMKaon;
	double Kaon0 = Const::fMKaon0;
	//To have multiple output, handled by usage
	std::ostream &Out = (OutFile.is_open()) ? OutFile : std::cout;

	Out << "#MS\tMuPi\tMuKa\tMuKaPi0\tMuKa0Pi\t";
	Out << "#ElPi\tElKa\tElKaPi0\tElKa0Pi" << std::endl;
	for (double t = 0; t < 0.5; t += 0.5/1000)
	{
		Out << t << "\t";
		Out << Kine::ShrockFactor(Pion, Muon, t) << "\t";
		Out << Kine::ShrockFactor(Kaon, Muon, t) << "\t";
		//Out << Kine::ShrockFactor(Kaon-Pion0, Muon, t) << "\t";
		//Out << Kine::ShrockFactor(Kaon0-Pion, Muon, t) << "\t";
		Out << Kine::ShrockFactor(Pion, Elec, t) << "\t";
		Out << Kine::ShrockFactor(Kaon, Elec, t) << "\t";
		//Out << Kine::ShrockFactor(Kaon-Pion0, Elec, t) << "\t";
		//Out << Kine::ShrockFactor(Kaon0-Pion, Elec, t) << "\t";
		Out << std::endl;
	}
	//Create object from classes
	//e.g.	Decay * SuperGamma = new Decay(M_Sterile, U_e, U_m, U_t);

	//Main body

	//Garbage collection

	return 0;
}
