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



	//To have multiple output, handled by usage
	std::ostream &Out = (OutFile.is_open()) ? OutFile : std::cout;

	Neutrino *N_up = new Neutrino(0, Neutrino::Dirac | Neutrino::Particle | Neutrino::Left );
	Neutrino *N_do = new Neutrino(0, Neutrino::Dirac | Neutrino::Particle | Neutrino::Right);

	//Out << "#MS\tElPi\tElKa\tElCh\tMuPi\tMuKa\tMuCh" << std::endl;
	for (double t = 0.0; t < 1.0; t += 0.01)
	{
		N_up->SetMass(t);
		N_do->SetMass(t);
		Out << t << "\t";
		//Out << N_up->ProductionScale("KaonM") << "\t";
		//Out << N_do->ProductionScale("KaonM") << "\t";
		//Out << N_up->ProductionScale("MuonE") << "\t";
		//Out << N_do->ProductionScale("MuonE") << "\t";
		Out << N_up->ProductionScale("KaonCM") << "\t";
		//Out << N_do->ProductionScale("KaonCM") << "\t";
		Out << std::endl;
	}
		//Out << Kine::ShrockFactor(Pion, Elec, t) << "\t";
		//Out << Kine::ShrockFactor(Kaon, Elec, t) << "\t";
		//Out << Kine::ShrockFactor(Charm, Elec, t) << "\t";
		//Out << Kine::ShrockFactor(Pion, Muon, t) << "\t";
		//Out << Kine::ShrockFactor(Kaon, Muon, t) << "\t";
		//Out << Kine::ShrockFactor(Charm, Muon, t) << "\t";
		//Out << std::endl;
	//Create object from classes
	//e.g.	Decay * SuperGamma = new Decay(M_Sterile, U_e, U_m, U_t);

	//Main body

	//Garbage collection

	return 0;
}
