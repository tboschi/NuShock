#include <iostream>
#include <fstream>
#include <vector>
#include <cstring>
#include <sstream>
#include <getopt.h>

#include "Tools.h"
#include "3Body.h"

int main(int argc, char** argv)
{
	const struct option longopts[] = 
	{
		{"output", 	required_argument, 	0, 'o'},
		{"smconfig", 	required_argument, 	0, 's'},
		{"muon", 	no_argument,	 	0, 'm'},
		{"kaon", 	no_argument, 		0, 'k'},
		{"kaon0", 	no_argument,	 	0, 'K'},
		{"help", 	no_argument,		0, 'h'},
		{0,	0, 	0,	0},
	};

	int index; 
	int iarg = 0;
	opterr = 1;
	
	std::ifstream ConfigFile;
	std::ofstream OutFile;
	bool IsElectron = false;
	bool IsMuon = false;
	bool MuonFlag = false;
	bool KaonFlag = false;
	bool Kaon0Flag = false;
	
	while((iarg = getopt_long(argc,argv, "o:s:EMmkKh", longopts, &index)) != -1)	
	{
		switch(iarg)
		{
			case 'o':
				OutFile.open(optarg);
				break;
			case 's':
				ConfigFile.open(optarg);
				break;
			case 'E':
				IsElectron = true;
				break;
			case 'M':
				IsMuon = true;
				break;
			case 'm':
				MuonFlag = true;
				break;
			case 'k':
				KaonFlag = true;
				break;
			case 'K':
				Kaon0Flag = true;
				break;
			case 'u':
				IsMuon = true;
				break;
			case 'h':
				std::cout << "Compute decay width fro three body decay" << std::endl;
				std::cout << "Usage : " << std::endl;
				std::cout << "decayplot [OPTIONS]" << std::endl;
				std::cout <<"\n  -s,  --smconfig" << std::endl;
				std::cout << "\t\tInput file, with the model parameter configuration" << std::endl;
				std::cout <<"\n  -o,  --output" << std::endl;
				std::cout << "\t\tOutput file to save plot; if not specified the standard output is used instead" << std::endl;
				std::cout <<"\n  -E" << std::endl;
				std::cout << "\t\tSpecifiy electron channel mode for kaon decay" << std::endl;
				std::cout <<"\n  -M" << std::endl;
				std::cout << "\t\tSpecifiy muon channel mode for kaon decay" << std::endl;
				std::cout <<"\n  -m,  --muon" << std::endl;
				std::cout << "\t\tMuon Michel decay" << std::endl;
				std::cout <<"\n  -k,  --kaon" << std::endl;
				std::cout << "\t\tKaon semileptonic decay" << std::endl;
				std::cout <<"\n  -K,  --kaon0" << std::endl;
				std::cout << "\t\tKaon0 semileptonic decay" << std::endl;
				std::cout <<"\n  -h,  --help" << std::endl;
				std::cout << "\t\tPrint this message and exit" << std::endl;
				return 1;
			default:
				break;
		}
	}

	std::ostream &Out = (OutFile.is_open()) ? OutFile : std::cout;

	double M_MSterile; 
	double U_e;
	double U_m;
	double U_t;

	std::string Line, Key, Name;
	std::stringstream ssL;
	double Element;

	while (std::getline(ConfigFile, Line))
	{
		if (Line[0] == '#') continue;

		ssL.str("");
		ssL.clear();
		ssL << Line;

		ssL >> Key >> Element;
		if (Key == "M_Sterile") M_MSterile = Element;
		if (Key == "U_e") U_e = Element;
		if (Key == "U_m") U_m = Element;
		if (Key == "U_t") U_t = Element;
	}
	ConfigFile.close();

	ThreeBody *Decay;
	if (MuonFlag)
		Decay = new ThreeBody("Muon", M_MSterile, U_e, U_m, U_t);
	else if (KaonFlag)
		Decay = new ThreeBody("Kaon", M_MSterile, U_e, U_m, U_t);
	else if (Kaon0Flag)
		Decay = new ThreeBody("Kaon0", M_MSterile, U_e, U_m, U_t);
	else Decay = 0;

	Decay->SetX(0.5);
	if (Decay->InLimX())
		std::cout << "Muon " << "\t" << Decay->M2MuonIntY() << std::endl;
	std::cout << Decay->GetEnergyX() << std::endl;
	std::cout << Const::fGF2 << std::endl;

	return 0;
}
