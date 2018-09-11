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

	Production *TheP = new Production();
	Production *The0 = new Production();
	double *mix = new double[3];
	bool Ferm = true;
	int Hel = 0;
	mix[0] = 1.0;
	mix[1] = 1.0;
	mix[2] = 1.0;
	The0->SetNeutrino(0.0, mix, Ferm, 1, Hel);

	TheP->SetMass(1);
	TheP->SetNeutrino(0.0, mix, Ferm, 1, Hel);
	std::cout << 1.0/TheP->I_MesonThree(0 ,0, 0, 0, 0) << std::endl;
	//std::cout << 1.0/TheP->MesonThreeDecay2(1 ,0, 0, 0, 0) << std::endl;
	//std::cout << Const::fPi3*TheP->dGammad2_3B(1) << std::endl;
	//for (double t = 0.0; t < 0.35; t += 0.001)
	//{
	//	TheP->SetNeutrino(t, mix, Ferm, 1, Hel);

	//	Out << t;
	//	//Out << "\t" << TheP->Gamma(Amplitude::_Kaon0E)/The0->Gamma(Amplitude::_Kaon0E);
	//	//Out << "\t" << TheP->MesonThreeDecay2(Const::fMKaon0, Const::fMPion, Const::fMElectron, Const::fK0L_, Const::fK0L0)/The0->MesonThreeDecay2(Const::fMKaon0, Const::fMPion, Const::fMElectron, Const::fK0L_, Const::fK0L0);
	//	Out << "\t" << TheP->Gamma(Amplitude::_Kaon0E);
	//	Out << "\t" << TheP->MesonThreeDecay2(Const::fMKaon0, Const::fMPion, Const::fMElectron, Const::fK0L_, Const::fK0L0);
	//	Out << std::endl;
	//}

	//Main body

	//Garbage collection

	return 0;
}
