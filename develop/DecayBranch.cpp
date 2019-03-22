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
	bool UeFlag = false, 
	     UmFlag = false, 
	     UtFlag = false;

	while((iarg = getopt_long(argc,argv, "c:o:EMTh", longopts, &index)) != -1)	
	{
		switch(iarg)
		{
			case 'c':
				Channel.assign(optarg);
				break;
			case 'o':
				OutFile.open(optarg);
				break;
			case 'E':
				UeFlag = true;
				break;
			case 'M':
				UmFlag = true;
				break;
			case 'T':
				UtFlag = true;
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

	Neutrino *NuD = new Neutrino(0, Neutrino::Dirac    | Neutrino::Unpolarised );
	//Neutrino *NuM = new Neutrino(0, Neutrino::Majorana | Neutrino::Unpolarised );

	NuD->SetDecayChannel(Channel);
	//NuM->SetDecayChannel(Channel);

	if (UeFlag)
	{
		NuD->SetMixings(1, NuD->Um(), NuD->Ut());
		//NuM->SetMixings(1, NuM->Um(), NuM->Ut());
	}
	if (UmFlag)
	{
		NuD->SetMixings(NuD->Ue(), 1, NuD->Ut());
		//NuM->SetMixings(NuM->Ue(), 1, NuM->Ut());
	}
	if (UtFlag)
	{
		NuD->SetMixings(NuD->Ue(), NuD->Um(), 1);
		//NuM->SetMixings(NuM->Ue(), NuM->Um(), 1);
	}

	Out << "#Mass\tDirac\tMajorana\n";

	for (double t = 0.0; t < 2.0; t += 0.0001)
	{
		NuD->SetMass(t);
		//NuM->SetMass(t);
		//std::cout << "mass " << t << std::endl;

		double db = NuD->DecayBranch();
		Out << t << "\t" << (db == 0 ? 1e-30 : db) << std::endl;//"\t" << NuM->DecayBranch() << std::endl;
	}

	return 0;
}
