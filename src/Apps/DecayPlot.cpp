#include <iostream>
#include <fstream>

#include "Tools.h"
#include "DecayRates.h"

int main(int argc, char** argv)
{

	const struct option longopts[] = 
	{
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
	
	double M_Sterile = 0.0350;	//350 MeV
	double U_e = 1.0/sqrt(3.0);		//All mixing enabled to be maximal
	double U_m = 1.0/sqrt(3.0);
	double U_t = 1.0/sqrt(3.0);
	
	while(iarg != -1)
	{
		iarg = getopt_long(argc,argv, "s:e:m:t:", longopts, &index);

		switch(iarg)
		{
			case 's':
				mS = strtod(optarg,NULL);
				break;
			case 'e':
				mS = strtod(optarg,NULL);
				break;
			case 'm':
				mS = strtod(optarg,NULL);
				break;
			case 't':
				mS = strtod(optarg,NULL);
				break;
			case 'h':
				std::cout << "Usage : " << std::endl;
				std::cout << "decayplot [OPTIONS]" std::endl;
				std::cout <<"\n  -s,  --sterile" << std::endl;
				std::cout << "\t\tSpeificy sterile mass in GeV, default 0.350 MeV" << std::endl;
				std::cout <<"\n  -e,  --ue4" << std::endl;
				std::cout << "\t\tSpecifiy electron-heavy mixing, default 1/sqrt(3)" << std::endl;
				std::cout <<"\n  -m,  --um4" << std::endl;
				std::cout << "\t\tSpecifiy muon-heavy mixing, default 1/sqrt(3)" << std::endl;
				std::cout <<"\n  -t,  --ut4" << std::endl;
				std::cout << "\t\tSpecifiy tau-heavy mixing, default 1/sqrt(3)" << std::endl;
				std::cout <<"\n  -h,  --help" << std::endl;
				std::cout << "\t\tPrint this message and exit" << std::endl;
				break;
			default:
				break;
		return 0;
		}
	
	}

	DecayRates * SuperGamma = new SuperGamma(M_Sterile, U_e, U_m, U_t);
	std::vector<std::string> vChannels = SuperGamma->ListChannels();

	for (int i = 0; i < vChannels.size(); ++i)
		std::cout << vChannel.at(i) << std::endl;

	delete SuperGamma;

	return 0;
}
	
