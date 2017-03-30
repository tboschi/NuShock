#include <iostream>
#include <fstream>
#include <vector>
#include <cstring>
#include <getopt.h>

#include "headers.h"

int main(int argc, char** argv)
{

	const struct option longopts[] = 
	{
		{"option1", 	required_argument, 	0, '1'},
		{"option2", 	no_argument,	 	0, '2'},
		{"help", 	no_argument,	 	0, 'h'},
		{0,	0, 	0,	0},
	};

	int index; 
	int iarg = 0;
	opterr = 1;
	
	std::ofstream OutFile;
	
//Initialize variables
	double M_Sterile = 0.350;	//350 MeV
	double U_e = 1.0/sqrt(3.0);	//All mixing enabled to be maximal
	double U_m = 1.0/sqrt(3.0);
	double U_t = 1.0/sqrt(3.0);
	
	while(iarg != -1)
	{
		iarg = getopt_long(argc,argv, "1:2h", longopts, &index);

		switch(iarg)
		{
			case '1':
				//do something and use optarg
				//e.g. U_e = strtod(optarg,NULL);
				break;
			case '2':
				//do something
				break;
			case 'h':
				std::cout << "Description" << std::endl;
				std::cout << "Usage : " << std::endl;
				std::cout << "name [OPTIONS]" << std::endl;
				std::cout <<"\n  -1,  --option1" << std::endl;
				std::cout << "\t\tDescription" << std::endl;
				std::cout <<"\n  -2,  --option2" << std::endl;
				std::cout << "\t\tDescription" << std::endl;
				std::cout <<"\n  -h,  --help" << std::endl;
				std::cout << "\t\tPrint this message and exit" << std::endl;
				return 1;
			default:
				break;
		}
	
	}

	//To have multiple output, handled by usage
	std::ostream &Out = (OutFile.is_open()) ? OutFile : std::cout;

	//Create object from classes
	//e.g.	Decay * SuperGamma = new Decay(M_Sterile, U_e, U_m, U_t);

	//Main body

	//Garbage collection

	return 0;
}
	
