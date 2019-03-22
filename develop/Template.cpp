#include <iostream>
#include <fstream>
#include <vector>
#include <cstring>
#include <getopt.h>

void Usage(char* Name);
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
	
	while((iarg = getopt_long(argc,argv, "1:2h", longopts, &index)) != -1)
	{
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
				Usage(argv[0]);
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

	std::cout << "Hello World!" << std::endl;
	//Garbage collection

	return 0;
}
	
void Usage(char *Name)
{
	std::cout << "Description" << std::endl;
	std::cout << "Usage : " << std::endl;
	std::cout << Name << " [OPTIONS]" << std::endl;
	std::cout <<"\n  -1,  --option1" << std::endl;
	std::cout << "\t\tDescription" << std::endl;
	std::cout <<"\n  -2,  --option2" << std::endl;
	std::cout << "\t\tDescription" << std::endl;
	std::cout <<"\n  -h,  --help" << std::endl;
	std::cout << "\t\tPrint this message and exit" << std::endl;
}
