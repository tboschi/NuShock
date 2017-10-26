#include <iostream>
#include <fstream>
#include <cstring>
#include <sstream>
#include <vector>
#include <cmath>
#include <getopt.h>
#include <cstdlib>

void Usage(char* argv0);
int main(int argc, char** argv)
{
	const struct option longopts[] = 
	{
		{"threshold", 	required_argument,	0, 't'},
		{"massdep", 	required_argument,	0, 'm'},
		{"start", 	required_argument,	0, 'A'},
		{"end", 	required_argument,	0, 'B'},
		{"input", 	no_argument,		0, 'i'},
		{"output", 	required_argument,	0, 'o'},
		{"help", 	no_argument,	 	0, 'h'},
		{0,	0, 	0,	0},
	};

	int index; 
	int iarg = 0;
	opterr = 1;
	
	//Initialize variables
	std::ofstream OutFile;
	std::ifstream InAFile, InBFile;
	
	while((iarg = getopt_long(argc,argv, "A:B:o:h", longopts, &index)) != -1)
	{
		switch(iarg)
		{
			case 'A':
				InAFile.open(optarg);
				break;
			case 'B':
				InBFile.open(optarg);
				break;
			case 'o':
				OutFile.open(optarg);
				break;
			case 'h':
				Usage(argv[0]);
				return 1;
			default:
				return 1;
				break;
		}
	}

	std::ostream &Out = (OutFile.is_open()) ? OutFile : std::cout;

	double Mass = 0, Uu = 0, nEvt = 0;
	std::string Line, Key, Name;
	std::stringstream ssL;
	std::vector<double> vnA, vnB, vMass, vUu;

	while (std::getline(InAFile, Line))
	{
		if (Line[0] == '#') continue;

		ssL.str("");
		ssL.clear();
		ssL << Line;

		ssL >> Mass >> Uu >> nEvt;
		vnA.push_back(nEvt);
	}	
	InAFile.close();

	while (std::getline(InBFile, Line))
	{
		if (Line[0] == '#') continue;

		ssL.str("");
		ssL.clear();
		ssL << Line;

		ssL >> Mass >> Uu >> nEvt;
		vMass.push_back(Mass);
		vUu.push_back(Uu);
		vnB.push_back(nEvt);
	}	
	InBFile.close();

	for (unsigned int i = 0; i < vMass.size(); ++i)
	{
		Out << vMass.at(i) << "\t";
		Out << vUu.at(i) << "\t";
		Out << vnA.at(i)+vnB.at(i) << std::endl;
	}

	return 0;
}
	
void Usage(char* argv0)
{
	std::cout << "Description" << std::endl;
	std::cout << "Usage : " << std::endl;
	std::cout << argv0 << " [OPTIONS]" << std::endl;
	std::cout <<"\n  -A,  --input" << std::endl;
	std::cout << "\t\tInput file A, tabulated as Mass\tUu\tnEvt" << std::endl;
	std::cout <<"\n  -B,  --input" << std::endl;
	std::cout << "\t\tInput file B, tabulated as Mass\tUu\tnEvt" << std::endl;
	std::cout <<"\n  -o,  --output" << std::endl;
	std::cout << "\t\tOutput file" << std::endl;
	std::cout <<"\n  -h,  --help" << std::endl;
	std::cout << "\t\tPrint this message and exit" << std::endl;
}
