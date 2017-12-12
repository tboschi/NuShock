#include <iostream>
#include <fstream>
#include <cstring>
#include <sstream>
#include <vector>
#include <cmath>
#include <getopt.h>
#include <cstdlib>
#include <set>

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
		{"background", 	required_argument,	0, 'W'},
		{"output", 	required_argument,	0, 'o'},
		{"help", 	no_argument,	 	0, 'h'},
		{0,	0, 	0,	0},
	};

	int index; 
	int iarg = 0;
	opterr = 1;
	
	//Initialize variables
	std::ofstream OutFile;
	std::ifstream InFile;
	double Threshold = 2.44;
	double MassDep = 0.0;
	double MassA = 0.0, MassB = 0.5;
	
	while((iarg = getopt_long(argc,argv, "i:o:t:m:A:B:h", longopts, &index)) != -1)
	{
		switch(iarg)
		{
			case 'i':
				InFile.open(optarg);
				break;
			case 'o':
				OutFile.open(optarg);
				break;
			case 't':
				Threshold = strtod(optarg, NULL);
				break;
			case 'm':
				MassDep = strtod(optarg, NULL);
				break;
			case 'A':
				MassA = strtod(optarg, NULL);
				break;
			case 'B':
				MassB = strtod(optarg, NULL);
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

	double rMass, rUu2, rEvt;
	double Uu2Min = 2.0, Uu2Max = 0.0;
	std::vector <double> vMass, vUu2Min, vUu2Max;

	std::string Line;
	std::stringstream ssL;
	while (std::getline(InFile, Line))
	{
		if (Line[0] == '#') continue;

		ssL.str("");
		ssL.clear();
		ssL << Line;

		ssL >> rMass >> rUu2 >> rEvt;

		if (rUu2 < Uu2Min)
			Uu2Min = rUu2;
		if (rUu2 > Uu2Max)
			Uu2Max = rUu2;

		if (vMass.empty())
			vMass.push_back(rMass);
		else if (rMass != vMass.back())
		{	
			vMass.push_back(rMass);
			vUu2Min.push_back(Uu2Min);
			vUu2Max.push_back(Uu2Max);
			Uu2Min = 2.0;
			Uu2Max = 0.0;
		}
	}
	//insertion at EOF
	vUu2Min.push_back(Uu2Min);
	vUu2Max.push_back(Uu2Max);

	std::cout << "size check " << vMass.size() << "\t" << vUu2Min.size() << "\t" << vUu2Max.size() << std::endl;
	//lower line
	for (unsigned int i = 0; i < vMass.size(); ++i)
		Out << vMass.at(i) << "\t" << vUu2Min.at(i) << std::endl;

	//upper line
	for (unsigned int i = vMass.size()-1; i > 0; --i)
		Out << vMass.at(i) << "\t" << vUu2Max.at(i) << std::endl;

	//close the line
	Out << vMass.front() << "\t" << vUu2Min.front() << std::endl;

	InFile.close();
	OutFile.close();

	return 0;
}

void Usage(char* argv0)
{
	std::cout << "Description" << std::endl;
	std::cout << "Usage : " << std::endl;
	std::cout << argv0 << " [OPTIONS]" << std::endl;
	std::cout <<"\n  -i,  --input" << std::endl;
	std::cout << "\t\tInput file, tabulated as Mass\tUu\tnEvt" << std::endl;
	std::cout <<"\n  -o,  --output" << std::endl;
	std::cout << "\t\tOutput file" << std::endl;
	std::cout <<"\n  -t,  --threshold" << std::endl;
	std::cout << "\t\tEvent threshold for signal. As mass dependent threshold, t is the y-intercept" << std::endl;
	std::cout <<"\n  -m,  --massdep" << std::endl;
	std::cout << "\t\tIf specified, a mass dependance is considered for threshold, [thr] = t + m * [Mass]" << std::endl;
	std::cout <<"\n  -h,  --help" << std::endl;
	std::cout << "\t\tPrint this message and exit" << std::endl;
}
