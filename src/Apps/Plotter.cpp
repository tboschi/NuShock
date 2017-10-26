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

	std::vector<double> vMass, vUu;
	for (double logMass = -2.0; logMass < -0.3; logMass += 0.0068)	//increase mass log
		vMass.push_back(pow(10.0, logMass));
	for (double logUu2 = -10.0; logUu2 < -0.0; logUu2 += 0.04)	//increase Uu logarithmically
		vUu.push_back(pow(10.0, logUu2));

	std::vector<std::vector<double> > mnEvt;
	std::vector<double> vnEvt;

	double Mass_ = -1, Uu_ = -1, nEvt_ = -1;
	double Mass = 0, Uu = 0, nEvt = 0;
	std::string Line, Key, Name;
	std::stringstream ssL;
	while (std::getline(InFile, Line))
	{
		if (Line[0] == '#') continue;

		ssL.str("");
		ssL.clear();
		ssL << Line;

		ssL >> Mass >> Uu >> nEvt;

		if (Mass != Mass_ && Mass_ >= 0)
		{
			mnEvt.push_back(vnEvt);
			vnEvt.clear();
		}

		vnEvt.push_back(nEvt);

		Mass_ = Mass;
		Uu_ = Uu;
	}	

	//last entry
	mnEvt.push_back(vnEvt);
	vnEvt.clear();

	for (int i = 0; i < vMass.size(); ++i)
	{
		double Thr = Threshold;
		if (MassDep != 0)
		{
			if (vMass.at(i) < MassA)
				Thr += MassDep*MassA;
			else if (vMass.at(i) > MassB)
				Thr += MassDep*MassB;
			else Thr += MassDep*vMass.at(i);
		}

		for (int j = 0; j < vUu.size(); ++j)
		{
			if (mnEvt.at(i).at(j) > Thr)
			{
				Out << vMass.at(i) << "\t" << vUu.at(j) << std::endl;
				break;
			}
		}
	}
	for (int i = vMass.size()-1; i > 0; --i)
	{
		double Thr = Threshold;
		if (vMass.at(i) < MassA)
			Thr = Threshold + MassDep*MassA;
		else if (vMass.at(i) > MassB)
			Thr = Threshold + MassDep*MassB;
		else Thr = Threshold + MassDep*vMass.at(i);

		for (int j = vUu.size()-1; j > 0; --j)
		{
			if (mnEvt.at(i).at(j) > Thr)
			{
				Out << vMass.at(i) << "\t" << vUu.at(j) << std::endl;
				break;
			}
		}
	}

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
