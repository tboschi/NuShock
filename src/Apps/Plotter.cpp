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
	double RedChi2 = 117.407/99.0;		//90 C.L. for 99 d.o.f.
	bool Efficiency = false;
	
	while((iarg = getopt_long(argc,argv, "i:o:t:Wh", longopts, &index)) != -1)
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
			case 'W':
				Efficiency = true;
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
	for (double logMass = -2.0; logMass < -0.3; logMass += 0.0034)	//increase mass log
		vMass.push_back(pow(10.0, logMass));
	for (double logUu2 = -10.0; logUu2 < -0.0; logUu2 += 0.02)	//increase Uu logarithmically
		vUu.push_back(pow(10.0, logUu2));

	std::vector<std::vector<double> > mnEvt, mSignal, mBackground, mChi2;;
	std::vector<double> vnEvt, vSignal, vBackground, vChi2;

	double Mass_ = -1, Uu_ = -1, nEvt_ = -1;
	double Mass = 0, Uu = 0, nEvt = 0, Sign = 0, Back = 0, Chi2 = 0;
	std::string Line, Key, Name;
	std::stringstream ssL;
	while (std::getline(InFile, Line))
	{
		if (Line[0] == '#') continue;

		ssL.str("");
		ssL.clear();
		ssL << Line;

		ssL >> Mass >> Uu >> nEvt >> Sign >> Back >> Chi2;

		if (Mass != Mass_ && Mass_ >= 0)
		{
			mnEvt.push_back(vnEvt);
			mSignal.push_back(vSignal);
			mBackground.push_back(vBackground);
			mChi2.push_back(vChi2);

			vnEvt.clear();
			vSignal.clear();
			vBackground.clear();
			vChi2.clear();
		}

		vnEvt.push_back(nEvt);
		vSignal.push_back(Sign);
		vBackground.push_back(Back);
		vChi2.push_back(Chi2);

		Mass_ = Mass;
		Uu_ = Uu;
	}	

	//last entry
	mnEvt.push_back(vnEvt);
	mSignal.push_back(vSignal);
	mBackground.push_back(vBackground);
	mChi2.push_back(vChi2);

	vnEvt.clear();
	vSignal.clear();
	vBackground.clear();
	vChi2.clear();

	for (int i = 0; i < vMass.size(); ++i)
	{
		for (int j = 0; j < vUu.size(); ++j)
		{
			if (!Efficiency)
			{
				if (mnEvt.at(i).at(j) > Threshold)
				{
					Out << vMass.at(i) << "\t" << vUu.at(j) << std::endl;
					break;
				}
			}
			else 
			{
				if (mChi2.at(i).at(j) > RedChi2)	//threshold becomes a redchi2
				{
					Out << vMass.at(i) << "\t" << vUu.at(j) << std::endl;
					break;
				}
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
	std::cout <<"\n  -t" << std::endl;
	std::cout << "\t\tevent trheshold for signal" << std::endl;
	std::cout <<"\n  -W" << std::endl;
	std::cout << "\t\treduced chi2 will be used at 90\% C.L." << std::endl;
	std::cout <<"\n  -h,  --help" << std::endl;
	std::cout << "\t\tPrint this message and exit" << std::endl;
}
