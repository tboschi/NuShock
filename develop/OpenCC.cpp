#include <iostream>
#include <fstream>
#include <cstring>
#include <vector>
#include <getopt.h>

#include "Hadron.h"

void Usage(char* argv0);

int main(int argc, char** argv)
{
	const struct option longopts[] = 
	{
		{"", 	required_argument, 	0, 'b'},
		{"confidence", 	required_argument, 	0, 'C'},
		{"help", 	no_argument,	 	0, 'h'},
		{0,	0, 	0,	0},
	};

	int index; 
	int iarg = 0;
	opterr = 1;

	std::string sProb, sTarg("H");
	std::ofstream OutFile;
	double SE = 1000;

	while((iarg = getopt_long(argc,argv, "s:t:o:h", longopts, &index)) != -1)
	{
		switch(iarg)
		{
			case 's':
				SE = strtod(optarg, NULL);
				break;
			case 't':
				sTarg.assign(optarg);
				break;
			case 'o':
				OutFile.open(optarg);
				break;
			case 'h':
				Usage(argv[0]);
				return 1;
			default:
				break;
		}
	}

	std::ostream &Out = (OutFile.is_open()) ? OutFile : std::cout;

	double mp = Const::fMProton;
	double mn = Const::fMNeutron;
	double Mass;
	unsigned int A;
	std::string tName;

	if (sTarg == "H")
	{
		A = 1;
		Mass = mp/A;
		tName.assign("nCTEQ15_1");
	}
	else if (sTarg == "C")
	{
		A = 12;
		Mass = (6*mp + 6*mn)/A;	//12C
		tName.assign("nCTEQ15FullNuc_12_6");
	}
	else
		std::cerr << "Sorry, don't know this mate" << std::endl;

	Hadron *CCxsec = new Hadron("nCTEQ15_1", tName, 4);	//4 is the charm
	CCxsec->SetTarg(0, 0, 0, Mass);	//12C

	for (double cmE = 40; cmE < 450; cmE += 40)
	{
		CCxsec->SetProb(sqrt(cmE*cmE - mp*mp), 0, 0, cmE);
		double Err, Chi;
		Out << cmE << "\t" << sqrt(CCxsec->s()) << "\t" << Const::fGeV2ub * A * CCxsec->Total(Err, Chi) << "\t" << Err << std::endl;
	}

	/*
	//for (double q2 = Min; q2 < Max; q2 += dM/100.0)
	CCxsec->SetParton(SE);
	for (unsigned int i = 0; i < 1000; ++i)
	{
		CCxsec->SetFromPS();
		Out << CCxsec->CosT_() << "\t";
		Out << CCxsec->dXSdOmega_qq() << "\t";
		Out << CCxsec->dXSdOmega_gg() << "\t";
		double pt2 = CCxsec->Pt(2);
		Out << CCxsec->aS(pt2) * log(CCxsec->Q2_() / pt2);
		Out << std::endl;
	}
	*/
	
	if (OutFile.is_open())
		OutFile.close();

	return 0;
}

void Usage(char* argv0)
{
	std::cout << "Description" << std::endl;
	std::cout << "Usage : " << std::endl;
	std::cout << argv0 << " [OPTIONS]" << std::endl;
	std::cout <<"\n  -t,  --target" << std::endl;
	std::cout << "\t\tThe element of the target (available 'H', 'C')" << std::endl;
	std::cout <<"\n  -o,  --output" << std::endl;
	std::cout << "\t\tOutput file" << std::endl;
	std::cout <<"\n  -h,  --help" << std::endl;
	std::cout << "\t\tPrint this message and exit" << std::endl;
}
