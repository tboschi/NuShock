#include <iostream>
#include <fstream>
#include <cstring>
#include <sstream>
#include <vector>
#include <cmath>
#include <getopt.h>
#include <cstdlib>

#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"

#include "Efficiency.h"

void Usage(char* argv0);
int main(int argc, char** argv)
{
	const struct option longopts[] = 
	{
		{"input", 	required_argument,	0, 'i'},
		{"output", 	required_argument,	0, 'o'},
		{"highenergy", 	no_argument,		0, 'E'},
		{"complete", 	no_argument,		0, 'C'},
		{"timing", 	no_argument,		0, 'T'},
		{"help", 	no_argument,	 	0, 'h'},
		{0,	0, 	0,	0},
	};

	int index; 
	int iarg = 0;
	opterr = 1;
	
	//Initialize variables
	std::string InFile;
	TFile *OutFile;
	bool Complete = false, Time = false;
	
	while((iarg = getopt_long(argc,argv, "i:o:CTh", longopts, &index)) != -1)
	{
		switch(iarg)
		{
			case 'i':
				InFile.assign(optarg);	//will contain all links to simulation files, mass, cuts..
				break;
			case 'o':
				OutFile = new TFile(optarg, "RECREATE");
				break;
			case 'C':
				Complete = true;
				break;
			case 'T':
				Time = true;
				break;
			case 'h':
				Usage(argv[0]);
				return 1;
			default:
				return 1;
				break;
		}
	}


	TH2D *TheFunction;
	TH1D *TheAll, *TheCut;
	Efficiency *MakeEff = new Efficiency(InFile, Time);
	
	MakeEff->InitFunc();
	MakeEff->LoopFile();

	if (Complete)
		MakeEff->CompleteFunction();

	OutFile->cd();
	TheFunction = MakeEff->GetFunction();
	TheAll = MakeEff->GetAll();
	TheCut = MakeEff->GetCut();

	TheFunction->Write();
	TheAll->Write();
	TheCut->Write();
	std::cout << "Survive " << TheCut->GetEntries() / double(TheAll->GetEntries()) << std::endl;

	OutFile->Close();

	return 0;
}
	
void Usage(char* argv0)
{
	std::cout << "Description" << std::endl;
	std::cout << "Usage : " << std::endl;
	std::cout << argv0 << " [OPTIONS]" << std::endl;
	std::cout <<"\n  -i,  --input" << std::endl;
	std::cout << "\t\tInput file, contains list of simulation files" << std::endl;
	std::cout <<"\n  -o,  --output" << std::endl;
	std::cout << "\t\tROOT output file, where to save the function" << std::endl;
	std::cout <<"\n  -C,  --complete" << std::endl;
	std::cout << "\t\tExtrapolate function and complete it for any mass value" << std::endl;
	std::cout <<"\n  -E,  --highenergy" << std::endl;
	std::cout << "\t\tConsiders only events higher than 5 GeV" << std::endl;
	std::cout <<"\n  -h,  --help" << std::endl;
	std::cout << "\t\tPrint this message and exit" << std::endl;
}
