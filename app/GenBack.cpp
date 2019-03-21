#include <iostream>
#include <fstream>
#include <vector>
#include <cstring>
#include <getopt.h>

#include "TTree.h"
#include "TFile.h"

#include "tools.h"
#include "detector.h"
#include "background.h"

void Usage(char *Name);
int main(int argc, char** argv)
{

	const struct option longopts[] = 
	{
		{"detconfig", 	required_argument,	0, 'd'},
		{"geniefile", 	required_argument,	0, 'i'},
		{"output", 	required_argument,	0, 'r'},
		{"channel", 	required_argument,	0, 'c'},
		{"help", 	no_argument,	 	0, 'h'},
		{0,	0, 	0,	0},
	};

	int index; 
	int iarg = 0;
	opterr = 1;
	
	std::string inFile, detConfig, channel;
	//std::string RootFile;
	TFile *RootFile;
	
	while((iarg = getopt_long(argc,argv, "d:i:r:c:h", longopts, &index)) != -1)
	{
		switch(iarg)
		{
			case 'd':
				detConfig.assign(optarg);
				break;
			case 'i':
				inFile.assign(optarg);
				break;
			case 'r':
				RootFile = new TFile(optarg, "RECREATE");
				//RootFile.assign(optarg);
				break;
			case 'c':
				channel.assign(optarg);
				break;
			case 'h':
				Usage(argv[0]);
				return 1;
			default:
				break;
		}
	
	}


	GenieBack *bkg = new GenieBack(inFile);
	Tracker *detector = new Tracker(detConfig);

	// here you need to implement your definition of process, i.e. what final state you are looking for
	// I define it as a map of 2 integers: a pdg code and the number of particles I want
	// e.g. final state 2 muons, 1 elec
	std::map<int, int> process;
	if (channel == "MME")
	{
		process[13] = 2;		//muon: PDG 13, number 2
		process[11] = 1;		//elec: PDG 11, number 1
	}

	TTree *genie;
	RootFile->cd();

	int save = 10;
	for (unsigned int i = 0; i < save; ++i)
	{
		genie = bkg->FindBackground(detector, process, save);	//loop over the events..
		genie->Write();
	}

	RootFile->Close();

	return 0;
}
	
void Usage(char *Name)
{
	std::cout << "Description" << std::endl;
	std::cout << "Usage : " << std::endl;
	std::cout << Name << " [OPTIONS]" << std::endl;
	std::cout <<"\n  -d,  --detconfig" << std::endl;
	std::cout << "\t\tDetector config file" << std::endl;
	std::cout <<"\n  -i,  --geniefile" << std::endl;
	std::cout << "\t\tGENIE input file, converted with gntpc" << std::endl;
	std::cout <<"\n  -o,  --output" << std::endl;
	std::cout << "\t\tLog output file" << std::endl;
	std::cout <<"\n  -r,  --root" << std::endl;
	std::cout << "\t\tROOT output file with background tree" << std::endl;
	std::cout <<"\n  -c,  --channel" << std::endl;
	std::cout << "\t\tChannel selection" << std::endl;
	std::cout <<"\n  -h,  --help" << std::endl;
	std::cout << "\t\tPrint this message and exit" << std::endl;
}
