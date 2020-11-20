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
	bool chargeID = false, verbose = false;
	
	std::string inFile, detConfig, module, channel, outName;
	
	while((iarg = getopt_long(argc,argv, "d:i:o:l:c:Cvh", longopts, &index)) != -1)
	{
		switch(iarg)
		{
			case 'd':
				detConfig.assign(optarg);
				break;
			case 'l':
				module.assign(optarg);
				break;
			case 'i':
				inFile.assign(optarg);
				break;
			case 'o':
				outName.assign(optarg);
				//RootFile.assign(optarg);
				break;
			case 'c':
				channel.assign(optarg);
				break;
			case 'C':
				chargeID = true;
				break;
			case 'v':
				verbose = true;
				break;
			case 'h':
				Usage(argv[0]);
				return 1;
			default:
				break;
		}
	
	}

	std::string catChannel = "_" + channel;
	std::string catModule  = "_" + module;
	outName.insert(outName.find(".root"), catChannel);
	outName.insert(outName.find(".root"), catModule);

	if (inFile.find("numu") != std::string::npos)
		outName.insert(outName.find(".root"), "_M");
	else if (inFile.find("nue") != std::string::npos)
		outName.insert(outName.find(".root"), "_E");
	else if (inFile.find("nubar") != std::string::npos)
		outName.insert(outName.find(".root"), "_B");

	GenieBack *bkg = new GenieBack(inFile, outName, verbose);
	Tracker *detector = new Tracker(detConfig, module);

	// here you need to implement your definition of process, i.e. what final state you are looking for
	// I define it as a map of 2 integers: a pdg code and the number of particles I want
	// e.g. final state 2 muons, 1 elec

	//sum of charges of the two particles is 0
	std::map<int, int> process;
	if (channel == "nEE")
		process[11] = 2;
	else if (channel == "nME" || channel == "nEM")
	{
		process[11] = 1;
		process[13] = 1;
	}
	else if (channel == "nMM")
		process[13] = 2;
	else if (channel == "nPI0")
		process[22] = 2;
	else if (channel == "EPI")
	{
		process[11] = 1;
		process[211] = 1;
	}
	else if (channel == "MPI")
	{
		process[13] = 1;
		process[211] = 1;
	}

	TTree *genie;

	int save = 10;
	genie = bkg->FindBackground(detector, process, save);	//loop over the events..

	delete bkg;
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
