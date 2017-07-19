#include <iostream>
#include <fstream>
#include <vector>
#include <cstring>
#include <getopt.h>

#include "TTree.h"
#include "TFile.h"

#include "Background.h"
#include "Particle.h"
#include "Detector.h"

void Usage(char *Name);
int main(int argc, char** argv)
{

	const struct option longopts[] = 
	{
		{"detconfig", 	required_argument,	0, 'd'},
		{"geniedb", 	required_argument,	0, 'i'},
		{"root", 	required_argument,	0, 'r'},
		{"output", 	required_argument,	0, 'o'},
		{"channel", 	required_argument,	0, 'c'},
		{"help", 	no_argument,	 	0, 'h'},
		{0,	0, 	0,	0},
	};

	int index; 
	int iarg = 0;
	opterr = 1;
	
	std::ofstream OutFile;
	std::string InFile, DetConfig, Channel;
	//std::string RootFile;
	TFile *RootFile;
	
	while((iarg = getopt_long(argc,argv, "d:i:r:o:c:h", longopts, &index)) != -1)
	{
		switch(iarg)
		{
			case 'd':
				DetConfig.assign(optarg);
				break;
			case 'i':
				InFile.assign(optarg);
				break;
			case 'r':
				RootFile = new TFile(optarg, "RECREATE");
				//RootFile.assign(optarg);
				break;
			case 'o':
				OutFile.open(optarg);
				break;
			case 'c':
				Channel.assign(optarg);
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

	TTree *Event = new TTree("Event", "Event");
	//Background *Bkg = new Background(InFile, DetConfig, RootFile, Channel);
	Background *Bkg = new Background(InFile, DetConfig, Channel);

	RootFile->cd();
	Bkg->InitTree();
	int NumWrite = 100;
	for (unsigned int i = 0; i < NumWrite; ++i)	//save at least 100 times
	{
		Bkg->Loop(NumWrite);	//loop over the events...
		//std::cout << "now writing" << std::endl;
		Event = Bkg->GetTree();
		Event->Write();
	}

	RootFile->Close();

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
