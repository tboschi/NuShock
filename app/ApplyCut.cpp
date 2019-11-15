#include <iostream>
#include <fstream>
#include <vector>
#include <cstring>
#include <iomanip>
#include <getopt.h>

#include "flux.h"
#include "detector.h"

#include "TFile.h"
#include "TTree.h"
#include "TDirectory.h"
#include "TKey.h"
#include "TH1D.h"
#include "TGraph.h"

void Usage(char* argv0);
int main(int argc, char** argv)
{
	const struct option longopts[] = 
	{
		{"detconfig", 	required_argument,	0, 'd'},
		{"xsection", 	required_argument, 	0, 'x'},
		{"flux", 	required_argument,	0, 'f'},
		{"output", 	required_argument,	0, 'o'},
		{"help", 	no_argument,	 	0, 'h'},
		{0,	0, 	0,	0},
	};

	int index; 
	int iarg = 0;
	opterr = 1;
	
	//Initialize variables
	std::ofstream out;
	std::string backFile, signFile, cutFile;
	
	while((iarg = getopt_long(argc,argv, "b:c:s:o:h", longopts, &index)) != -1)
	{
		switch(iarg)
		{
			case 'b':
				backFile.assign(optarg);
				break;
			case 's':
				signFile.assign(optarg);
				break;
			case 'c':
				cutFile.assign(optarg);
				break;
			case 'o':
				out.open(optarg);
				break;
			case 'h':
				Usage(argv[0]);
				return 1;
			default:
				break;
		}
	}

	TFile bf(backFile.c_str());
	TTree *data = dynamic_cast<TTree*>(bf.Get("Data"));
	std::cout << "background before cause  " << data->GetEntries() << std::endl;
	bf.Close();

	Efficiency effBack(backFile);
	Efficiency effSign(signFile);
	effBack.LoadCut(cutFile);
	effBack.ApplyCut();
	effSign.LoadCut(cutFile);
	effSign.ApplyCut();
	std::cout << "background events left are " << effBack.EntriesLeft()
		  << " (" << effBack.ReductionFactor() << ")\n";
	std::cout << "simulation events left are " << effSign.EntriesLeft()
		  << " (" << effSign.ReductionFactor() << ")\n";

	return 0;
}

void Usage(char* argv0)
{
	std::cout << "Description" << std::endl;
	std::cout << "Usage : " << std::endl;
	std::cout << argv0 << " [OPTIONS]" << std::endl;
	std::cout <<"\n  -d,  --detconfig" << std::endl;
	std::cout << "\t\tDetector configuration file" << std::endl;
	std::cout <<"\n  -x,  --xsection" << std::endl;
	std::cout << "\t\tDetector configuration file" << std::endl;
	std::cout <<"\n  -f,  --flux" << std::endl;
	std::cout << "\t\tFlux configuration file" << std::endl;
	std::cout <<"\n  -o,  --output" << std::endl;
	std::cout << "\t\tOutput file, to save message" << std::endl;
	std::cout <<"\n  -h,  --help" << std::endl;
	std::cout << "\t\tPrint this message and exit" << std::endl;
}
