#include <iostream>
#include <fstream>
#include <vector>
#include <cstring>
#include <iomanip>
#include <getopt.h>

#include "flux.h"
#include "detector.h"
#include "analysis.h"

#include "TFile.h"
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
	std::string listConfig, simPath, cutPath, outPath;
	bool charge = false, only = false, special = false;
	
	while((iarg = getopt_long(argc,argv, "i:s:c:o:CSOh", longopts, &index)) != -1)
	{
		switch(iarg)
		{
			case 'i':
				listConfig.assign(optarg);
				break;
			case 's':
				simPath.assign(optarg);
				break;
			case 'c':
				cutPath.assign(optarg);
				break;
			case 'o':
				outPath.assign(optarg);
				break;
			case 'C':
				charge = true;
				break;
			case 'S':
				special = true;
				break;
			case 'O':
				only = true;
				break;
			case 'h':
				Usage(argv[0]);
				return 1;
			default:
				break;
		}
	}

	CardDealer lc(listConfig);

	int maxchl = charge ? 2 : 6;
	int maxlnv = 1 + charge;

	std::string channel[6] = {"EPI", "MPI", "nPI0", "nEE", "nEM", "nMM"};
	std::string module[2] = {"LAr", "FGT"};
	std::string fermion[2] = {"dirac", "major"};
	std::string lnv[2] = {"LNV", "LNC"};
	for (int c = 0; c < maxchl; ++c)
		for (int m = 0; m < 2; ++m)
			for (int f = 0; f < 2; ++f)
				for (int l = 0; l < maxlnv; ++l)
	{
		int LN = charge ? 2*l - 1 : 0; 	//is ±1

		std::cout << "Channel " << channel[c] << " (" << lnv[l] << "), "
			  << module[m] << ", " << fermion[f]  << std::endl;

		//std::string key = channel[c] + "_" +
		//		  module[m] + "_" + fermion[f];
		std::string key = channel[c] + "_" +
				  (charge ? lnv[l] + "_" : "") +
				  module[m] + "_" + fermion[f] + "_";
		std::string sim = channel[c] + "_" + 
				  module[m] + "_" + fermion[f] + "_";
		std::map<std::string, std::string> names;
		std::map<std::string, std::string>::iterator in;

		std::cout << "continue" << std::endl;
		if (!lc.Get(key, names))
			continue;
		std::cout << "no " << names.size() << std::endl;

		std::string outName = outPath;

		if (outName.find_last_of('.') != std::string::npos)
			outName.insert(outName.find_last_of('.'), key, 0, key.size()-1);
		else
			outName.insert(outName.size(), key, 0, key.size()-1);

		TFile *outFile = new TFile(outName.c_str(), "RECREATE");
		Efficiency *eff = new Efficiency("", LN, special);

		eff->MakeFunction();

		for (in = names.begin(); in != names.end(); ++in)
		{
			std::string fullName = sim + in->first;
			std::cout << "Loading " << fullName << std::endl;

			std::string simName = simPath;
			std::string cutName = cutPath;

			if (simName.find_last_of('.') != std::string::npos)
				simName.insert(simName.find_last_of('.'), fullName);
			else
				simName.append(fullName + ".root");

			//if (cutName.find_last_of('.') != std::string::npos)
			//	cutName.insert(cutName.find_last_of('.'), fullName);
			//else
			//	cutName.append(fullName);


			std::cout << "simulation " << simName << std::endl;
			std::cout << "cut file   " << in->second << std::endl;
			eff->LoadFile(simName);

			//std::string mass = in->first.substr(in->first.find_last_of('_')+1);
			double mass = std::strtod(in->first.c_str(), NULL) / 1e3;
			std::cout << "with mass " << mass << std::endl;

			if (eff->ValidEntries() > 0)
			{
				if (!only)
					eff->LoadCut(in->second);
				eff->ApplyCut(mass);
			}
			else	//if no valid entries set minimum efficiency at this mass
				eff->MinEfficiency(mass);
		}

		std::cout << "Completing function" << std::endl;
		TH2D *hhf  = eff->CompleteFunction();

		outFile->cd();
		hhf->Write();

		outFile->Close();
		delete eff;
	}

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
