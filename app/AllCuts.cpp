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
	std::string backPath, detConfig, signFile, cutFile, scaleFile, mass;
	
	while((iarg = getopt_long(argc,argv, "b:c:d:s:e:m:o:h", longopts, &index)) != -1)
	{
		switch(iarg)
		{
			case 'm':
				mass.assign(optarg);
				break;
			case 'b':
				backPath.assign(optarg);
				break;
			case 'd':
				detConfig.assign(optarg);
				break;
			case 's':
				signFile.assign(optarg);
				break;
			case 'e':
				scaleFile.assign(optarg);
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

	Detector *nd = new Detector(detConfig);
	std::map<std::string, double> mw;
	mw["LAr"] = nd->WeightLAr() / nd->Weight();
	mw["FGT"] = nd->WeightFGT() / nd->Weight();

	CardDealer cut(cutFile);
	CardDealer sig(signFile);

	std::map<std::string, double> sfact;
	CardDealer sca(scaleFile);
	sca.GetAll(sfact);

	std::string channel[6] = {"EPI", "MPI", "nPI0", "nEE", "nEM", "nMM"};
	std::string module[2] = {"LAr", "FGT"};
	std::string fermion[2] = {"dirac", "major"};
	std::string type[3] = {"E", "M", "B"};
	for (int c = 0; c < 6; ++c)
		for (int f = 0; f < 2; ++f)
	{
		std::cout << channel[c] << ", " << fermion[f] << std::endl;

		double mx0 = 0, mx1 = 0, totx = 0;
		std::map<std::string, double> mn0, mn1;
		std::map<std::string, double>::iterator im;
		for (int t = 0; t < 3; ++t)
		{
			mn0[type[t]] = 0;
			mn1[type[t]] = 0;
			for (int m = 0; m < 2; ++m)
			{
				std::string skey = module[m] + "_" + type[t];
				std::string backfile = backPath + channel[c] + "_" + skey + ".root";

				std::string cutfile;
				std::string cutkey = channel[c] + "_" + module[m] + "_" +
						     fermion[f] + "_" + mass;
				cut.Get(cutkey, cutfile);

				Efficiency effBack(backfile);
				effBack.LoadCut(cutfile);
				effBack.ApplyCut();

				mn0[type[t]] += mw[module[m]] * effBack.GetEntries() / 1000.0;
				mn1[type[t]] += mw[module[m]] * effBack.EntriesLeft() / 1000.0;

				mx0 += sfact[skey] * mw[module[m]] * effBack.GetEntries() / 1000.0;
				mx1 += sfact[skey] * mw[module[m]] * effBack.EntriesLeft() / 1000.0;
				totx += sfact[skey];
			}

			std::cout << "nu_" << type[t] << " : "
				  << mn0[type[t]] << " --> " << mn1[type[t]] << std::endl;
		}
		std::cout << "______________________________________\n";
		std::cout << "avg  : " << mx0 / totx << " --> " << mx1 / totx << std::endl;

		double red = 0;
		for (int m = 0; m < 2; ++m)
		{
			std::string key = channel[c] + "_" + module[m] + "_" +
					  fermion[f] + "_" + mass;
			std::string sigfile, cutfile;
			cut.Get(key, cutfile);
			sig.Get(key, sigfile);

			Efficiency effSign(sigfile);
			effSign.LoadCut(cutfile);
			effSign.ApplyCut();

			red += mw[module[m]] * effSign.ReductionFactor() * 100.0;
		}

		std::cout << "efficiency : " << red << "%" << std::endl;
		std::cout << std::endl;
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
