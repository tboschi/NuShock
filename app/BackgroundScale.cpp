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
	std::string backConfig, backPath, scaleConfig;
	double CL = 0.9;
	bool charge = false, special = false;
	
	while((iarg = getopt_long(argc,argv, "b:s:r:o:l:CSh", longopts, &index)) != -1)
	{
		switch(iarg)
		{
			case 'b':
				backConfig.assign(optarg);
				break;
			case 's':
				scaleConfig.assign(optarg);
				break;
			case 'r':
				backPath.assign(optarg);
				break;
			case 'o':
				out.open(optarg);
				break;
			case 'l':
				CL = std::strtod(optarg, NULL);
				break;
			case 'C':
				charge = true;
				break;
			case 'S':
				special = true;
				break;
			case 'h':
				Usage(argv[0]);
				return 1;
			default:
				break;
		}
	}

	int maxchl = charge ? 2 : 6;
	int maxlnv = 1 + charge;
	std::cout << "max channel " << maxchl << std::endl;

	CardDealer bc(backConfig);
	CardDealer sc(scaleConfig);

	std::map<std::string, double> scaleFacts;
	sc.GetAll(scaleFacts);

	std::string channel[6] = {"EPI", "MPI", "nPI0", "nEE", "nEM", "nMM"};
	std::string module[2] = {"LAr", "FGT"};
	std::string fermion[2] = {"dirac", "major"};
	std::string lnv[2] = {"LNV", "LNC"};
	for (int c = 0; c < maxchl; ++c)
		for (int m = 0; m < 2; ++m)
	{
		std::string cmd = "ls " + backPath + "*" + channel[c]
				+ "_" + module[m] + "*.root > .tmp_backfiles";
		system(cmd.c_str());

		std::string head = backPath.substr(backPath.find_last_of('/')+1);
		std::map<std::string, std::string> backFiles;
		std::map<std::string, std::string>::iterator ib;
		std::string line;
		std::ifstream inback(".tmp_backfiles");
		while (std::getline(inback, line))
		{
			int p0 = line.find(head) + head.length();
			int p1 = line.find(".root");
			backFiles[line.substr(p0, p1-p0)] = line;
		}

		for (ib = backFiles.begin(); ib != backFiles.end(); ++ib)
		{
			std::cout << "Loading " << ib->first << "\t"
				  << ib->second << std::endl;
			std::string nd = ib->first.substr(ib->first.find_first_of('_')+1);
			double scale = scaleFacts[nd];

			for (int f = 0; f < 2; ++f)
				for (int l = 0; l < maxlnv; ++l)
			{
				int LN = charge ? 2*l - 1 : 0; 	//is ±1

				//load MC in efficiency object
				Efficiency eff(ib->second, LN, special);

				std::cout << "Channel " << channel[c] << ", "
					<< module[m] << ", " << fermion[f]  << std::endl;

				std::string key = channel[c] + "_" + (charge ? lnv[l] + "_" : "") +
						  module[m] + "_" + fermion[f];
				std::map<std::string, std::string> paths;
				std::map<std::string, std::string>::iterator ip;

				if (!bc.Get(key, paths))
					continue;

				//apply cuts to count background
				std::vector<double> masses, events, signal;
				for (ip = paths.begin(); ip != paths.end(); ++ip)
				{
					eff.Reset();
					std::string mass = ip->first.substr(ip->first.find_last_of('_')+1);

					double evt = 0;
					std::cout << "Applying " << mass << "\t"
						<< ip->second << std::endl;
					if (eff.ValidEntries() > 0)
					{
						std::cout << "events left " << eff.EventsLeft() << std::endl;
						eff.LoadCut(ip->second);
						eff.ApplyCut();
						std::cout << "events left " << eff.EventsLeft() << std::endl;
						evt = scale * eff.EventsLeft();
					}
					//double s = Belt(evt, CL);

					//events.push_back(evt);
					//signal.push_back(s);
	
					std::cout << "events " << evt << std::endl;
					if (evt == 0)
						evt = scale/sqrt(12.0);

					std::cout << "events " << evt << std::endl;
					signal.push_back(evt);
					masses.push_back(std::strtod(mass.c_str(), NULL) / 1e3);
				}

				//offset on y direction to improve fit
				double offset = signal.front();
				for (int e = 0; e < signal.size(); ++e)
					signal[e] -= offset;

				//fit with polynominal line and add offset back again
				int ord = 2;
				PolyFit pf(masses, signal, ord);
				masses = pf.GetAxis();
				signal = pf.GetData();
				std::vector<double> beta = pf.LeastSquare(0);
				beta[0] += offset;

				//out << "# ";
				out << key + ib->first.substr(ib->first.find_last_of('_'));
				for (int b = 0; b < ord+1; ++b)
					out << "\t" << beta[b];
				out << std::endl;
				//for (int m = 0; m < masses.size(); ++m)
				//	out << masses[m] << "\t" << signal[m] + offset << "\t"
				//		<< pf.ff(masses[m]) + offset << std::endl;
				//out << std::endl << std::endl;
			}
		}
	}

	out.close();

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
