#include <iostream>
#include <fstream>
#include <vector>
#include <cstring>
#include <iomanip>
#include <getopt.h>

#include "flux.h"
#include "detector.h"

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
	std::string cutConfig, backPath, scaleConfig, outName;
	double CL;
	
	while((iarg = getopt_long(argc,argv, "s:b:o:c:h", longopts, &index)) != -1)
	{
		switch(iarg)
		{
			case 'c':
				cutConfig.assign(optarg);
				break;
			case 'b':
				backPath.assign(optarg);
				break;
			case 's':
				scaleConfig.assign(optarg);
				break;
			case 'o':
				outName.assign(optarg);
				break;
			case 'C':
				CL = std::strtod(optarg, NULL);
				break;
			case 'h':
				Usage(argv[0]);
				return 1;
			default:
				break;
		}
	}

	CardDealer bc(cutConfig);
	CardDealer sc(scaleConfig);

	std::map<std::string, double> scaleFacts;
	sc.GetAll(scaleFacts);
	std::cout << "scale " << scaleFacts.size() << std::endl;
	

	std::string head = backPath.substr(backPath.find_last_of('/')+1);
	std::cout << "Background head " << head << std::endl;
	std::string channel[6] = {"EPI", "MPI", "nPI0", "nEE", "nEM", "nMM"};
	for (int ch = 0; ch < 6; ++ch)
	{
		std::cout << "Channel " << channel[ch] << std::endl;

		std::map<std::string, std::string> pathMC;
		std::map<std::string, std::string>::iterator ip;
		if (!bc.Get(channel[ch], pathMC))
			continue;

		//loop on background files
		for (ip = pathMC.begin(); ip != pathMC.end(); ++ip)
		{
			std::cout << "Loading " << ip->second << "\t";
			//load MC in efficiency object
			Efficiency eff(ip->second);

			std::vector<std::string> cuts = eff.AvailableCuts(channel[ch]);

			std::string name = channel[ch] + ip->first.substr(0, ip->first.find_first_of('_')+4);
			std::string mass = ip->first.substr(ip->first.find_last_of('_') + 1);
			double m = std::strtod(mass.c_str(), NULL) / 1e3;
			std::cout << "for mass " << m << "\t";

			/*get background files
			 */
			std::string cmd = "ls " + backPath + "*" + name + "*.root > .tmp_backfiles";
			std::cout << "cmd " << cmd << std::endl;
			system(cmd.c_str());
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
			///

			/* create cut file
			 */
			std::string nout = outName + channel[ch] + ip->first;
			std::cout << "and saving in " << nout << std::endl;
			double alpha = 0, cLo, cUp;	// <-- CL
			int totalBackgrounds = 5e6;	//can't be bigger than this
			//try different alphas to improve overall efficiencyes and leaving background
			double step = 0.05;
			for (double CL = 0; CL <= 1; CL += step)
			{
				eff.Reset();
				std::ofstream out(nout.c_str());
				for (int i = 0; i < cuts.size(); ++i)
				{
					if (channel[ch] == "nPI0" && cuts[i] == "E_0")
						eff.FindCut(cuts[i], cLo, cUp, -4.5*(1-CL)+CL, m);
					//else if (channel[ch] == "nPI0" && cuts[i] == "E_A")
						//eff.FindCut(cuts[i], cLo, cUp, 2.0/3.0*(1-CL)+CL, m);
						//eff.FindCut(cuts[i], cLo, cUp, CL, m);
					else if ((channel[ch] == "nEE" ||
						  channel[ch] == "nEM" ||
						  channel[ch] == "nMM") && cuts[i] == "M_0")
						eff.FindCut(cuts[i], cLo, cUp, 0.95*(1-CL)+CL, m/2.0);
					else
						eff.FindCut(cuts[i], cLo, cUp, CL, m);

					out << cuts[i] << "\t" << cLo << "\t" << cUp << std::endl;

					if ((channel[ch] == "nPI0" ||	//particle A & B are identical
					     channel[ch] == "nEE"  ||
					     channel[ch] == "nMM") && cuts[i].find("A") == cuts[i].size()-1)
					{
						std::string var = cuts[i];
						var.replace(var.size()-1, 1, "B");
						eff.SetCut(var, cLo, cUp);

						out << var << "\t" << cLo << "\t" << cUp << std::endl;
					}
				}
				out.close();

				double totBkg = 0;
				double oneBkg = 0;
				for (ib = backFiles.begin(); ib != backFiles.end(); ++ib)
				{
					//std::cout << ib->second << std::endl;
					//std::string nd = ib->first.substr(ib->first.find_first_of('_')+1);
					std::string nd = ib->first.substr(ib->first.find(channel[ch])+channel[ch].length()+1);
					//std::cout << "applying on background " << ib->second << "\n";
					//load MC in efficiency object
					Efficiency bkg(ib->second);
					bkg.LoadCut(nout);
					bkg.ApplyCut();
					std::cout << "\t" << nd << " : " << bkg.EntriesLeft() << ",";
					totBkg += scaleFacts[nd] * bkg.EntriesLeft();
					oneBkg += bkg.EntriesLeft();
				}


				std::cout << std::endl;

				std::cout << "at " << CL << " calculated : " << totBkg << "(" << oneBkg << ")" << std::endl;
				if (totBkg > totalBackgrounds)
				{
					alpha = CL - step;
					break;
				}
				else
					totalBackgrounds = totBkg;
			}
			std::cout << "stopped at " << alpha << " and " << totalBackgrounds << std::endl;

			eff.Reset();
			std::ofstream out(nout.c_str());
			for (int i = 0; i < cuts.size(); ++i)
			{
				if (channel[ch] == "nPI0" && cuts[i] == "E_0")
						eff.FindCut(cuts[i], cLo, cUp, -4.5*(1-alpha)+alpha, m);
				//else if (channel[ch] == "nPI0" && cuts[i] == "E_A")
						//eff.FindCut(cuts[i], cLo, cUp, 2.0/3.0*(1-alpha)+alpha, m);
						//eff.FindCut(cuts[i], cLo, cUp, alpha, m);
				else if ((channel[ch] == "nEE" ||
					  channel[ch] == "nEM" ||
					  channel[ch] == "nMM") && cuts[i] == "M_0")
					eff.FindCut(cuts[i], cLo, cUp, 0.95*(1-alpha)+alpha, m/2.0);
				else
					eff.FindCut(cuts[i], cLo, cUp, alpha, m);

				out << cuts[i] << "\t" << cLo << "\t" << cUp << std::endl;

				if ((channel[ch] == "nPI0" ||	//particle A & B are identical
				     channel[ch] == "nEE"  ||
				     channel[ch] == "nMM") && cuts[i].find("A") == cuts[i].size()-1)
				{
					std::string var = cuts[i];
					var.replace(var.size()-1, 1, "B");
					eff.SetCut(var, cLo, cUp);

					out << var << "\t" << cLo << "\t" << cUp << std::endl;
				}
			}

			eff.ApplyCut();
			out << "\n#efficiency\t" << eff.ReductionFactor() << std::endl;
			out << "\n#entries\t\t" << eff.EntriesLeft() << std::endl;
			std::cout << "efficiency " << eff.EntriesLeft() * 1e-4 << std::endl;
			out.close();
		}
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
