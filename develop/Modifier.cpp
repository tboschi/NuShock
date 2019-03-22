#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <algorithm>
#include <getopt.h>

#include "TFile.h"
#include "TH1D.h"

class Sorter
{
	private:
		std::vector<double> vInd;
	public:
		Sorter(std::vector<double>& vvect) : vInd(vvect) {}
		bool operator()(int i, int j) const { return vInd.at(i) < vInd.at(j); }
};

int main(int argc, char** argv)
{
	const struct option longopts[] = 
	{
		{"listfile", 	required_argument, 	0, 'l'},
		{"output", 	required_argument, 	0, 'o'},
		{"help", 	no_argument,	 	0, 'h'},
		{0,	0, 	0,	0},
	};

	int index; 
	int iarg = 0;
	opterr = 1;

	//std::string sProb, sTarg("H"), OutName;
	std::ifstream ListFile;
	std::string OutBase;

	while((iarg = getopt_long(argc,argv, "l:o:h", longopts, &index)) != -1)
	{
		switch(iarg)
		{
			case 'l':
				ListFile.open(optarg);
				break;
			case 'o':
				OutBase.assign(optarg);
				break;
			case 'h':
				return 1;
			default:
				break;
		}
	}

	std::stringstream ssL;
	double Mass;
	std::string Line, File;
	
	std::vector<double> vM;
	std::vector<std::string> vF;
	std::vector<unsigned int> vI;
	while (std::getline(ListFile, Line))
	{
		if (Line[0] == '#') continue;

		ssL.str("");
		ssL.clear();
		ssL << Line;

		if (ssL >> File)
		{
			//vM.push_back(Mass);
			vF.push_back(File);
			std::size_t Pos_ = File.find("tau_");
			std::size_t Pos0 = File.find("_0.root");
			std::size_t PosB = File.find("_B.root");
			std::size_t End_;
			if (Pos0 != std::string::npos)
				End_ = Pos0;
			else if (PosB != std::string::npos)
				End_ = PosB;

			std::string mass = File.substr(Pos_+4, End_-Pos_-4);
			Mass = std::strtod(mass.c_str(), NULL) / 1000.0;
			vM.push_back(Mass);
			vI.push_back(vI.size());
		}
	}

	std::sort(vI.begin(), vI.end(), Sorter(vM));

	std::ofstream Charm(std::string(OutBase + "_Charm.dat").c_str());
	std::ofstream TauE (std::string(OutBase + "_TauE.dat").c_str());
	std::ofstream TauM (std::string(OutBase + "_TauM.dat").c_str());
	std::ofstream Pion (std::string(OutBase + "_Pion.dat").c_str());
	std::ofstream PPion(std::string(OutBase + "_2Pion.dat").c_str());

	std::ofstream *Out;

	double Start = 0, End = 0, Peak = 0, Area = 0;
	std::vector<double> vPeak(5), vArea(5);
	std::string Name;
	TH1D *hFlux = 0;
	TObject *X = 0;
	TFile *InFile;

	for (unsigned int j = 0; j < vM.size(); ++j)
	{
		unsigned int i = vI.at(j);
		InFile = new TFile(vF.at(i).c_str(), "OPEN");

		for (unsigned int s = 0; s < 5; ++s)
		{
			switch (s)
			{
				case 0:
					Name = "hcharm";
					Out = &Charm;
					break;
				case 1:
					Name = "htaue";
					Out = &TauE;
					break;
				case 2:
					Name = "htaum";
					Out = &TauM;
					break;
				case 3:
					Name = "hpion";
					Out = &Pion;
					break;
				case 4:
					Name = "h2pion";
					Out = &PPion;
					break;
			}

			Start = 0;
			End   = 0;
			Peak  = 0;
			Area  = 0;

			if (X = InFile->Get(Name.c_str()))
			{
				hFlux = dynamic_cast<TH1D*>(X);
				if (hFlux->GetEntries())
				{
					Start = hFlux->GetXaxis()->GetBinCenter(hFlux->FindFirstBinAbove(0));
					End   = hFlux->GetXaxis()->GetBinCenter(hFlux->FindLastBinAbove(0));
					Peak  = hFlux->GetMaximum();
					Area  = hFlux->Integral("WIDTH");
				}
				else
				{
					Start = sqrt(-1);
					End = sqrt(-1);
					Peak = sqrt(-1);
					Area = sqrt(-1);
				}

				if (vM.at(i) == 0.0)
				{
					vPeak.at(s) = Peak;
					vArea.at(s) = Area;
				}

				*Out << vM.at(i) << "\t";
				*Out << Start << "\t";
				*Out << End << "\t";
				*Out << Peak/vPeak.at(s) << "\t";
				*Out << Area/vArea.at(s) << std::endl;
			}
		}
			
		InFile->Close();
	}

	return 0;
}
