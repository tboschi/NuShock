#include <iostream>
#include <fstream>
#include <cstring>
#include <vector>
#include <getopt.h>

//#include "Hadron.h"
#include "Tools.h"
#include "Detector.h"
#include "Physics.h"
#include "Flux.h"

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TRandom3.h"
#include "TGenPhaseSpace.h"

void Usage(char* argv0);
double lDparam(double xf, double pt)
{
	//from 1708.08700, 250GeV proton beam E796
	double b = 1.08;
	double n = 6.1;
	return n*(1-std::abs(xf)) - b * pt*pt;
}

int main(int argc, char** argv)
{
	const struct option longopts[] = 
	{
		{"", 	required_argument, 	0, 'b'},
		{"confidence", 	required_argument, 	0, 'C'},
		{"help", 	no_argument,	 	0, 'h'},
		{0,	0, 	0,	0},
	};

	int index; 
	int iarg = 0;
	opterr = 1;

	//std::string sProb, sTarg("H"), OutName;
	std::string OutName, DetConfig;
	std::ofstream OutFile;
	double BeamE = 800;	//DONUT
	unsigned int nMAX = 1e5;
	double NuMass = 0.0;	//neutrino mass

	while((iarg = getopt_long(argc,argv, "r:s:m:E:t:o:I:d:h", longopts, &index)) != -1)
	{
		switch(iarg)
		{
			case 'r':
				OutName.assign(optarg);
				break;
			case 'm':
				NuMass = strtod(optarg, NULL);
				break;
			case 'E':
				BeamE = strtod(optarg, NULL);
				break;
			case 'o':
				OutFile.open(optarg);
				break;
			case 'I':
				nMAX = strtol(optarg, NULL, 10);
				break;
			case 'd':
				DetConfig.assign(optarg);
				break;
			case 'h':
				Usage(argv[0]);
				return 1;
			default:
				break;
		}
	}

	std::ostream &Out = (OutFile.is_open()) ? OutFile : std::cout;

	std::string Out_ = OutName + ".root";
	std::string OutB = OutName + "Bar.root";

	TFile *FileOut_ = new TFile(Out_.c_str(), "RECREATE");
	TFile *FileOutB = new TFile(OutB.c_str(), "RECREATE");

	TH1D * hTotal_ = new TH1D("htotal", "total",  100, 0, 20);
	TH1D * hTotalB = new TH1D("htotal", "total",  100, 0, 20);

	//neutrino
	TH1D * hCharm = new TH1D("hcharm", "charm",  100, 0, 20);

	//antineutrino
	TH1D * hTauE = new TH1D("htaue",  "taue",   100, 0, 20);
	TH1D * hTauM  = new TH1D("htaum",  "taum",   100, 0, 20);
	TH1D * hPion  = new TH1D("hpion",  "pion",   100, 0, 20);
	TH1D * h2Pion = new TH1D("h2pion", "2 pion", 100, 0, 20);

	Detector *TheBox = new Detector(DetConfig);
	double Th0 = atan2(2*TheBox->Get("Height"), TheBox->Get("Baseline"));	//angular acceptance of detector

	TRandom3 *Gen = new TRandom3(0);

	TLorentzVector Beam(0, 0, sqrt(pow(BeamE, 2) - pow(Const::fMProton, 2)), BeamE);
	TLorentzVector Targ(0, 0, 0, Const::fMProton);
	TLorentzVector S = Beam+Targ;
	double sqrts = S.M();		//CM energy

	Neutrino *Nu_ = new Neutrino(NuMass, Neutrino::Dirac | Neutrino::Unpolarised);
	Neutrino *NuB = new Neutrino(NuMass, Neutrino::Dirac | Neutrino::Unpolarised | Neutrino::Antiparticle);
	
	//find max of Ds param
	double ptmax = sqrts;
	double maxF = lDparam(-1, 0);
	double minF = lDparam(1, ptmax);

	std::vector<Particle*> vProducts;
	unsigned int nIter = 0;
	while (nIter < nMAX)
	{
		std::cout << "Simulating " << nIter << std::endl;
		vProducts.clear();
		double pt = Gen->Uniform(0, ptmax);
		double xf = Gen->Uniform(-1.0, 1.0);
		std::cout << minF << "\t" << maxF << "\t" << lDparam(xf, pt) << std::endl;
		std::cout << "H0" << std::endl;
		if (Gen->Uniform(0, pow(10,maxF)) < pow(10, lDparam(xf, pt)))
		{
		std::cout << "H1" << std::endl;
			double px, py;
			double pz = sqrts*xf*0.5;
			Gen->Circle(px, py, pt);

			TLorentzVector Ds_vec(px, py, pz, sqrt(pt*pt + pz*pz + pow(Const::fMDs, 2)));
			Ds_vec.Boost(S.BoostVector());	//parent lab frame

		std::cout << "H2" << std::endl;
			vProducts = Nu_->ProductionPS(Amplitude::_CharmT, Ds_vec);
		std::cout << "H3" << std::endl;

			if (vProducts.at(0)->Theta() < Th0)	//neutrino
				hCharm->Fill(vProducts.at(0)->Energy());

		std::cout << "H4" << std::endl;
			TLorentzVector Tau_vec = vProducts.at(1)->FourVector();

		std::cout << "H5" << std::endl;
			for (unsigned int i = 0; i < 4; ++i)
			{
				Amplitude::Channel Name;
				TH1D* hFill;
				switch (i)
				{
					case 0:
		std::cout << "C0" << std::endl;
						Name = Amplitude::_TauET;
						hFill = hTauE;
						break;
					case 1:
		std::cout << "C1" << std::endl;
						Name = Amplitude::_TauMT;
						hFill = hTauM;
						break;
					case 2:
		std::cout << "C2" << std::endl;
						Name = Amplitude::_TauPI;
						hFill = hPion;
						break;
					case 3:
		std::cout << "C3" << std::endl;
						Name = Amplitude::_Tau2PI;
						hFill = h2Pion;
						break;
					default:
						break;
				}

		std::cout << "H6" << std::endl;
				vProducts.clear();
				vProducts = NuB->ProductionPS(Name, Tau_vec);

		std::cout << "H7" << std::endl;
				if (vProducts.at(0)->Theta() < Th0)	//neutrino
					hFill->Fill(vProducts.at(0)->Energy());
			}

		std::cout << "H8" << std::endl;
			++nIter;
		}
	}
	//		cc	pC	fDs	Br	Tot evts  
	double SF = 12.1e-3 / 331.4 * 0.077 * 0.0548 / double(nMAX);

	hCharm->Scale(SF);
	hTauE->Scale(SF * 0.1785);	//tau->e (17.85 %)
	hTauM->Scale(SF * 0.1736);	//tau->mu (17.36 %)
	hPion->Scale(SF * 0.1082);	//tau->pi (10.82 %)
	h2Pion->Scale(SF * 0.2551);	//tau->2pi (25.62 %)	Phase space only!!

	hTotal_->Add(hCharm);

	hTotalB->Add(hTauE);
	hTotalB->Add(hTauM);
	hTotalB->Add(hPion);
	hTotalB->Add(h2Pion);

	FileOut_->cd();

	hTotal_->Write();
	hCharm->Write();

	FileOutB->cd();

	hTauE->Write();
	hTauM->Write();
	hPion->Write();
	h2Pion->Write();

	FileOut_->Close();
	FileOutB->Close();

	return 0;
}

void Usage(char* argv0)
{
	std::cout << "Description" << std::endl;
	std::cout << "Usage : " << std::endl;
	std::cout << argv0 << " [OPTIONS]" << std::endl;
	std::cout <<"\n  -t,  --target" << std::endl;
	std::cout << "\t\tThe element of the target (available 'H', 'C')" << std::endl;
	std::cout <<"\n  -o,  --output" << std::endl;
	std::cout << "\t\tOutput file" << std::endl;
	std::cout <<"\n  -h,  --help" << std::endl;
	std::cout << "\t\tPrint this message and exit" << std::endl;
}
