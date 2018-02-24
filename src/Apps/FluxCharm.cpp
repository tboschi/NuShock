#include <iostream>
#include <fstream>
#include <cstring>
#include <vector>
#include <getopt.h>

//#include "Hadron.h"
#include "Tools.h"
#include "ThreeBody.h"

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

	std::string sProb, sTarg("H");
	std::ofstream OutFile;
	std::string NuEFile, NuMFile;
	TFile *OutF;
	double SE = 1000;
	double Eb = 800;
	unsigned int nMAX = 1e5;

	while((iarg = getopt_long(argc,argv, "r:s:E:t:o:I:A:B:h", longopts, &index)) != -1)
	{
		switch(iarg)
		{
			case 'r':
				OutF = new TFile(optarg, "RECREATE");
				break;
			case 's':
				SE = strtod(optarg, NULL);
				break;
			case 'E':
				Eb = strtod(optarg, NULL);
				break;
			case 't':
				sTarg.assign(optarg);
				break;
			case 'o':
				OutFile.open(optarg);
				break;
			case 'A':
				NuEFile.assign(optarg);
				break;
			case 'B':
				NuMFile.assign(optarg);
				break;
			case 'I':
				nMAX = strtol(optarg, NULL, 10);
				break;
			case 'h':
				Usage(argv[0]);
				return 1;
			default:
				break;
		}
	}

	std::ostream &Out = (OutFile.is_open()) ? OutFile : std::cout;

	TH1D * hCharmE = new TH1D("hcharme", "charm",  100, 0, 20);
	TH1D * hCharmM = new TH1D("hcharmm", "charm",  100, 0, 20);

	double Th0 = atan2(5, 575);

	double N1E, N1Pt, N1th;
	double N2E, N2Pt, N2th;
	double N3E, N3Pt, N3th;
	bool Decay0, Decay1, Decay2, Decay3;

	TRandom3 *Gen = new TRandom3(0);

	double mp = Const::fMProton;
	double mc = Const::fMQuarkC;		//c quark
	double mDs = Const::fMDs;		//Ds meson
	double mt = Const::fMTau;
	double mm = Const::fMMuon;
	double me = Const::fMElectron;
	double mpi = Const::fMPion;
	double mpi0 = Const::fMPion0;

	TLorentzVector Beam(0, 0, sqrt(Eb*Eb - mp*mp), Eb);
	TLorentzVector Targ(0, 0, 0, mp);
	TLorentzVector S = Beam+Targ;
	double sqrts = S.M();		//CM energy

	double psint, pcost;
	double massdecayE[2] = {0.0, me};
	double massdecayM[2] = {0.0, mm};

	TGenPhaseSpace event;
	
	//find max of Ds param
	double ptmax = sqrts;
	double maxF = lDparam(-1, 0);
	double minF = lDparam(1, ptmax);

	ThreeBody * Space = new ThreeBody("");

	unsigned int nIter = 0;
	while (nIter < nMAX)
	{
		double pt = Gen->Uniform(0, ptmax);
		double xf = Gen->Uniform(-1.0, 1.0);
		//if (Gen->Uniform(minF, maxF) < lDparam(xf, pt))
		TLorentzVector nut_vec;
		if (Gen->Uniform(0, pow(10,maxF)) < pow(10, lDparam(xf, pt)))
		{
			double px, py;
			Gen->Circle(px, py, pt);
			//scale properly
			double pz = sqrts*xf*0.5;

			TLorentzVector Ds_vec(px, py, pz, sqrt(pt*pt + pz*pz + mDs*mDs));
			Ds_vec.Boost(S.BoostVector());

			event.SetDecay(Ds_vec, 2, massdecayE);
			event.Generate();
			nut_vec = *(event.GetDecay(0));

			if (nut_vec.Theta() < Th0)
				hCharmE->Fill(nut_vec.E());

			event.SetDecay(Ds_vec, 2, massdecayM);
			event.Generate();
			nut_vec = *(event.GetDecay(0));

			if (nut_vec.Theta() < Th0)
				hCharmM->Fill(nut_vec.E());

			++nIter;
		}
	}
	//		cc	pC	fDs	Tot evts  
	double SF = 12.1e-3 / 331.4 * 0.077 / double(nMAX);

	hCharmE->Scale(SF * 8.3e-5);
	hCharmM->Scale(SF * 5.5e-3);

	TFile FE(NuEFile.c_str(), "RECREATE");
	hCharmE->Write("hcharm");
	FE.Close();

	TFile FM(NuMFile.c_str(), "RECREATE");
	hCharmM->Write("hcharm");
	FM.Close();

//	OutF->Close();
//	if (OutFile.is_open())
//		OutFile.close();

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
