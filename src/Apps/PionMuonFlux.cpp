#include <iostream>
#include <fstream>
#include <vector>
#include <cstring>
#include <getopt.h>

#include "TFile.h"
#include "TH1.h"
#include "TRandom3.h"

#include "Tools.h"

double DecayProb(double Tau, double Energy, double Mass);
int main(int argc, char** argv)
{

	const struct option longopts[] = 
	{
		{"output", 	required_argument, 	0, 'o'},
		{"sigma", 	required_argument, 	0, 's'},
		{"mean", 	required_argument, 	0, 'm'},
		{"sterile", 	required_argument, 	0, 'S'},
		{"help", 	no_argument,	 	0, 'h'},
		{0,	0, 	0,	0},
	};

	int index; 
	int iarg = 0;
	opterr = 1;
	
	std::ofstream OutFile;
	
//Initialize variables
	
	double Sigma = 1.0, Mean = 0.0, Sterile = 0.35;
	TFile *OutF;
	while((iarg = getopt_long(argc,argv, "o:s:m:S:h", longopts, &index)) != -1)
	{
		switch(iarg)
		{
			case 'o':
				OutF = new TFile(optarg, "RECREATE");
				break;
			case 's':
				Sigma = strtod(optarg, NULL);
				break;
			case 'm':
				Mean = strtod(optarg, NULL);
				break;
			case 'S':
				Sterile = strtod(optarg, NULL);
				break;
			case 'h':
				std::cout << "Description" << std::endl;
				std::cout << "Usage : " << std::endl;
				std::cout << "name [OPTIONS]" << std::endl;
				std::cout <<"\n  -o,  --output" << std::endl;
				std::cout << "\t\tROOT file" << std::endl;
				std::cout <<"\n  -h,  --help" << std::endl;
				std::cout << "\t\tPrint this message and exit" << std::endl;
				return 1;
			default:
				return 1;
		}
	}

	TRandom3 *RanGen = new TRandom3(0);
	TH1D * PionSpectrum = new TH1D("pion", "pion spectrum", 500, 0, 5);
	TH1D * pPion = new TH1D("ppion", "pion spectrum", 500, 0, 5);
	TH1D * KaonSpectrum = new TH1D("kaon", "kaon spectrum", 500, 0, 5);
	TH1D * pKaon = new TH1D("pkaon", "kaon spectrum", 500, 0, 5);
	TH1D * piNeutSpectrum = new TH1D("pineut", "pionneutrino spectrum", 500, 0, 5);
	TH1D * kaNeutSpectrum = new TH1D("kaneut", "kaonneutrino spectrum", 500, 0, 5);
	TH1D * SterileSpectrum = new TH1D("sterile", "neutrino spectrum", 500, 0, 5);
	
	double Mp = Const::fMPion;
	double Mk = Const::fMKaon;
	double Pp, Pk;
	double Ev, Ep, Ek;
	for (int i = 0; i < 1000000; ++i)
	{
		Pp = RanGen->Gaus(2,1);
		Pk = RanGen->Gaus(2,1);
		Ep = sqrt(Mp*Mp+Pp*Pp);
		Ek = sqrt(Mk*Mk+Pk*Pk);

		pPion->Fill(Pp);
		pKaon->Fill(Pk);
		PionSpectrum->Fill(Ep, DecayProb(26.03e-9, Ep, Mp));
		KaonSpectrum->Fill(Ek, DecayProb(12.385e-9, Ek, Mk));
	}

	double Mu = Const::fMMuon;
	double Cos, MM, EM;
	std::cout << "PionRest frame: " << (Mp*Mp - Mu*Mu)/(2*Mp) << std::endl;
	std::cout << "KaonRest frame: " << (Mk*Mk - Mu*Mu)/(2*Mk) << std::endl;
	for (int i = 0; i < 1000000; ++i)
	{
		Cos = RanGen->Uniform(-1,1);
		if (RanGen->Rndm() < 0.95)
		{
			MM = Mp;
			EM = PionSpectrum->GetRandom();
			Ev = (MM*MM - Mu*Mu)/(2*MM*MM)*(EM + Cos*sqrt(EM*EM - MM*MM));
			piNeutSpectrum->Fill(Ev);
		}
		else
		{
			MM = Mk;
			EM = KaonSpectrum->GetRandom();
			Ev = (MM*MM - Mu*Mu)/(2*MM*MM)*(EM + Cos*sqrt(EM*EM - MM*MM));
			kaNeutSpectrum->Fill(Ev);
		}
	}

	pPion->Write();
	pKaon->Write();
	PionSpectrum->Write();
	KaonSpectrum->Write();
	piNeutSpectrum->Write();
	kaNeutSpectrum->Write();

	piNeutSpectrum->Scale(Kine::ShrockFactor(Mp, Mu, Sterile));
	kaNeutSpectrum->Scale(Kine::ShrockFactor(Mk, Mu, Sterile));

	SterileSpectrum->Add(piNeutSpectrum);
	SterileSpectrum->Add(kaNeutSpectrum);
	SterileSpectrum->Write();

	//Create object from classes
	//e.g.	Decay * SuperGamma = new Decay(M_Sterile, U_e, U_m, U_t);

	//Main body

	//Garbage collection

	OutF->Close();

	return 0;
}
	
double DecayProb(double Tau, double Energy, double Mass)
{
	Tau *= Const::fS2GeV;
	double Length = 50*Const::fM2GeV;
	double Lorentz = Mass/sqrt(Energy*Energy-Mass*Mass);
	return exp(-Length*Lorentz/Tau);
}
