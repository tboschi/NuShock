#include<iostream>
#include <fstream>

#include "Tools.h"
#include "TFile.h"
#include "TTree.h"

int main(int argc, char** argv)
{
	const struct option longopts[] = 
	{
		{"help", 	no_argument,	 	0, 'h'},
		{0,	0, 	0,	0},
	};

	int index; 
	int iarg = 0;
	opterr = 1;
	
	bool Back = false;
	bool Dirac = true;
	bool LNC = true;
	std::string base, outfile;
	
	while((iarg = getopt_long(argc,argv, "i:o:sbdmcv", longopts, &index)) != -1)
	{
		switch(iarg)
		{
			case 'i':
				base.assign(optarg);
				break;
			case 'o':
				outfile.assign(optarg);
				break;
			case 's':
				Back = false;
				break;
			case 'b':
				Back = true;
				break;
			case 'd':
				Dirac = true;
				break;
			case 'm':
				Dirac = false;
				break;
			case 'c':
				LNC = true;
				break;
			case 'v':
				LNC = false;
				break;
			default:
				break;
		}
	}
	//bool Back = true;

	TH1D* hmass = new TH1D("mass", "mass", 100, 0, 2);
	TH1D* hangl = new TH1D("angl", "angl", 100, 0, 60);
	TH1D* hcos0 = new TH1D("cos0", "cos0", 100, 0.95, 1);
	TH1D* hener = new TH1D("ener", "ener", 100, 0, 20);
	hmass->SetDirectory(0);
	hangl->SetDirectory(0);
	hener->SetDirectory(0);

	double bw_m = 2.0/100;
	double bw_a = 60.0/100.0;
	double bw_c = 0.05/100.0;
	double bw_e = 20.0/100.0;

	std::ofstream Out(outfile.c_str());

	if (Back)
	{
		std::string EE = base + "_E.root";
		std::string EM = base + "_M.root";
		std::string EB = base + "_B.root";
		double ew = (3.0+1.0);
		double mw = (240.0+79.0);
		double bw = (17.8+7.3);
		double tot = ew+mw+bw;
		ew /= tot;
		mw /= tot;
		bw /= tot;

		TFile *TF = new TFile(EE.c_str(), "OPEN");
		TTree *Data = dynamic_cast<TTree*> (TF->Get("Data"));

		double W, M_0, Angle, E_A;
		TBranch  *b_fW, *b_fAngle, *b_fEB, *b_fMass0;   //!
		//Data->SetBranchAddress("W",     &W,     &b_fW);
		Data->SetBranchStatus("*", 1);		// disable all branches
		Data->SetBranchAddress("TheB", &Angle, &b_fAngle);
		Data->SetBranchAddress("E_A", &E_A, &b_fEB);
		Data->SetBranchAddress("M_0",   &M_0,   &b_fMass0);

		unsigned int n = Data->GetEntries();
		for (unsigned int i = 0; i < n; ++i)
		{
			Data->GetEntry(i);
			hmass->Fill(M_0, ew);
			hcos0->Fill(cos(Angle), ew);
			hangl->Fill(Const::fDeg*Angle, ew);
			hener->Fill(E_A, ew);
		}

		TF->Close();
		TF = new TFile(EM.c_str(), "OPEN");
		Data = dynamic_cast<TTree*> (TF->Get("Data"));

		n = Data->GetEntries();
		Data->SetBranchStatus("*", 1);		// disable all branches
		Data->SetBranchAddress("TheB", &Angle, &b_fAngle);
		Data->SetBranchAddress("E_A", &E_A, &b_fEB);
		Data->SetBranchAddress("M_0",   &M_0,   &b_fMass0);
		
		for (unsigned int i = 0; i < n; ++i)
		{
			Data->GetEntry(i);
			hmass->Fill(M_0, mw);
			hcos0->Fill(cos(Angle), mw);
			hangl->Fill(Const::fDeg*Angle, mw);
			hener->Fill(E_A, bw);
		}

		TF->Close();
		TF = new TFile(EB.c_str(), "OPEN");
		Data = dynamic_cast<TTree*> (TF->Get("Data"));

		n = Data->GetEntries();
		Data->SetBranchStatus("*", 1);		// disable all branches
		Data->SetBranchAddress("TheB", &Angle, &b_fAngle);
		Data->SetBranchAddress("E_A", &E_A, &b_fEB);
		Data->SetBranchAddress("M_0",   &M_0,   &b_fMass0);
		
		for (unsigned int i = 0; i < n; ++i)
		{
			Data->GetEntry(i);
			hmass->Fill(M_0, bw);
			hcos0->Fill(cos(Angle), bw);
			hangl->Fill(Const::fDeg*Angle, bw);
			hener->Fill(E_A, bw);
		}

		hmass->Scale(1.0/bw_m);
		hangl->Scale(1.0/bw_a);
		hcos0->Scale(1.0/bw_c);
		hener->Scale(1.0/bw_e);
		//std::cout << "int " << hmass->Integral("WIDTH") << "\t" << hangl->Integral("WIDTH") << std::endl;
	}
	else	//signal
	{
		TFile *TF = new TFile(base.c_str(), "OPEN");
		TTree *Data = dynamic_cast<TTree*> (TF->Get("Data"));

		double W, M_0, Angle, E_A;
		bool P;
		TBranch  *b_fW, *b_fAngle, *b_fEB, *b_fMass0, *b_fP;   //!
		//Data->SetBranchAddress("W",     &W,     &b_fW);
		Data->SetBranchStatus("*", 1);		// disable all branches
		//Data->SetBranchAddress("R", &W, &b_fW);
		Data->SetBranchAddress("P", &P, &b_fP);
		Data->SetBranchAddress("W", &W, &b_fW);
		Data->SetBranchAddress("TheB", &Angle, &b_fAngle);
		Data->SetBranchAddress("E_0", &E_A, &b_fEB);
		Data->SetBranchAddress("M_0",   &M_0,   &b_fMass0);

		double tw;
		unsigned int n = Data->GetEntries();
		for (unsigned int i = 0; i < n; ++i)
		{
			Data->GetEntry(i);

			if (Dirac)
			{
				if (LNC)
					tw = W * P;
				else
					tw = W * !P;
			}
			else 	//Majorana
				tw = W / 2.0;

			hmass->Fill(M_0, tw);
			hcos0->Fill(cos(Angle), tw);
			hangl->Fill(Const::fDeg*Angle, tw);
			hener->Fill(E_A, tw);
		}

		TF->Close();
		hmass->Scale(1.0/bw_m);///n);
		hangl->Scale(1.0/bw_a);///n);
		hcos0->Scale(1.0/bw_c);///n);
		hener->Scale(1.0/bw_e);///n);
	}


	for (unsigned int i = 0; i < 100; ++i)
	{
		Out << hmass->GetBinCenter(i+1) << "\t";	//1
		Out << hmass->GetBinContent(i+1) << "\t";	//2
		Out << hangl->GetBinCenter(i+1) << "\t";	//3
		Out << hangl->GetBinContent(i+1) << "\t";	//4
		Out << hcos0->GetBinCenter(i+1) << "\t";	//5
		Out << hcos0->GetBinContent(i+1) << "\t";	//6
		Out << hener->GetBinCenter(i+1) << "\t";	//7
		Out << hener->GetBinContent(i+1) << "\t";	//8
		Out << std::endl;
	}

	return 0;
}
