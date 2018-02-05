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
	TFile *OutF;
	double SE = 1000;
	double Eb = 800;

	while((iarg = getopt_long(argc,argv, "r:s:E:t:o:h", longopts, &index)) != -1)
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
			case 'h':
				Usage(argv[0]);
				return 1;
			default:
				break;
		}
	}

	std::ostream &Out = (OutFile.is_open()) ? OutFile : std::cout;
	TTree *tMeson = new TTree("tmeson", "Mesons");

	double dE, dPt, dth, dxf;
	double DE, DPt, Dth;
	double TE, TPt, Tth;
	//0 promtp nu tau (fomr Ds), 1 from tau->e, 2 fomr tau->mu, 3 from tau-> 1pi
	double N0E, N0Pt, N0th;
	double N1E, N1Pt, N1th;
	double N2E, N2Pt, N2th;
	double N3E, N3Pt, N3th;
	bool Decay0, Decay1, Decay2, Decay3;

	tMeson->Branch("dxf",	&dxf,	"fdxf/D");
	tMeson->Branch("dE",	&dE,	"fdE/D");
	tMeson->Branch("dPt",	&dPt,	"fdPt/D");
	tMeson->Branch("dth",	&dth,	"fdth/D");
	tMeson->Branch("DE",	&DE,	"fDE/D");
	tMeson->Branch("DPt",	&DPt,	"fDPt/D");
	tMeson->Branch("Dth",	&Dth,	"fDth/D");
	tMeson->Branch("TE",	&TE,	"fTE/D");
	tMeson->Branch("TPt",	&TPt,	"fTPt/D");
	tMeson->Branch("Tth",	&Tth,	"fTth/D");
	tMeson->Branch("N0E",	&N0E,	"fN0E/D");
	tMeson->Branch("N0Pt",	&N0Pt,	"fN0Pt/D");
	tMeson->Branch("N0th",	&N0th,	"fN0th/D");
	tMeson->Branch("N1E",	&N1E,	"fN1E/D");
	tMeson->Branch("N1Pt",	&N1Pt,	"fN1Pt/D");
	tMeson->Branch("N1th",	&N1th,	"fN1th/D");
	tMeson->Branch("N2E",	&N2E,	"fN2E/D");
	tMeson->Branch("N2Pt",	&N2Pt,	"fN2Pt/D");
	tMeson->Branch("N2th",	&N2th,	"fN2th/D");
	tMeson->Branch("N3E",	&N3E,	"fN3E/D");
	tMeson->Branch("N3Pt",	&N3Pt,	"fN3Pt/D");
	tMeson->Branch("N3th",	&N3th,	"fN3th/D");
	tMeson->Branch("Decay0", &Decay0, "bDecay0/O");
	tMeson->Branch("Decay1", &Decay1, "bDecay1/O");
	tMeson->Branch("Decay2", &Decay2, "bDecay2/O");
	tMeson->Branch("Decay3", &Decay3, "bDecay3/O");

	TRandom3 *Gen = new TRandom3(0);

	double mp = Const::fMProton;
	double mc = Const::fMQuarkC;		//c quark
	double mDs = Const::fMDs;		//Ds meson
	double mt = Const::fMTau;
	double mm = Const::fMMuon;
	double me = Const::fMElectron;
	double mpi = Const::fMPion;

	TLorentzVector Beam(0, 0, sqrt(Eb*Eb - mp*mp), Eb);
	TLorentzVector Targ(0, 0, 0, mp);
	TLorentzVector S = Beam+Targ;
	double sqrts = S.M();		//CM energy

	double psint, pcost;
	double massdecay0[2] = {0.0, mt};
	double massdecay1[3] = {0.0, 0.0, mm};
	double massdecay2[3] = {0.0, 0.0, me};
	double massdecay3[2] = {0.0, mpi};

	TGenPhaseSpace event;
	
	//find max of Ds param
	double ptmax = sqrts;
	double maxF = lDparam(-1, 0);
	double minF = lDparam(1, ptmax);

	ThreeBody * Space = new ThreeBody("");

	unsigned int nIter = 0;
	while (nIter < 1e6)
	{
		Decay0 = false;
		Decay1 = false;
		Decay2 = false;
		Decay3 = false;
		TE  = -1.0;
                TPt = -1.0;
                Tth = -1.0;
                N0E  = -1.0;
                N0Pt = -1.0;
                N0th = -1.0;
                N1E  = -1.0;
                N1Pt = -1.0;
                N1th = -1.0;
                N2E  = -1.0;
                N2Pt = -1.0;
                N2th = -1.0;
                N3E  = -1.0;
                N3Pt = -1.0;
                N3th = -1.0;

		double pt = Gen->Uniform(0, ptmax);
		double xf = Gen->Uniform(-1.0, 1.0);
		//if (Gen->Uniform(minF, maxF) < lDparam(xf, pt))
		if (Gen->Uniform(0, pow(10,maxF)) < pow(10, lDparam(xf, pt)))
		{
			std::cout << xf << "\t" << pt << "\t" << lDparam(xf, pt) << std::endl;
			double px, py;
			Gen->Circle(px, py, pt);
			//scale properly
			double pz = sqrts*xf*0.5;

			TLorentzVector Ds_vec(px, py, pz, sqrt(pt*pt + pz*pz + mDs*mDs));
			dxf = xf;
			dE = Ds_vec.E();
			dPt = Ds_vec.Pt();
			dth = Ds_vec.Theta();

			Ds_vec.Boost(S.BoostVector());
			DE = Ds_vec.E();
			DPt = Ds_vec.Pt();
			Dth = Ds_vec.Theta();

			if (Gen->Rndm() < 0.0548)
			{
				Decay0 = true;
				event.SetDecay(Ds_vec, 2, massdecay0);
				event.Generate();

				TLorentzVector tau_vec = *(event.GetDecay(0));
				TLorentzVector nut_vec = *(event.GetDecay(1));
				TE  = tau_vec.E();
				TPt = tau_vec.Pt();
				Tth = tau_vec.Theta();
				N0E  = nut_vec.E();
				N0Pt = nut_vec.Pt();
				N0th = nut_vec.Theta();

				double Branch = Gen->Rndm();
				if (Branch < 0.1785)		//tau->e (17.85 %)
				{
					Decay1 = true;
					Space->SetParent("TauE");

					event.SetDecay(tau_vec, 3, massdecay1);
					double Weight = 1.0;
					while (Gen->Rndm() < Weight)
					{
						event.Generate();
						Space->SetEnergyX(event.GetDecay(0)->E());	//don't remember this!!
						Space->SetEnergyY(event.GetDecay(1)->E());
						Weight = Space->ddGamma()/Space->MaxGamma();
					}

					TLorentzVector nut_vec = *(event.GetDecay(0));
					N1E  = nut_vec.E();
					N1Pt = nut_vec.Pt();
					N1th = nut_vec.Theta();
				}
				else if (Branch < 0.3521)	//tau->mu (17.36 %)
				{
					Decay2 = true;
					Space->SetParent("TauM");

					event.SetDecay(tau_vec, 3, massdecay2);
					double Weight = 1.0;
					while (Gen->Rndm() < Weight)
					{
						event.Generate();
						Space->SetEnergyX(event.GetDecay(0)->E());	//don't remember this!!
						Space->SetEnergyY(event.GetDecay(1)->E());
						Weight = Space->ddGamma()/Space->MaxGamma();
					}

					TLorentzVector nut_vec = *(event.GetDecay(0));
					N2E  = nut_vec.E();
					N2Pt = nut_vec.Pt();
					N2th = nut_vec.Theta();
				}
				else if (Branch < 0.4612)	//tau->mu (17.36 %)
				{
					Decay3 = true;

					event.SetDecay(tau_vec, 2, massdecay3);
					event.Generate();

					TLorentzVector nut_vec = *(event.GetDecay(0));
					N3E  = nut_vec.E();
					N3Pt = nut_vec.Pt();
					N3th = nut_vec.Theta();
				}
			}

			tMeson->Fill();
			++nIter;
		}
		
		/*
		event.Generate();
		TLorentzVector D1vec = *(event.GetDecay(0));
		D1E = D1vec.E();

		event.SetDecay(D1vec, 2, massdecay);
		event.Generate();
		n1E = event.GetDecay(1)->E();

		TLorentzVector nuvec = *(event.GetDecay(1))+Targ;
		event.SetDecay(nuvec, 2, massscatt);
		event.Generate();
		t1E = event.GetDecay(0)->E();
		*/
	}

	tMeson->Write();
	OutF->Close();
	if (OutFile.is_open())
		OutFile.close();

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
