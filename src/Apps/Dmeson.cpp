#include <iostream>
#include <fstream>
#include <cstring>
#include <vector>
#include <getopt.h>

#include "Hadron.h"

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
	return n*(1-xf) - b * pt*pt;
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

	while((iarg = getopt_long(argc,argv, "r:s:t:o:h", longopts, &index)) != -1)
	{
		switch(iarg)
		{
			case 'r':
				OutF = new TFile(optarg, "RECREATE");
				break;
			case 's':
				SE = strtod(optarg, NULL);
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

	double nE, nPt, nM, nTh;
	double tE, tPt, tM, tTh;
	double vE, vPt, vM, vTh;
	double uE, uPt, uM, uTh;
	double DE, DPt, Dth;
	double dE, dPt, dth;
	tMeson->Branch("dE",	&dE,	"fdE/D");
	tMeson->Branch("dPt",	&dPt,	"fdPt/D");
	tMeson->Branch("dth",	&dth,	"fdth/D");
	tMeson->Branch("DE",	&DE,	"fDE/D");
	tMeson->Branch("DPt",	&DPt,	"fDPt/D");
	tMeson->Branch("Dth",	&Dth,	"fDth/D");
	//tMeson->Branch("n1", &n1E, "fn1E/D");
	//tMeson->Branch("n1prop", &n1E, "fn1E/D");
	//tMeson->Branch("t1", &t1E, "ft1E/D");
	/*
	tMeson->Branch("tE" , &tE , "ftE/D");
	tMeson->Branch("tPt", &tPt, "ftPt/D");
	tMeson->Branch("tM" , &tM , "ftM/D");
	tMeson->Branch("tTh", &tTh, "ftTheta/D");
	tMeson->Branch("vE" , &vE , "fvE/D");
	tMeson->Branch("vPt", &vPt, "fvPt/D");
	tMeson->Branch("vM" , &vM , "fvM/D");
	tMeson->Branch("vTh", &vTh, "fvTheta/D");
	tMeson->Branch("uE" , &uE , "fuE/D");
	tMeson->Branch("uPt", &uPt, "fuPt/D");
	tMeson->Branch("uM" , &uM , "fuM/D");
	tMeson->Branch("uTh", &uTh, "fuTheta/D");
	*/

	TRandom3 *Gen = new TRandom3(0);

	double mp = Const::fMProton;
	double mc = 1.27;		//c quark
	double mDs = 1.96847;		//Ds meson
	double mt = 1.77692;		//tau


	double Eb = 800;
	TLorentzVector Beam(0, 0, sqrt(Eb*Eb - mp*mp), Eb);
	TLorentzVector Targ(0, 0, 0, mp);
	TLorentzVector S = Beam+Targ;
	double sqrts = S.M();		//CM energy

	double psint, pcost;
	double masscharm[4] = {mDs, mDs, mDs, mDs};
	double massdecay[2] = {1.77692, 0.0};
	double massscatt[2] = {1.77692, mp};
	TGenPhaseSpace event;
	
	//find max of Ds param
	double ptmax = sqrts;
	double maxF = lDparam(-1, 0);
	double minF = lDparam(0.999, ptmax);

	for (unsigned int i = 0; i < 1e4; ++i)
	{
		double pt = Gen->Uniform(0, ptmax);
		double xf = Gen->Uniform(-1.0, 1.0);
		if (Gen->Uniform(minF, maxF) < lDparam(xf, pt))
		{
			std::cout << xf << "\t" << pt << "\t" << lDparam(xf, pt) << std::endl;
			double px, py;
			Gen->Circle(px, py, pt);
			//scale properly
			double pz = sqrts*xf*0.5;

			TLorentzVector Ds_vec(px, py, pz, sqrt(pt*pt + pz*pz - mDs*mDs));
			dxf = xf;
			dE = Ds_vec.E();
			dPt = Ds_vec.Pt();
			dth = Ds_vec.Theta();

			Ds_vec.Boost(S.BoostVector());
			DE = Ds_vec.E();
			DPt = Ds_vec.Pt();
			Dth = Ds_vec.Theta();

			tMeson->Fill();
		}
		
		/*
		event.SetDecay(S, 2, masscharm);
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
