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
	TFile *OutF;
	double SE = 1000;
	double Eb = 800;
	unsigned int nMAX = 1e5;

	while((iarg = getopt_long(argc,argv, "r:s:E:t:o:I:h", longopts, &index)) != -1)
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

	TH1D * hTotal = new TH1D("htotal", "htotal", 400, 0, 40);
	TH1D * hCharm = new TH1D("hcharm", "hcharm", 400, 0, 40);
	TH1D * hTauE = new TH1D("htaue", "htaue", 400, 0, 40);
	TH1D * hTauM = new TH1D("htaum", "htaum", 400, 0, 40);
	TH1D * hPion = new TH1D("hpion", "hpion", 400, 0, 40);

	TH1D * pTotal = new TH1D("ptotal", "total geo", 400, 0, 40);
	TH1D * pCharm = new TH1D("pcharm", "charm geo", 400, 0, 40);
	TH1D * pTauE = new TH1D("ptaue", "taue geo", 400, 0, 40);
	TH1D * pTauM = new TH1D("ptaum", "taum geo", 400, 0, 40);
	TH1D * pPion = new TH1D("ppion", "pion geo", 400, 0, 40);

	double Th0 = atan2(3.5, 574);

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
	while (nIter < nMAX)
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
			//std::cout << nIter << "\t";
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

			//if (Gen->Rndm() < 0.0548)
			if (true)		//every Ds decays
			{
				Decay0 = true;
				event.SetDecay(Ds_vec, 2, massdecay0);
				event.Generate();

				TLorentzVector tau_vec = *(event.GetDecay(1));
				TLorentzVector nut_vec = *(event.GetDecay(0));
				TE  = tau_vec.E();
				TPt = tau_vec.Pt();
				Tth = tau_vec.Theta();
				N0E  = nut_vec.E();
				N0Pt = nut_vec.Pt();
				N0th = nut_vec.Theta();

				hCharm->Fill(N0E);
				if (N0th < Th0)
					pCharm->Fill(N0E);

				double Branch = Gen->Rndm();
				if (Branch < 0.1785)		//tau->e (17.85 %)
				{
					Decay1 = true;
					Space->SetParent("TauE");
					Space->TauChannel();
					Space->SetUt(1.0);

					TLorentzVector tau_cm(0, 0, 0, mt);
					event.SetDecay(tau_cm, 3, massdecay1);
					double Weight = 1.0;
					while (Gen->Rndm() < Weight)	//works in CM only
					{
						event.Generate();
						Space->SetEnergyX(event.GetDecay(1)->E());
						Space->SetEnergyY(event.GetDecay(2)->E());
						Weight = Space->ddGamma()/Space->MaxGamma();
					}

					TLorentzVector nut_vec = *(event.GetDecay(0));
					nut_vec.Boost(tau_vec.BoostVector());
					N1E  = nut_vec.E();
					N1Pt = nut_vec.Pt();
					N1th = nut_vec.Theta();

					hTauE->Fill(N1E);
					if (N1th < Th0)
						pTauE->Fill(N1E);
				}
				else if (Branch < 0.3521)	//tau->mu (17.36 %)
				{
					Decay2 = true;
					Space->SetParent("TauM");
					Space->TauChannel();
					Space->SetUt(1.0);

					TLorentzVector tau_cm(0, 0, 0, mt);
					event.SetDecay(tau_cm, 3, massdecay2);
					double Weight = 1.0;
					while (Gen->Rndm() < Weight)	//works in CM only
					{
						event.Generate();
						Space->SetEnergyX(event.GetDecay(1)->E());
						Space->SetEnergyY(event.GetDecay(2)->E());
						Weight = Space->ddGamma()/Space->MaxGamma();
					}

					TLorentzVector nut_vec = *(event.GetDecay(0));
					nut_vec.Boost(tau_vec.BoostVector());
					N2E  = nut_vec.E();
					N2Pt = nut_vec.Pt();
					N2th = nut_vec.Theta();

					hTauM->Fill(N2E);
					if (N2th < Th0)
						pTauM->Fill(N2E);
				}
				else if (Branch < 0.4612)	//tau->pi (10.82 %)
				{
					Decay3 = true;

					event.SetDecay(tau_vec, 2, massdecay3);
					event.Generate();

					TLorentzVector nut_vec = *(event.GetDecay(0));
					N3E  = nut_vec.E();
					N3Pt = nut_vec.Pt();
					N3th = nut_vec.Theta();

					hPion->Fill(N3E);
					if (N3th < Th0)
						pPion->Fill(N3E);
				}
			}

			tMeson->Fill();
			++nIter;
		}
	}
	//		cc	pC	fDs	Tot evts  
	double SF = 12.1e-3 / 331.4 * 0.077 / double(nMAX);

	hCharm->Scale(SF);
	hTauE->Scale(SF);
	hTauM->Scale(SF);
	hPion->Scale(SF);

	pCharm->Scale(SF);
	pTauE->Scale(SF);
	pTauM->Scale(SF);
	pPion->Scale(SF);

	hTotal->Add(hCharm);
	hTotal->Add(hTauE);
	hTotal->Add(hTauM);
	hTotal->Add(hPion);

	pTotal->Add(pCharm);
	pTotal->Add(pTauE);
	pTotal->Add(pTauM);
	pTotal->Add(pPion);

	tMeson->Write();

	hTotal->Write();
	hCharm->Write();
	hTauE->Write();
	hTauM->Write();
	hPion->Write();

	pTotal->Write();
	pCharm->Write();
	pTauE->Write();
	pTauM->Write();
	pPion->Write();

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
