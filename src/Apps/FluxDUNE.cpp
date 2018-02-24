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

	TH1D * hTotal = new TH1D("htotal", "total",  100, 0, 20);
	TH1D * hCharm = new TH1D("hcharm", "charm",  100, 0, 20);
	TH1D * hTauE  = new TH1D("htaue",  "taue",   100, 0, 20);
	TH1D * hTauM  = new TH1D("htaum",  "taum",   100, 0, 20);
	TH1D * hPion  = new TH1D("hpion",  "pion",   100, 0, 20);
	TH1D * h2Pion = new TH1D("h2pion", "2 pion", 100, 0, 20);

	TH1D * pTotal = new TH1D("ptotal", "total geo", 100, 0, 20);
	TH1D * pCharm = new TH1D("pcharm", "charm geo", 100, 0, 20);
	TH1D * pTauE  = new TH1D("ptaue",  "taue geo",  100, 0, 20);
	TH1D * pTauM  = new TH1D("ptaum",  "taum geo",  100, 0, 20);
	TH1D * pPion  = new TH1D("ppion",  "pion geo",  100, 0, 20);
	TH1D * p2Pion = new TH1D("p2pion", "2pion geo", 100, 0, 20);

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
	double massdecay0[2] = {0.0, mt};
	double massdecay1[3] = {0.0, 0.0, mm};
	double massdecay2[3] = {0.0, 0.0, me};
	double massdecay3[2] = {0.0, mpi};
	double massdecay4[3] = {0.0, mpi, mpi0};

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
		if (Gen->Uniform(0, pow(10,maxF)) < pow(10, lDparam(xf, pt)))
		{
			double px, py;
			Gen->Circle(px, py, pt);
			//scale properly
			double pz = sqrts*xf*0.5;

			TLorentzVector Ds_vec(px, py, pz, sqrt(pt*pt + pz*pz + mDs*mDs));
			Ds_vec.Boost(S.BoostVector());

			Decay0 = true;
			event.SetDecay(Ds_vec, 2, massdecay0);
			event.Generate();

			TLorentzVector tau_vec = *(event.GetDecay(1));
			TLorentzVector nut_vec = *(event.GetDecay(0));

			pCharm->Fill(nut_vec.E());
			if (nut_vec.Theta() < Th0)
				hCharm->Fill(nut_vec.E());

			unsigned int nChannel = 0;
			while (++nChannel < 5)
			{
				double Weight;
				TLorentzVector tau_cm(0, 0, 0, mt);
				switch(nChannel)
				{
					case 1:		//tau -> nu_t nu_e e
						Space->SetParent("TauE");
						Space->TauChannel();
						Space->SetUt(1.0);

						event.SetDecay(tau_cm, 3, massdecay1);
						Weight = 1.0;
						while (Gen->Rndm() < Weight)	//works in CM only
						{
							event.Generate();
							Space->SetEnergyX(event.GetDecay(1)->E());
							Space->SetEnergyY(event.GetDecay(2)->E());
							Weight = Space->ddGamma()/Space->MaxGamma();
						}
	
						nut_vec = *(event.GetDecay(0));
						nut_vec.Boost(tau_vec.BoostVector());
		
						pTauE->Fill(nut_vec.E());
						if (nut_vec.Theta() < Th0)
							hTauE->Fill(nut_vec.E());
						break;

					case 2:		//tau -> nu_t nu_m mu
						Space->SetParent("TauM");
						Space->TauChannel();
						Space->SetUt(1.0);
		
						event.SetDecay(tau_cm, 3, massdecay2);
						Weight = 1.0;
						while (Gen->Rndm() < Weight)	//works in CM only
						{
							event.Generate();
							Space->SetEnergyX(event.GetDecay(1)->E());
							Space->SetEnergyY(event.GetDecay(2)->E());
							Weight = Space->ddGamma()/Space->MaxGamma();
						}
		
						nut_vec = *(event.GetDecay(0));
						nut_vec.Boost(tau_vec.BoostVector());
		
						pTauM->Fill(nut_vec.E());
						if (nut_vec.Theta() < Th0)
							hTauM->Fill(nut_vec.E());
						break;
	
					case 3:		//tau -> nu_t pion
						event.SetDecay(tau_vec, 2, massdecay3);
						event.Generate();

						nut_vec = *(event.GetDecay(0));

						pPion->Fill(nut_vec.E());
						if (nut_vec.Theta() < Th0)
							hPion->Fill(nut_vec.E());
						break;

					case 4:		//tau -> nu_t pi pi0
						event.SetDecay(tau_vec, 3, massdecay3);
						event.Generate();

						nut_vec = *(event.GetDecay(0));

						p2Pion->Fill(nut_vec.E());
						if (nut_vec.Theta() < Th0)
							h2Pion->Fill(nut_vec.E());
						break;

					default:
						break;
				}
			}
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

	pCharm->Scale(SF);
	pTauE->Scale(SF * 0.1785);
	pTauM->Scale(SF * 0.1736);
	pPion->Scale(SF * 0.1082);
	p2Pion->Scale(SF * 0.2551);

	hTotal->Add(hCharm);
	hTotal->Add(hTauE);
	hTotal->Add(hTauM);
	hTotal->Add(hPion);
	hTotal->Add(h2Pion);

	pTotal->Add(pCharm);
	pTotal->Add(pTauE);
	pTotal->Add(pTauM);
	pTotal->Add(pPion);
	pTotal->Add(p2Pion);

	hTotal->Write();
	hCharm->Write();
	hTauE->Write();
	hTauM->Write();
	hPion->Write();
	h2Pion->Write();

	pTotal->Write();
	pCharm->Write();
	pTauE->Write();
	pTauM->Write();
	pPion->Write();
	p2Pion->Write();

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
