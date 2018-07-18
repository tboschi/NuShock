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
#include "TH2D.h"
#include "TRandom3.h"
#include "TGenPhaseSpace.h"

void Usage(char* argv0);
double Ds(double xf, double pt)
{
	//from 1708.08700, 250GeV proton beam E796
	double b = 1.08;
	double n = 6.1;
	//return n*(1-std::abs(xf)) - b * pt*pt;
	return std::exp(n * std::log(1 - std::abs(xf)) - b * pow(pt, 2));
}

int main(int argc, char** argv)
{
	const struct option longopts[] = 
	{
		{"nue", 	required_argument, 	0, 'e'},
		{"numu", 	required_argument, 	0, 'u'},
		{"help", 	no_argument,	 	0, 'h'},
		{0,	0, 	0,	0},
	};

	int index; 
	int iarg = 0;
	opterr = 1;

	//std::string sProb, sTarg("H"), OutName;
	std::string OutName, DetConfig, NuEFile, NuMuFile;
	std::ofstream OutFile;
	double BeamE = 800;	//DONUT
	unsigned int Nevent = 1e5;
	double NuMass = 0.0;	//neutrino mass

	while((iarg = getopt_long(argc,argv, "r:s:m:E:t:o:I:d:e:u:h", longopts, &index)) != -1)
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
				Nevent = strtol(optarg, NULL, 10);
				break;
			case 'd':
				DetConfig.assign(optarg);
				break;
			case 'e':
				NuEFile.assign(optarg);
				break;
			case 'u':
				NuMuFile.assign(optarg);
				break;
			case 'h':
				Usage(argv[0]);
				return 1;
			default:
				break;
		}
	}

	std::ostream &Out = (OutFile.is_open()) ? OutFile : std::cout;

	std::string Out_ = OutName + "_0.root";
	std::string OutB = OutName + "_B.root";

	TFile *FileOut_ = new TFile(Out_.c_str(), "RECREATE");
	TFile *FileOutB = new TFile(OutB.c_str(), "RECREATE");

	TH1D * hTotal_ = new TH1D("htotal1", "total",  100, 0, 20);
	TH1D * hTotalB = new TH1D("htotal2", "total",  100, 0, 20);

	//neutrino
	TH1D * hCharmE = new TH1D("hcharme", "charm",  100, 0, 20);
	TH1D * hCharmM = new TH1D("hcharmm", "charm",  100, 0, 20);
	TH1D * hCharmT = new TH1D("hcharm", "charm",  100, 0, 20);

	//antineutrino
	TH1D * hTauE  = new TH1D("htaue",  "taue",   100, 0, 20);
	TH1D * hTauM  = new TH1D("htaum",  "taum",   100, 0, 20);
	TH1D * hPion  = new TH1D("hpion",  "pion",   100, 0, 20);
	TH1D * h2Pion = new TH1D("h2pion", "2 pion", 100, 0, 20);


	hTotal_->SetDirectory(0);
	hCharmT->SetDirectory(FileOut_);

	hCharmE->SetDirectory(0);
	hCharmM->SetDirectory(0);

	hTotalB->SetDirectory(0);
	hTauE->SetDirectory(FileOutB);
	hTauM->SetDirectory(FileOutB);
	hPion->SetDirectory(FileOutB);
	h2Pion->SetDirectory(FileOutB);

	//generous angular acceptance of detector
	Detector *TheBox = new Detector(DetConfig);
	double Th0 = TheBox->AngularAcceptance();

	TRandom3 *Gen = new TRandom3(0);

	TLorentzVector Beam(0, 0, sqrt(pow(BeamE, 2) - pow(Const::fMProton, 2)), BeamE);
	TLorentzVector Targ(0, 0, 0, Const::fMProton);
	TLorentzVector S = Beam+Targ;
	double ptmax = S.M();		//CM energy

	Neutrino *Nu_ = new Neutrino(NuMass, Neutrino::Dirac | Neutrino::Left );
	Neutrino *NuB = new Neutrino(NuMass, Neutrino::Dirac | Neutrino::Right | Neutrino::Antiparticle);
	
	//Normalisation
	double SF = 12.1e-3 / 331.4 * 0.077 / Nevent;

	//std::vector<Particle*> vProductDs, vProductTau;
	//std::vector<Particle*>::iterator iP;
	std::vector<Particle> vProductDs, vProductTau;
	std::vector<Particle>::iterator iP;

	unsigned int DecayCount = 0, InNDCount = 0;
	for (unsigned int ID = 0; ID < Nevent; ++ID)
	{
		double pt, xf;
		do
		{
			pt = Gen->Uniform(0, ptmax);
			xf = Gen->Uniform(-1.0, 1.0);
		}
		while (Gen->Rndm() > Ds(xf, pt));

		double px, py;
		double pz = ptmax * xf / 2.0;
		Gen->Circle(px, py, pt);

		TLorentzVector Ds_vec(px, py, pz, sqrt(pt*pt + pz*pz + pow(Const::fMDs, 2)));
		Ds_vec.Boost(S.BoostVector());	//parent lab frame

		//Ds decay into electrons
		if (!NuEFile.empty())
		{
			vProductDs = Nu_->ProductionPS(Amplitude::_CharmE, Ds_vec);
			if (vProductDs.at(0).Theta() < Th0)
				hCharmE->Fill(vProductDs.at(0).Energy(), SF * 8.3e-5);

			//for (iP = vProductDs.begin(); iP != vProductDs.end(); ++iP)
			//	delete *iP;
			vProductDs.clear();
		}

		//Ds decay into muons
		if (!NuMuFile.empty())
		{
			vProductDs = Nu_->ProductionPS(Amplitude::_CharmM, Ds_vec);
			if (vProductDs.at(0).Theta() < Th0)
				hCharmM->Fill(vProductDs.at(0).Energy(), SF * 5.5e-3);

			//for (iP = vProductDs.begin(); iP != vProductDs.end(); ++iP)
			//	delete *iP;
			vProductDs.clear();
		}

		//Ds decay into taus
		vProductDs = Nu_->ProductionPS(Amplitude::_CharmT, Ds_vec);
		//std::cout << "fluxds ";
		//for (unsigned int i = 0; i < vProductDs.size(); ++i)
		//	std::cout << vProductDs.at(i) << "\t";
		//std::cout << std::endl;

		if(!vProductDs.size())
			continue;
		else
			++DecayCount;

		if (vProductDs.at(0).Theta() <= Th0)	//neutrino
			hCharmT->Fill(vProductDs.at(0).Energy(), SF * 0.0548);
		else
			++InNDCount;

		//tau decay from Ds
		TLorentzVector Tau_vec(vProductDs.at(1).FourVector());
		for (unsigned int i = 0; i < 4; ++i)
		{
			Amplitude::Channel Name;
			TH1D* hFill;
			double Br;
			switch (i)
			{
				case 0:
					Name = Amplitude::_TauET;
					hFill = hTauE;
					Br = SF * 0.0548 * 0.1785;	//tau->e (17.85 %)
					break;
				case 1:
					Name = Amplitude::_TauMT;
					hFill = hTauM;
					Br = SF * 0.0548 * 0.1736;	//tau->mu (17.36 %)
					break;
				case 2:
					Name = Amplitude::_TauPI;
					hFill = hPion;
					Br = SF * 0.0548 * 0.1082;	//tau->pi (10.82 %)
					break;
				case 3:
					Name = Amplitude::_Tau2PI;
					hFill = h2Pion;
					Br = SF * 0.0548 * 0.2551;	//tau->2pi (25.62 %)
					break;			//Phase space only!!
				default:
					break;
			}

			vProductTau = NuB->ProductionPS(Name, Tau_vec);
			//std::cout << "fluxtau ";
			//for (unsigned int i = 0; i < vProductTau.size(); ++i)
			//	std::cout << vProductTau.at(i) << "\t";
			//std::cout << std::endl;
			if (!vProductTau.size())
				continue;

			if (vProductTau.at(0).Theta() <= Th0)	//neutrino
				hFill->Fill(vProductTau.at(0).Energy(), Br);

			//std::cout << "delete fluxT ";
			//for (iP = vProductTau.begin(); iP != vProductTau.end(); ++iP)
			//{
			//	//std::cout << *iP << "\t";
			//	delete *iP;
			//}
			//std::cout << std::endl;
			vProductTau.clear();
		}

		//if (ID % 10000 == 0)	//saving
		if (false)	//saving
		{
			FileOut_->Write("", TObject::kOverwrite);
			FileOutB->Write("", TObject::kOverwrite);
		}

		//std::cout << "delete ds ";
		for (iP = vProductDs.begin(); iP != vProductDs.end(); ++iP)
		{
			//std::cout << *iP << "\t";
			//delete *iP;
		}
		//std::cout << std::endl;
		vProductDs.clear();
	}

	hTotal_->Add(hCharmT);

	hTotalB->Add(hTauE);
	hTotalB->Add(hTauM);
	hTotalB->Add(hPion);
	hTotalB->Add(h2Pion);

	FileOut_->cd();
	FileOut_->Write();

	hTotal_->Write("htotal");
	//hCharm->Write();

	FileOutB->cd();
	FileOutB->Write();

	hTotalB->Write("htotal");

	std::cout << "Ds meson decays are " << 100.0 * DecayCount / double(Nevent) << " %\n";
	std::cout << "Products in ND are " << 100.0 * (1.0 - InNDCount / double(Nevent)) << " %\n";
	std::cout  << "Neutrinos simulated " << hTotal_->GetEntries();
	std::cout << " (" << hTotal_->GetEntries()*100.0/double(Nevent) << " %)";
	std::cout << ", saved in " << FileOut_->GetName() << std::endl;
	std::cout  << "Antineutrinos simulated " << hTotalB->GetEntries();
	std::cout << " (" << hTotalB->GetEntries()*100.0/double(Nevent) << " %)";
	std::cout << ", saved in " << FileOutB->GetName() << std::endl;

	FileOut_->Close();
	FileOutB->Close();

	if (!NuEFile.empty())
	{
		FileOut_ = new TFile(NuEFile.c_str(), "UPDATE");
		hTotal_ = dynamic_cast<TH1D*>(FileOut_->Get("htotal"));
		hTotal_->Add(hCharmE);
		hTotal_->Write("", TObject::kOverwrite);
		hCharmE->Write("hcharm", TObject::kOverwrite);
		FileOut_->Close();
	}

	if (!NuEFile.empty())
	{
		FileOut_ = new TFile(NuMuFile.c_str(), "UPDATE");
		hTotal_ = dynamic_cast<TH1D*>(FileOut_->Get("htotal"));
		hTotal_->Add(hCharmE);
		hTotal_->Write("", TObject::kOverwrite);
		hCharmM->Write("hcharm", TObject::kOverwrite);
		FileOut_->Close();
	}

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
