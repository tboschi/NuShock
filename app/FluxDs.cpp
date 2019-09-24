#include <iostream>
#include <fstream>
#include <cstring>
#include <vector>
#include <getopt.h>

#include "tools.h"
#include "detector.h"
#include "physics.h"
#include "flux.h"

#include "TKey.h"
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
	std::string outName, detConfig, nuEFile, nuMFile;
	std::ofstream outFile;
	double beamE = 800;	//(GeV)	//800 GeV for DONUT for comparison
	int Nevent = 1e5;
	double mass = 0.0;	//neutrino mass in MeV

	while((iarg = getopt_long(argc,argv, "r:s:m:E:t:o:I:d:e:u:h", longopts, &index)) != -1)
	{
		switch(iarg)
		{
			case 'r':
				outName.assign(optarg);
				break;
			case 'm':
				mass = strtod(optarg, NULL);
				mass /= 1000.0;
				break;
			case 'E':
				beamE = strtod(optarg, NULL);
				break;
			case 'o':
				outFile.open(optarg);
				break;
			case 'I':
				Nevent = strtol(optarg, NULL, 10);
				break;
			case 'd':
				detConfig.assign(optarg);
				break;
			case 'e':
				nuEFile.assign(optarg);
				break;
			case 'u':
				nuMFile.assign(optarg);
				break;
			case 'h':
				Usage(argv[0]);
				return 1;
			default:
				break;
		}
	}

	std::ostream &out = (outFile.is_open()) ? outFile : std::cout;

	std::string out0 = outName + "_0.root";
	std::string outB = outName + "_B.root";

	TFile *FileOut0 = new TFile(out0.c_str(), "RECREATE");
	TFile *FileOutB = new TFile(outB.c_str(), "RECREATE");

	TH1D * hTotal0 = new TH1D("htotal1", "total",  100, 0, 20);
	TH1D * hTotalB = new TH1D("htotal2", "total",  100, 0, 20);

	//neutrino
	TH1D * hCharmE = new TH1D("hcharme", "charm",  100, 0, 20);
	TH1D * hCharmM = new TH1D("hcharmm", "charm",  100, 0, 20);
	TH1D * hCharmT = new TH1D("hcharm", "charm",  100, 0, 20);
	TH1D * pCharmT = new TH1D("pcharm", "charm",  100, 0, 20);

	//antineutrino
	TH1D * hTauE  = new TH1D("htaue",  "taue",   100, 0, 20);
	TH1D * hTauM  = new TH1D("htaum",  "taum",   100, 0, 20);
	TH1D * hPion  = new TH1D("hpion",  "pion",   100, 0, 20);
	TH1D * h2Pion = new TH1D("h2pion", "2 pion", 100, 0, 20);


	hTotal0->SetDirectory(0);
	hCharmT->SetDirectory(FileOut0);
	pCharmT->SetDirectory(FileOut0);

	hCharmE->SetDirectory(0);
	hCharmM->SetDirectory(0);

	hTotalB->SetDirectory(0);
	hTauE->SetDirectory(FileOutB);
	hTauM->SetDirectory(FileOutB);
	hPion->SetDirectory(FileOutB);
	h2Pion->SetDirectory(FileOutB);

	//generous angular acceptance of detector
	Detector *theBox = new Detector(detConfig);
	double th0 = theBox->AngularAcceptance();

	TRandom3 *mt = new TRandom3(0);

	TLorentzVector beam(0, 0, sqrt(pow(beamE, 2) - pow(Const::MProton, 2)), beamE);
	TLorentzVector targ(0, 0, 0, Const::MProton);
	TLorentzVector S = beam+targ;
	double ptmax = S.M();		//CM energy

	//neutrinos
	//Nu0 and NuB are to directly produce HNL flux from tau mixing
	//if mass == 0 then it is light neutrino flux
	//no difference between Majorana or Dirac
	//
	Neutrino Nu0(mass, Neutrino::Dirac | Neutrino::Left );
	Neutrino NuB(mass, Neutrino::Dirac | Neutrino::Right | Neutrino::Antiparticle);
	
	//Normalisation from open charm calculation
	//	    cc xsec /  pA xsec * fragmentation / Nevent
	double SF = 12.1e-3 / 331.4 * 0.077 / Nevent;
	//normalisation of baseline and area

	std::vector<Particle> vProductDs, vProductTau;
	std::vector<Particle>::iterator iP;

	int DecayCount = 0, InNDCount = 0, ID;
	for (ID = 0; ID < Nevent; ++ID)
	{
		double pt, xf;
		do
		{
			pt = mt->Uniform(0, ptmax);
			xf = mt->Uniform(-1.0, 1.0);
		}
		while (mt->Rndm() > Ds(xf, pt));

		double px, py;
		double pz = ptmax * xf / 2.0;
		mt->Circle(px, py, pt);

		TLorentzVector Ds_vec(px, py, pz, sqrt(pt*pt + pz*pz + pow(Const::MDs, 2)));
		Ds_vec.Boost(S.BoostVector());	//parent lab frame

		//Ds decay into electrons
		if (!nuEFile.empty())
		{
			std::vector<Particle> vProductDs = Nu0.ProductionPS(Ds_vec, "CharmE");
			if (vProductDs.size() && vProductDs[0].Theta() <= th0)
				hCharmE->Fill(vProductDs[0].Energy(), SF * 8.3e-5);
		}							//BR(Ds -> e nu)

		//Ds decay into muons
		if (!nuMFile.empty())
		{
			std::vector<Particle> vProductDs = Nu0.ProductionPS(Ds_vec, "CharmM");
			if (vProductDs.size() && vProductDs[0].Theta() <= th0)
				hCharmM->Fill(vProductDs[0].Energy(), SF * 5.5e-3);	//corrected BR
		}							//BR(Ds -> mu nu)

		//Ds decay into taus
		std::vector<Particle> vProductDs = Nu0.ProductionPS(Ds_vec, "CharmT");

		if (vProductDs.size())
		{
			++DecayCount;

			if (vProductDs[0].Theta() <= th0)	//neutrino
				hCharmT->Fill(vProductDs[0].Energy(), SF * 0.0548);
			else
				++InNDCount;
		}

		if (mass > 0)	//use light neutrino to have tau flux 
		{
			Neutrino lightNu(0.0,  Neutrino::Dirac | Neutrino::Left );
			vProductDs = lightNu.ProductionPS(Ds_vec, "CharmT");
		}


		if (vProductDs.size() < 2)
			continue;	//generation failed

		//tau decay from Ds
		TLorentzVector Tau_vec(vProductDs[1].FourVector());
		for (int i = 0; i < 4; ++i)
		{
			std::string channel;
			TH1D* hFill;
			double Br;
			switch (i)
			{
				case 0:
					channel = "TauET";
					hFill = hTauE;
					Br = SF * 0.0548 * 0.1785;	//tau->e (17.85 %)
					break;
				case 1:
					channel = "TauMT";
					hFill = hTauM;
					Br = SF * 0.0548 * 0.1736;	//tau->mu (17.36 %)
					break;
				case 2:
					channel = "TauPI";
					hFill = hPion;
					Br = SF * 0.0548 * 0.1082;	//tau->pi (10.82 %)
					break;
				case 3:
					channel = "Tau2PI";
					hFill = h2Pion;
					Br = SF * 0.0548 * 0.2551;	//tau->2pi (25.62 %)
					break;			//Phase space only!!
				default:
					break;
			}

			std::vector<Particle> vProductTau = NuB.ProductionPS(Tau_vec, channel);
			if (vProductTau.size())
				if (vProductTau.at(0).Theta() <= th0)	//neutrino
					hFill->Fill(vProductTau.at(0).Energy(), Br);
		}

		if (ID % 10000 == 0)	//saving
		{
			FileOut0->Write("", TObject::kOverwrite);
			FileOutB->Write("", TObject::kOverwrite);
		}

	}

	hTotal0->Add(hCharmT);

	hTotalB->Add(hTauE);
	hTotalB->Add(hTauM);
	hTotalB->Add(hPion);
	hTotalB->Add(h2Pion);

	FileOut0->cd();
	FileOut0->Write();

	hTotal0->Write("htotal");
	//hCharm->Write();

	FileOutB->cd();
	FileOutB->Write();

	hTotalB->Write("htotal");

	std::cout << "Ds meson decays are " << 100.0 * DecayCount / double(Nevent) << " %\n";
	//std::cout << "Products in ND are " << 100.0 * (1.0 - InNDCount / double(Nevent)) << " %\n";
	std::cout << "Neutrinos in ND are " << 100.0 * (hCharmT->GetEntries() / double(Nevent)) << " %\n";
	std::cout << "Antineuts in ND are " << 100.0 * (hPion->GetEntries() / double(DecayCount)) << " %\n";
	std::cout  << "Neutrinos simulated " << hTotal0->GetEntries();
	std::cout << " (" << hTotal0->GetEntries()*100.0/double(Nevent) << " %)";
	std::cout << ", saved in " << FileOut0->GetName() << std::endl;
	std::cout  << "Antineutrinos simulated " << hTotalB->GetEntries();
	std::cout << " (" << hTotalB->GetEntries()*100.0/double(Nevent) << " %)";
	std::cout << ", saved in " << FileOutB->GetName() << std::endl;

	FileOut0->Close();
	FileOutB->Close();

	if (!nuEFile.empty())	//creating new files so that it doesn't screw up existing files
	{
		TFile fileInE(nuEFile.c_str(), "READ");

		if (nuEFile.find(".root") != std::string::npos)
			nuEFile.insert(nuEFile.find(".root"), "_new");
		else
			nuEFile += "_new.root";
		TFile fileOutE(nuEFile.c_str(), "RECREATE");

		fileInE.cd();
		TIter next(fileInE.GetListOfKeys());
		TKey *kkk;
		hTotal0 = 0;
		while (kkk = static_cast<TKey*> (next()))
		{
			if (kkk->GetName() == "htotal" || kkk->GetName() == "hcharm")
				continue;
			TH1D *hflux = static_cast<TH1D*> (kkk->ReadObj());
			if (!hTotal0)
				hTotal0 = hflux;
			else
				hTotal0->Add(hflux);
		}

		fileOutE.cd();
		if (hTotal0)
		{
			hTotal0->Add(hCharmE);
			hTotal0->Write("", TObject::kOverwrite);
		}
		hCharmE->Write("hcharm", TObject::kOverwrite);

		fileInE.Close();
		fileOutE.Close();
	}

	if (!nuMFile.empty())	//creating new files so that it doesn't screw up existing files
	{
		TFile fileInM(nuMFile.c_str(), "READ");

		if (nuMFile.find(".root") != std::string::npos)
			nuMFile.insert(nuMFile.find(".root"), "_new");
		else
			nuMFile += "_new.root";
		TFile fileOutM(nuMFile.c_str(), "RECREATE");

		fileInM.cd();
		TIter next(fileInM.GetListOfKeys());
		TKey *kkk;
		hTotal0 = 0;
		while (kkk = static_cast<TKey*> (next()))
		{
			if (kkk->GetName() == "htotal" || kkk->GetName() == "hcharm")
				continue;
			TH1D *hflux = static_cast<TH1D*> (kkk->ReadObj());
			if (!hTotal0)
				hTotal0 = hflux;
			else
				hTotal0->Add(hflux);
		}

		fileOutM.cd();
		if (hTotal0)
		{
			hTotal0->Add(hCharmM);
			hTotal0->Write("", TObject::kOverwrite);
		}
		hCharmM->Write("hcharm", TObject::kOverwrite);

		fileInM.Close();
		fileOutM.Close();
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
