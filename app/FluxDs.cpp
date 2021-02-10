#include <iostream>
#include <fstream>
#include <cstring>
#include <vector>
#include <random>
#include <memory>
#include <iomanip>
#include <getopt.h>

#include "tools/CardDealer.h"

#include "physics/Const.h"
#include "physics/ProductionSpace.h"
#include "physics/Particle.h"
#include "physics/Neutrino.h"
#include "physics/OpenQQ.h"

#include "detector/Driver.h"
#include "detector/Detector.h"

#include "TKey.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"

//from 1708.08700, 250GeV proton beam E796
const double b = 1.08;
const double n = 6.1;

double Dsdx(double xf, double pt)
{
	return std::exp(n * std::log(1 - std::abs(xf)) - b * std::pow(pt, 2));
}

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

	//std::string sProb, sTarg("H"), OutName;
	std::string det_card, nuEfile, nuMfile;
	std::ofstream outfile;
	double mass = -1.;

	while((iarg = getopt_long(argc, argv, "m:", longopts, &index)) != -1)
	{
		switch(iarg)
		{
			case 'm':
				mass = std::strtod(optarg, NULL);
				break;
			default:
				break;
		}
	}

	CardDealer cd(argv[optind]);

	// geometrical configuration
	std::string config;
	if (!cd.Get("detector_configuration", config))
		throw std::logic_error("There is no detector configuration\n");
	Detector box(config);


	// physics configuration
	size_t nEvents;
	if (!cd.Get("number_events", nEvents))
		nEvents = 1e5;

	double Eb = box.BeamEnergy();
	//if (!cd.Get("beam_energy", Eb))
		//Eb = 80;
	if (Eb <= Const::MProton)
		throw std::invalid_argument("Beam energy is too low!\n");

	double cme = std::sqrt(2 * Const::MProton * (Const::MProton + Eb));
	TLorentzVector S_vec(0, 0, std::sqrt(std::pow(Eb,2) - std::pow(Const::MProton,2)),
				   Eb + Const::MProton);

	double ccxsec;
	if (!cd.Get("cc_xsec", ccxsec)) {
		if (!cd.Get("opencc_configuration", config)) {
			ccxsec = 12;
			std::cerr << "No open charm configuration file, using CC xsec = " << ccxsec
				  << " at beam energy of " << Eb << " GeV\n";
		}
		else {	// compute on the fly
			std::unique_ptr<OpenQQ> xsec(new OpenQQ(config));
			xsec->SetCMEnergy(cme);
			double err, chi2;
			ccxsec = xsec->Integrate(err, chi2);
		}
	}
	std::cout << "FluxDs: center of mass energy " << cme << "\n";
	std::cout << "FluxDs: using open cc xsec " << ccxsec << "\n";

	// neutrino mass for tau flux
	if (mass < 0 && !cd.Get("neutrino_mass", mass))
		mass = 0.;

	//std::unique_ptr<TFile> fileOut0(new TFile(out0.c_str(), "RECREATE"));
	//std::unique_ptr<TFile> fileOutB(new TFile(outB.c_str(), "RECREATE"));

	std::vector<double> bins;
	if (!cd.Get("binning", bins)) {
		bins.reserve(126);
		std::generate_n(std::back_inserter(bins), 126,
				[n=0]() mutable { return 0.2 * n++; } );
	}

	//neutrino
	std::shared_ptr<TH1D> hCharmE(new TH1D("hcharme", "hcharm", bins.size()-1, &bins[0]));

	std::shared_ptr<TH1D> hCharmM(new TH1D("hcharmm", "hcharm", bins.size()-1, &bins[0]));
	std::shared_ptr<TH1D> hCharmT(new TH1D("hcharm",  "hcharm", bins.size()-1, &bins[0]));

	//antineutrino
	std::shared_ptr<TH1D> hTauE  (new TH1D("htaue",  "htaue",   bins.size()-1, &bins[0]));
	std::shared_ptr<TH1D> hTauM  (new TH1D("htaum",  "htaum",   bins.size()-1, &bins[0]));
	std::shared_ptr<TH1D> hPion  (new TH1D("hpion",  "hpion",   bins.size()-1, &bins[0]));
	std::shared_ptr<TH1D> h2Pion (new TH1D("h2pion", "h2pion",  bins.size()-1, &bins[0]));

	std::shared_ptr<TH1D> hTauEE (new TH1D("htauee", "htaue",   bins.size()-1, &bins[0]));
	std::shared_ptr<TH1D> hTauMM (new TH1D("htaumm", "htaum",   bins.size()-1, &bins[0]));

	// Normalisation from open charm calculation and proper
	// histogram units -> nu / cm2 / GeV @ 1 m from source
	//	    cc xsec /  pA xsec * fragmentation / nEvents
	double SF = ccxsec * 1e-3 / 331.4 * 0.077 / nEvents
	//	    baseline / surface in cm2 / GeV
		    / box.Scaling();
	double ptmax = cme / 2.;	// half of sqrt(s)
	
	// simulation config
	std::uniform_real_distribution<> rndm;

	//neutrinos
	//Nu0 and NuB are to directly produce HNL flux from tau mixing
	//if mass == 0 then it is light neutrino flux
	//no difference between Majorana or Dirac

	std::cout << "FluxDs: creating neutrinos\n";

	std::cout << "FluxDs: creating phasespace terms\n";
	ProductionSpace ps0(Neutrino(mass, Neutrino::dirac | Neutrino::left | Neutrino::particle));
	ProductionSpace psB(Neutrino(mass, Neutrino::dirac | Neutrino::right | Neutrino::antiparticle));
	ProductionSpace light(Neutrino(0., Neutrino::dirac | Neutrino::left));

	std::vector<Particle> vProductDs, vProductTau;
	std::vector<Particle>::iterator iP;

	std::cout << "FluxDs: running simualation for " << nEvents << " events\n";
	for (size_t ID = 0; ID < nEvents; ++ID)
	{
		double pt = rndm(RNG::_mt) * ptmax;
		double xf = rndm(RNG::_mt) * 2. - 1.;
		while (rndm(RNG::_mt) > Dsdx(xf, pt)) {
			pt = rndm(RNG::_mt) * ptmax;
			xf = rndm(RNG::_mt) * 2. - 1.;
		}

		double theta = rndm(RNG::_mt) * 2. * Const::pi;
		double px = std::cos(theta) * pt;
		double py = std::sin(theta) * pt;
		double pz = cme * xf / 2.0;

		Particle Ds(431, std::sqrt(pt*pt + pz*pz + std::pow(Const::MDs, 2)), px, py, pz);
		Ds.Boost(S_vec.BoostVector());	//parent lab frame

		//Ds decay into electrons
		if (cd.Get("nuE0_file"))
		{
			// pair of particle vector and event weight
			const auto prodDs = ps0.Generate(Production::Channel::CharmE, Ds);
			if (prodDs.first.size() && box.AngularAccept(prodDs.first[0]))
				hCharmE->Fill(prodDs.first[0].E(), prodDs.second * SF * 8.3e-5);
		}							//BR(Ds -> e nu)

		//Ds decay into muons
		if (cd.Get("nuM0_file"))
		{
			const auto prodDs = ps0.Generate(Production::Channel::CharmM, Ds);
			if (prodDs.first.size() && box.AngularAccept(prodDs.first[0]))
				hCharmM->Fill(prodDs.first[0].E(), prodDs.second * SF * 5.5e-3);
		}							//BR(Ds -> mu nu)

		//Ds decay into taus
		auto prodDs = ps0.Generate(Production::Channel::CharmT, Ds);
		//std::cout << "ps0 generated: ";
		//for (const auto &p : prodDs.first)
		//	std::cout << p << "\t";
		//std::cout << "\nwith weight " << prodDs.second << "\n";

		// Ds decay to tau successful
		if (prodDs.first.size() && box.AngularAccept(prodDs.first[0]))
			hCharmT->Fill(prodDs.first[0].E(), prodDs.second * SF * 0.0548);

		if (mass > 0)	//use light neutrino to make tau flux 
			prodDs = light.Generate(Production::Channel::CharmT, Ds);

		if (prodDs.first.size() < 2)
			continue;

		//tau decay from Ds
		Particle tau(prodDs.first[1]);

		//tau decay into electrons via electron mixing
		if (cd.Get("nuEB_file"))
		{	//tau->e vie e mix (17.85 %)
			const auto prodTau = light.Generate(Production::Channel::TauEE, tau);
			if (prodTau.first.size() && box.AngularAccept(prodTau.first[0]))
				hTauEE->Fill(prodTau.first[0].E(), prodTau.second * SF
								  * 0.0548 * 0.1785);
		}

		//Ds decay into muons
		if (cd.Get("nuMB_file"))
		{	//tau->mu vie mu mix (17.36 %)
			const auto prodTau = psB.Generate(Production::Channel::TauMM, tau);
			if (prodTau.first.size() && box.AngularAccept(prodTau.first[0]))
				hTauMM->Fill(prodTau.first[0].E(), prodTau.second * SF
								  * 0.0548 * 0.1736);
		}

		for (size_t t = 0; t < 4; ++t)
		{
			Production::Channel chan;
			std::shared_ptr<TH1D> hist;
			double br;
			switch (t)
			{
				case 0:
					chan = Production::Channel::TauET;
					hist = hTauE;
					br = 0.1785;	//tau->e (17.85 %)
					break;
				case 1:
					chan = Production::Channel::TauMT;
					hist = hTauM;
					br = 0.1736;	//tau->mu (17.36 %)
					break;
				case 2:
					chan = Production::Channel::TauPI;
					hist = hPion;
					br = 0.1082;	//tau->pi (10.82 %)
					break;
				case 3:
					chan = Production::Channel::Tau2PI;
					hist = h2Pion;
					br = 0.2551;	//tau->2pi (25.62 %)
					break;		//Phase space only!!
				default:
					break;
			}

			br *= prodDs.second * SF * 0.0548;

			const auto prodTau = psB.Generate(chan, tau);
			if (prodTau.first.size() && box.AngularAccept(prodTau.first[0]))
				hist->Fill(prodTau.first[0].E(), br);
		}
	}


	// output files
	std::string outname;
	if (!cd.Get("output", outname))
		outname = "nu";

	std::string out0 = outname + "0.root";
	std::string outB = outname + "B.root";
	if (mass > 0) {
		std::stringstream ssm;
		ssm << "_" << std::setw(4) << std::fixed << std::setfill('0')
		    << std::setprecision(0) << mass * 1000;
		out0.insert(out0.find(".root"), ssm.str());
		outB.insert(outB.find(".root"), ssm.str());
	}

	TFile fileOut0(out0.c_str(), "RECREATE");
	hCharmT->Scale(1, "WIDTH");
	hCharmT->Write("", TObject::kOverwrite);
	fileOut0.Close();

	TFile fileOutB(outB.c_str(), "RECREATE");
	hTauE->Scale(1, "WIDTH");
	hTauE->Write("", TObject::kOverwrite);
	hTauM->Scale(1, "WIDTH");
	hTauM->Write("", TObject::kOverwrite);
	hPion->Scale(1, "WIDTH");
	hPion->Write("", TObject::kOverwrite);
	h2Pion->Scale(1, "WIDTH");
	h2Pion->Write("", TObject::kOverwrite);
	fileOutB.Close();


	std::cout << "Neutrinos in ND are\t" << 100.0 * (hCharmT->GetEntries() / double(nEvents)) << " %\n";
	std::cout << "Antineuts in ND are\t" << 100.0 * (hPion->GetEntries() / double(nEvents)) << " %\n";


	// extra files
	for (size_t i = 0; i < 4; ++i) {
		std::string extra;
		std::shared_ptr<TH1D> hist;
		switch (i) {
			case 0:
				extra = "nuE0_file";
				hist = hCharmE;
				break;
			case 1:
				extra = "nuM0_file";
				hist = hCharmM;
				break;
			case 2:
				extra = "nuEB_file";
				hist = hTauEE;
				break;
			case 3:
				extra = "nuMB_file";
				hist = hTauMM;
				break;
		}
		std::string extra_file;
		if (!cd.Get(extra, extra_file))
			continue;

		TFile fileIn(extra_file.c_str(), "READ");

		if (extra_file.find(".root") != std::string::npos)
			extra_file.insert(extra_file.find(".root"), "_new");
		else
			extra_file += "_new.root";

		TFile fileOut(extra_file.c_str(), "RECREATE");

		TIter next(fileIn.GetListOfKeys());
		TKey *kkk;
		fileOut.cd();
		while ((kkk = static_cast<TKey*> (next())))
		{
			// skip htaus or hcharm as we are overwriting them
			if ((std::strcmp(kkk->GetName(), "hcharm") == 0)
			 || (std::strcmp(kkk->GetName(), "htaue") == 0)
			 || (std::strcmp(kkk->GetName(), "htaum") == 0))
				continue;

			// copy old histogram
			std::shared_ptr<TH1D> hflux(static_cast<TH1D*>(kkk->ReadObj()));
			//hflux->Scale(1, "WIDTH");  DO NOT SCALE!!!
			hflux->Write(kkk->GetName(), TObject::kOverwrite);
		}

		// write last ones
		hist->Scale(1, "WIDTH");
		hist->Write(hist->GetTitle(), TObject::kOverwrite);

		fileIn.Close();
		fileOut.Close();
	}

	return 0;
}
