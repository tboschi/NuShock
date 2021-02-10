#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <cstring>
#include <getopt.h>

#include "TTree.h"
#include "TList.h"
#include "TFile.h"
#include "TGraph.h"

#include "tools/CardDealer.h"

#include "physics/Flavours.h"
#include "physics/Neutrino.h"
#include "physics/Decays.h"
#include "physics/ProductionRate.h"

#include "detector/Driver.h"
#include "detector/Tracker.h"
#include "montecarlo/GENIEback.h"

int main(int argc, char** argv)
{

	const struct option longopts[] = 
	{
		{"channel", 	required_argument,	0, 'c'},
		{"chargeID", 	required_argument,	0, 'I'},
		{"verbose", 	required_argument,	0, 'v'},
		{"help", 	no_argument,	 	0, 'h'},
		{0,	0, 	0,	0},
	};

	int index; 
	int iarg = 0;
	opterr = 1;
	bool chargeID = false, verbose = false;
	
	std::string channel;
	while((iarg = getopt_long(argc,argv, "c:Iv", longopts, &index)) != -1)
	{
		switch(iarg)
		{
			case 'c':
				channel.assign(optarg);
				break;
			case 'I':
				chargeID = true;
				break;
			case 'v':
				verbose = true;
				break;
			default:
				break;
		}
	}

	// detection channel
	std::vector<Decay::Channel> chans;
	std::stringstream sschan(channel);
	std::string llchan;
	while (std::getline(sschan, llchan, ','))
		chans.push_back(Decay::fromString(llchan));

	// for flux building
	std::string flux(argv[optind]);
	Driver drive(flux);

	// for detector configuration
	std::string det(argv[optind+1]);
	Tracker trk(det);

	// for background info
	CardDealer cd(argv[optind+2]);

	// contains xsec and gst per neutrino flavour
	std::map<std::string, std::vector<std::string> > file_backs;
	if (!cd.Get("back_", file_backs))
		throw std::invalid_argument("Couldn't find xsec names and backgrounds for fluxes\n");

	std::string file_xsec;
	if (!cd.Get("xsec_file", file_xsec))
		throw std::invalid_argument("Couldn't find xsec files\n");
	TFile inxsec(file_xsec.c_str());

	ProductionRate phnl;    // for flux making of a light neutrino

	std::string outname;
	if (!cd.Get("output", outname))
		throw std::logic_error("GenBack: no output \".root\" file specified in card");

	if (outname.find(".root") == std::string::npos)
		outname += ".root";
	if (chargeID)
		outname.insert(outname.find(".root"), "_I");

	std::map<Decay::Channel, std::string> files;

	for (const auto &fb : file_backs) {
		std::cout << "Doing flavour " << fb.first << "\n";
		Nu::Flavour flv = Nu::fromString(fb.first);
		// make flux
		//drive.MakeFlux(n, mix);
		// get full flux
		std::shared_ptr<TH1D> hist = drive.MakeComponent(phnl, flv, 0.);
		hist->Scale(std::pow(1./trk.Baseline(),2) * 1e20); //trk.POTs() * trk.Weight());

		std::string cc = fb.second[0] + "/tot_cc";
		std::string nc = fb.second[0] + "/tot_nc";
		std::shared_ptr<TGraph> xsecCC(static_cast<TGraph*>(inxsec.Get(cc.c_str())));
		std::shared_ptr<TGraph> xsecNC(static_cast<TGraph*>(inxsec.Get(nc.c_str())));

		double events = 0;
		for (int bin = 1; bin < hist->GetNbinsX(); ++bin) {
			double en = hist->GetBinCenter(bin);
			events += hist->GetBinContent(bin) * hist->GetBinWidth(bin)
			    * 1.e-38 * (xsecCC->Eval(en) + xsecNC->Eval(en));
		}
		// multiply per number of targets per ton
		// to get number of events per ton per 1e20 POTs
		events *= Const::Na * 1e6 / Material::A(Material::LAr);

		std::cout << "Event rate for " << fb.first << " are " << events << "\n";
		auto d0 = GENIEback::GenerateBackground(trk, chans, fb.second[1], events,
							chargeID, verbose);
		std::cout << "Background events:\n";
		size_t sum = 0;
		for (const auto & dc : d0) {
			sum += dc.second->GetEntries();

			std::string name = outname;
			name.insert(name.find(".root"), "_" + Decay::toString(dc.first));
			name.insert(name.find(".root"), "_" + fb.first);
			files[dc.first] += " " + name;

			std::cout << "\t" << Decay::toString(dc.first) << "   "
				  << dc.second->GetEntries() << " ("
				  << dc.second->GetEntries() / 1.e4 << " %)"
				  << "\t-> " << name << "\n";

			TFile out(name.c_str(), "RECREATE");
			out.cd();
			dc.second->Write();
		}
		std::cout << "     Total\t" << sum << "\n\n";
	}

	for (const auto & ff : files) { 
		std::cout << "solving for " << Decay::toString(ff.first) << "\n";
		std::string name = outname;
		name.insert(name.find(".root"), "_" + Decay::toString(ff.first));
		std::string cmd = "hadd -f " + name + ff.second;
		std::cout << cmd << "\n";
		system(cmd.c_str());
		cmd = "rm " + ff.second;
		std::cout << cmd << "\n";
		system(cmd.c_str());
	}

	return 0;
}
