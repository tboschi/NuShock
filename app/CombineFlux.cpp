#include <iostream>
#include <fstream>
#include <vector>
#include <cstring>
#include <memory>

#include "TH1D.h"
#include "TFile.h"

#include "detector/Flux.h"
#include "tools/CardDealer.h"

int main(int argc, char** argv)
{
	// for background info
	CardDealer cd(argv[1]);

	std::map<std::string, std::vector<std::string> > fluxes;
	if (!cd.Get("flux_", fluxes))
		throw std::invalid_argument("Couldn't find xsec and backgrounds for fluxes\n");

	for (const auto &fb : fluxes) {
		TFile fhc(fb.second[0].c_str());
		TFile rhc(fb.second[1].c_str());
		TFile out(fb.second[2].c_str(), "RECREATE");
		out.cd();

		std::shared_ptr<TH1D> hall = nullptr;
		for (const auto &parent : Flux::Parents()) {
			std::string name = "h" + Flux::toString(parent);

			std::shared_ptr<TH1D> hfhc(static_cast<TH1D*>(fhc.Get(name.c_str())));
			std::shared_ptr<TH1D> hrhc(static_cast<TH1D*>(rhc.Get(name.c_str())));
			
			if (!hfhc && !hrhc)
				continue;

			std::shared_ptr<TH1D> hsum = hfhc ? hfhc : hrhc;
			if (hfhc && hrhc)
				hsum->Add(hrhc.get());
			hsum->Scale(0.5);

			if (!hall)
				hall = std::shared_ptr<TH1D>(static_cast<TH1D*>(hsum->Clone()));
			else
				hall->Add(hsum.get());
			hsum->Write(name.c_str());
		}
		hall->Write("htotal");

		//out.Close();
		//fhc.Close();
		//rhc.Close();
	}

	return 0;
}
