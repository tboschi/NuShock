#include <iostream>
#include <string>
#include <vector>
#include <map>

#include "physics/Flavours.h"
#include "detector/Flux.h"
#include "detector/Detector.h"

#include "flux/dk2nuTree.h"

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1D.h"

class bsim;

int main(int argc, char **argv) {

	std::cout << "Opening files " << argv[1] << "\n";

	TFile *ff = new TFile(argv[1]);
	dk2nuTree *dk = new dk2nuTree(static_cast<TTree*>(ff->Get("dk2nuTree")));
	//TTree *dk = static_cast<TTree*>(ff->Get("dk2nuTree"));
	//TChain *dk = new TChain("dk2nuTree");
	//size_t pots = 100000 * dk->Add(argv[1]);	// add files to chain
	double pots = 25 * 100000;	// add files to chain

	std::cout << "scanning over " << dk->fChain->GetEntries() << "\n";

	/*
	int nuray_ = 0, ancestor_ = 0, traj_ = 0;

	int decay_ntype = 0;
	int decay_ptype = 0;
	double decay_nimpwt = 0.;
	double nuray_E[3] = {0., 0., 0.};
	double nuray_wgt[3] = {0., 0., 0.};

	//dk2nuTree *dk = new dk2nuTree(tt);

	dk->SetBranchAddress("nuray", &nuray_);
	dk->SetBranchAddress("ancestor", &ancestor_);
	dk->SetBranchAddress("traj", &traj_);
	dk->SetBranchAddress("decay.ntype", &decay_ntype);
	dk->SetBranchAddress("decay.ptype", &decay_ptype);
	dk->SetBranchAddress("decay.nimpwt", &decay_nimpwt);
	dk->SetBranchAddress("nuray.E", nuray_E);
	dk->SetBranchAddress("nuray.wgt", nuray_wgt);

	dk->SetBranchStatus("*", 0);
	dk->SetBranchStatus("nuray", 1);
	dk->SetBranchStatus("decay.ntype", 1);
	dk->SetBranchStatus("decay.ptype", 1);
	dk->SetBranchStatus("decay.nimpwt", 1);
	dk->SetBranchStatus("nuray.E", 1);
	dk->SetBranchStatus("nuray.wgt", 1);
	*/

	//return 1;

	// total number of POTs
	//size_t pots = 100000;
	//size_t pots = 100000 * dk->Add(argv[1]);	// add files to chain
	//std::shared_ptr<dk2nuTree> dk2nu(dk);
	//
	double binning[] = {0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1, 1.125, 1.25, 1.375, 1.5, 1.625, 1.75, 1.875, 2, 2.125, 2.25, 2.375, 2.5, 2.625, 2.75, 2.875, 3, 3.125, 3.25, 3.375, 3.5, 3.625, 3.75, 3.875, 4, 4.125, 4.25, 4.375, 4.5, 4.625, 4.75, 4.875, 5, 5.125, 5.25, 5.375, 5.5, 5.625, 5.75, 5.875, 6, 6.125, 6.25, 6.375, 6.5, 6.625, 6.75, 6.875, 7, 7.125, 7.25, 7.375, 7.5, 7.625, 7.75, 7.875, 8, 8.5, 9, 9.5, 10, 10.5, 11, 11.5, 12, 12.5, 13, 13.5, 14, 14.5, 15, 15.5, 16, 16.5, 17, 17.5, 18, 18.5, 19, 19.5, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 40};
	size_t bins = sizeof(binning)/sizeof(binning[0]) - 1;
	
	std::map<int, Flux::Parent> pdgs = { {211, Flux::Parent::Pion},
					     {321, Flux::Parent::Kaon},
					     {130, Flux::Parent::Kaon0},
					     {13,  Flux::Parent::Muon} };
	std::map<Nu::Flavour, Flux::Component> all;
	for (const auto flv : Nu::All()) {
		Flux::Component nu_comp;
		for (const auto &pp : pdgs) {
			std::string name = "h" + Flux::toString(pp.second) + "_" + Nu::toString(flv);
			nu_comp[pp.second] = std::shared_ptr<TH1D>(new TH1D(name.c_str(),
						//name.c_str(), 125, 0., 25.));
						name.c_str(), bins, binning));
		}
		all[flv] = std::move(nu_comp);
	}

	for (size_t i = 0; i < dk->fChain->GetEntries(); ++i) {
		dk->fChain->GetEntry(i);

		// get neutrino flavour
		Nu::Flavour flv = Nu::fromPDG(dk->decay_ntype);
		if (!all.count(flv))	// odd
			continue;
		if (!pdgs.count(std::abs(dk->decay_ptype)))	// odd, again
			continue;
		Flux::Parent pp = pdgs[std::abs(dk->decay_ptype)];
		if (!all[flv].count(pp))	// mega odd!!
			continue;
		/*
		std::cout << "event " << i << "\t" << dk->nuray_ << "\t"
			  << Nu::toString(flv) << " from " << Flux::toString(pp)
			  << ", energy " << dk->nuray_E[1] << " wgt "
			  << dk->decay_nimpwt * dk->nuray_wgt[1] << "\n";
		*/

		all[flv][pp]->Fill(dk->nuray_E[1], dk->decay_nimpwt * dk->nuray_wgt[1]);
	}

	// base name
	std::string out(argv[2]);
	for (const auto & af : all) {
		std::string file = out + "_" + Nu::toString(af.first) + ".root";
		TFile outf(file.c_str(), "RECREATE");
		for (const auto & comp : af.second) {
			double ff = 1. / pots * (574. * 574.) / (Const::pi * 1.e4);
			comp.second->Scale(ff, "WIDTH");
			std::string name = "h" + Flux::toString(comp.first);
			comp.second->Write(name.c_str());
		}
		outf.Close();
	}

	//delete dk;

	return 0;
}
