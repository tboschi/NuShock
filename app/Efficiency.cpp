#include <cstdlib>
#include <iostream>
#include <unordered_map>
#include <string>

#include "tools/CardDealer.h"

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"

#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"
#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"

int main( int argc, char** argv )
{
	TMVA::Tools::Instance();

	// Register the training and test trees

	CardDealer cd(argv[1]);
	std::unordered_map<std::string, std::vector<std::string> > alls;
	if (!cd.Get("eff_", alls))
		throw std::invalid_argument("Card file is useless\n");

	std::string outname;
	if (!cd.Get("output", outname))
		outname = "tmva.root";
	if (outname.find(".root") == std::string::npos)
		outname += ".root";

	double eff;
	if (!cd.Get("efficiency", eff))
		eff = 0.7;

	std::vector<double> bins;
	if (!cd.Get("binning", bins)) {
		bins.reserve(126);
		std::generate_n(std::back_inserter(bins), 126,
				[n=0]() mutable { return 0.2 * n++; } );
	}

	for (const auto & chan : alls) {
		std::cout << "Training selection for " << chan.first << "\n";
		TFile* inSignal = TFile::Open(chan.second[0].c_str());
		TFile* inBackgr = TFile::Open(chan.second[1].c_str());

		TTree *signaltree = static_cast<TTree*>(inSignal->Get("hnl"));
		TTree *background = static_cast<TTree*>(inBackgr->Get("hnl"));

		std::string output = outname;
		output.insert(output.find(".root"), "_" + chan.first);

		// Create a ROOT output file where TMVA will store ntuples, histograms, etc.
		TFile* outTMVA = TFile::Open(output.c_str(), "RECREATE");

		std::string tmvaname = "cuts_ " + chan.first;
		TMVA::Factory *factory = new TMVA::Factory(tmvaname.c_str(), outTMVA,
				"Silent:DrawProgressBar:AnalysisType=Classification");
				//:Transformations=I;D;P;G,D ??

		TMVA::DataLoader *dl = new TMVA::DataLoader("GA");

		// Define the input variables that shall be used for the MVA training
		// note that you may also use variable expressions, such as: "3*var1/var2*abs(var3)"
		// [all types of expressions that can also be parsed by TTree::Draw( "expression" )]
		dl->AddVariable("E_A", 'F');
		dl->AddVariable("E_B", 'F');
		dl->AddVariable("E_0", 'F');
		dl->AddVariable("M_0", 'F');
		dl->AddVariable("T_A", 'F');
		dl->AddVariable("T_B", 'F');
		dl->AddVariable("T_0", 'F');
		dl->AddVariable("TheA", 'F');
		dl->AddVariable("TheB", 'F');
		dl->AddVariable("The0", 'F');
		dl->AddVariable("Angle", 'F');

		// true energy is spectator
		// from which efficiency will be calculated
		dl->AddSpectator("True", 'F');

		// You can add an arbitrary number of signal or background trees
		dl->AddSignalTree    (signaltree);
		dl->AddBackgroundTree(background);

		dl->SetWeightExpression("W");
		//dataloader->SetBackgroundWeightExpression( "W" );

		// Tell the dataloader how to use the training and testing events
		//
		// If no numbers of events are given, half of the events in the tree are used
		// for training, and the other half for testing:
		//
		//    dataloader->PrepareTrainingAndTestTree( mycut, "SplitMode=random:!V" );
		//dataloader->PrepareTrainingAndTestTree( mycuts, mycutb,
				//"nTrain_Signal=1000:nTrain_Background=1000:SplitMode=Random:NormMode=NumEvents:!V" );

		factory->BookMethod(dl, TMVA::Types::kCuts, "CutsGA"); //
				//"!H:!V:FitMethod=MC:EffSel:SampleSize=200000:VarProp=FSmart" );

		// Now you can tell the factory to train, test, and evaluate the MVAs
		//
		// Train MVAs using the set of training events
		factory->TrainAllMethods();

		// Evaluate all MVAs using the set of test events
		factory->TestAllMethods();

		// Evaluate and compare performance of all configured MVAs
		//factory->EvaluateAllMethods();

		// Save the output
		outTMVA->Close();

		std::cout << "==> Wrote root file: " << outTMVA->GetName() << std::endl;
		std::cout << "==> TMVAClassification is done!" << std::endl;

		delete factory;
		delete dl;



		// create reader to read weights
		TMVA::Reader *reader = new TMVA::Reader( "!Color:!Silent" );

		float E_A, E_B, E_0, M_0, T_A, T_B, T_0, TheA, TheB, The0, Angle;
		float True, W;
		reader->AddVariable("E_A", &E_A);
		reader->AddVariable("E_B", &E_B);
		reader->AddVariable("E_0", &E_0);
		reader->AddVariable("M_0", &M_0);
		reader->AddVariable("T_A", &T_A);
		reader->AddVariable("T_B", &T_B);
		reader->AddVariable("T_0", &T_0);
		reader->AddVariable("TheA", &TheA);
		reader->AddVariable("TheB", &TheB);
		reader->AddVariable("The0", &The0);
		reader->AddVariable("Angle", &Angle);

		reader->AddSpectator("True", &True);
		reader->AddSpectator("W", &W);

		// get cuts
		std::string bookname = "dataset/weights/" + tmvaname + ".weights.xml";
		reader->BookMVA("CutsGA method", bookname.c_str());

		// tree for computing efficiency
		signaltree->SetBranchAddress("E_A", &E_A);
		signaltree->SetBranchAddress("E_B", &E_B);
		signaltree->SetBranchAddress("E_0", &E_0);
		signaltree->SetBranchAddress("M_0", &M_0);
		signaltree->SetBranchAddress("T_A", &T_A);
		signaltree->SetBranchAddress("T_B", &T_B);
		signaltree->SetBranchAddress("T_0", &T_0);
		signaltree->SetBranchAddress("TheA", &TheA);
		signaltree->SetBranchAddress("TheB", &TheB);
		signaltree->SetBranchAddress("The0", &The0);
		signaltree->SetBranchAddress("Angle", &Angle);
		// --------------------------------------------------------------
		
		TH1D * hAll = new TH1D("hall", "hall", bins.size()-1, &bins[0]);
		TH1D * hFit = new TH1D("hfit", "hfit", bins.size()-1, &bins[0]);

		for (long int ievt = 0; ievt < signaltree->GetEntries(); ++ievt) {
			signaltree->GetEntry(ievt);

			if (reader->EvaluateMVA("CutsGA method", eff))
				hFit->Fill(True);
			hAll->Fill(True);
		}

		for (int bin = 1; bin <= hAll->GetNbinsX(); ++bin)
			hFit->SetBinContent(bin,
				hFit->GetBinContent(bin) / hAll->GetBinContent(bin));

	}

	return 0;
}
