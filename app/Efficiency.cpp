#include <cstdlib>
#include <iostream>
#include <unordered_map>
#include <string>
#include <getopt.h>

#include "tools/CardDealer.h"
#include "montecarlo/hnl.h"

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TArrayD.h"

#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"
#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"


void trainTMVA(std::string tmvaout, std::string tmvaname,
		TTree* signaltree, TTree* background,
		bool chID = false) {
	/* this routines leaks memory!!
	 * Maybe TMVA is not supposed to be called over and over again...
	 * unfortunately the documentation is not so great
	 */

	std::cout << "==> Training TMVA with name " << tmvaname << " to file "
		<< tmvaout << "\n";

	// Create a ROOT output file where TMVA will store ntuples, histograms, etc.
	TFile* outTMVA = TFile::Open(tmvaout.c_str(), "RECREATE");
	TMVA::Factory* ft = new TMVA::Factory(tmvaname.c_str(), outTMVA,
		"Silent:!V:!DrawProgressBar:AnalysisType=Classification");
	//":Transformations=I;G;U;P;D;U,G;U,P");

	TMVA::DataLoader* dl = new TMVA::DataLoader("efftmva");

	// [all types of expressions that can also be parsed by TTree::Draw( "expression" )]
	//dl->AddVariable("E_A", 'F');
	//dl->AddVariable("E_B", 'F');
	dl->AddVariable("E_0", 'F');
	dl->AddVariable("M_0", 'F');
	dl->AddVariable("log(T_A)", 'F');
	dl->AddVariable("log(T_B)", 'F');
	dl->AddVariable("log(T_0)", 'F');
	dl->AddVariable("log(TheA)", 'F');
	dl->AddVariable("log(TheB)", 'F');
	dl->AddVariable("log(The0)", 'F');
	dl->AddVariable("Angle", 'F');
	//dl->AddVariable("CosAB := abs(cos(PhiA-PhiB))", 'F');
	//dl->AddVariable("T_AB0 := T_0 - abs(T_A-T_B)", 'F');

	// true energy is spectator
	// from which efficiency will be calculated

	// You can add an arbitrary number of signal or background trees
	dl->AddSignalTree    (signaltree);
	dl->AddBackgroundTree(background);

	//dl->SetWeightExpression("W");

	TCut cut_sig = "";
	TCut cut_bkg = "";
	if (chID) {
		if (signaltree->GetEntries("ChID") < 1 || background->GetEntries("ChID") < 1) {
			std::cout << "\tNo events satisfying charge-ID cuts\n";
			outTMVA->Close();
			delete ft;
			delete dl;

			return;
		}
		cut_sig = cut_bkg = "ChID";
	}
	std::cout << "Cutting signal and background with " << cut_sig << " and " << cut_bkg << "\n";


	//int nTrainSig = 0.5 * signaltree->GetEntries();
	//int nTrainBkg = 0.5 * background->GetEntries();

	// build the string options for DataLoader::PrepareTrainingAndTestTree
	//TString prepareOpts = TString::Format("nTrain_Signal=%d:nTrain_Background=%d:SplitMode=Random:", nTrainSig, nTrainBkg);
	// If no numbers of events are given, half of the events in the tree are used
	// for training, and the other half for testing:
	dl->PrepareTrainingAndTestTree(cut_sig, cut_bkg, "SplitMode=Random:NormMode=NumEvents" );

	ft->BookMethod(dl, TMVA::Types::kFisher, "Fisher"); //, "VarTransform=U");

	// Now you can tell the factory to train, test, and evaluate the MVAs
	// Train MVAs using the set of training events
	ft->TrainAllMethods();

	// Evaluate all MVAs using the set of test events
	ft->TestAllMethods();

	// Evaluate and compare performance of all configured MVAs
	ft->EvaluateAllMethods();

	std::cout << "==> Wrote root file: " << outTMVA->GetName() << std::endl;
	std::cout << "==> TMVAClassification is done!" << std::endl;

	// Save the output
	outTMVA->Close();

	delete ft;	// destructor calls DeleteAllMethods()
	delete dl;
}

void evalTMVA(std::string tmvaout, std::string tmvaname, std::string outname,
		const std::vector<double> &e_bins,
		TTree* signaltree, TTree* background,
		bool chID = false,
		double min_eff = -1., double max_back = -1.) {

	std::cout << "==> Evaluating TMVA with name " << tmvaname << " from file "
		<< tmvaout << "\n";
	TFile *info = TFile::Open(tmvaout.c_str());
	// obtain
	TH1D* sig_eff = dynamic_cast<TH1D*>(info->Get("efftmva/Method_Fisher/Fisher/MVA_Fisher_effS"));
	TH1D* bkg_eff = dynamic_cast<TH1D*>(info->Get("efftmva/Method_Fisher/Fisher/MVA_Fisher_effB"));

	// Create a ROOT output file where TMVA will store ntuples, histograms, etc.
	TFile* outEff = TFile::Open(outname.c_str(), "RECREATE");
	std::cout << "+++ saving efficiency information in " << outname << " +++\n";

	TH1D *hEff = new TH1D("heff", "energy efficiency", e_bins.size()-1, &e_bins[0]);
	TH1D* hCut = new TH1D("hcut", "cut", e_bins.size()-1, &e_bins[0]);

	TArrayD evts(3);	// store efficiency info here

	double eff = 0.0;
	TMVA::Reader* reader = NULL;
	float E_0, M_0, T_A, T_B, T_0, TheA, TheB, The0, Angle;

	if (sig_eff && bkg_eff) {
		double maxp = -1;
		for (int bin = 1; bin < sig_eff->GetNbinsX()+1; ++bin) {
			double S = sig_eff->GetBinContent(bin) * signaltree->GetEntries();
			double B = bkg_eff->GetBinContent(bin) * background->GetEntries();

			double pure = 0.;
			//std::cout << "with " << S << " and " << B << "\t";
			if (B > max_back && sig_eff->GetBinContent(bin) > min_eff) {
				//std::cout << "first\n";
				pure = (S+B == 0) ? 0 : S / (S + B);
			}
			else {
				//std::cout << "second\n";
				pure = (S+B == 0) ? 0 : S / std::sqrt(S + B);
			}

			if (pure > maxp) {
				maxp = pure;
				eff = sig_eff->GetBinCenter(bin);
			}
		}
		std::cout << "==> Optimal cut position is at " << eff << "\n";

		// create reader to read weights
		reader = new TMVA::Reader("!Color:Silent");

		//reader->AddVariable("E_A", &E_A);
		//reader->AddVariable("E_B", &E_B);
		reader->AddVariable("E_0", &E_0);
		reader->AddVariable("M_0", &M_0);
		reader->AddVariable("log(T_A)", &T_A);
		reader->AddVariable("log(T_B)", &T_B);
		reader->AddVariable("log(T_0)", &T_0);
		reader->AddVariable("log(TheA)", &TheA);
		reader->AddVariable("log(TheB)", &TheB);
		reader->AddVariable("log(The0)", &The0);
		reader->AddVariable("Angle", &Angle);
		//reader->AddVariable("CosAB := abs(cos(PhiA-PhiB))", &CosAB);
		//reader->AddVariable("T_AB0 := T_0 - abs(T_A-T_B)", &T_AB0);

		// get cuts
		std::string book = "efftmva/weights/" + tmvaname + "_Fisher.weights.xml";
		reader->BookMVA("Fisher method", book.c_str());
		std::cout << "==> Reading with method " << book << "\n";
	}
	else
		std::cerr << "### No training file, maximum efficiency ###\n";


	std::cout << "H0\n";
	hnl* sig = new hnl(signaltree);
	std::cout << "H1\n";

	//int m_bin = mass_ener->GetXaxis()->FindBin(mass);

	std::cout << "Loading signal entries\n";
	double tot = 0, cut = 0;
	for (size_t ievt = 0; ievt < sig->GetEntries(); ++ievt) {
		sig->GetEntry(ievt);
		tot += sig->W;
		hEff->Fill(sig->True, sig->W);

		if (chID && !sig->ChID)
			continue;

		//E_A   = sig->E_A;
		//E_B   = sig->E_B;
		E_0   = sig->E_0;
		M_0   = sig->M_0;
		T_A   = std::log(sig->T_A);
		T_B   = std::log(sig->T_B);
		T_0   = std::log(sig->T_0);
		TheA  = std::log(sig->TheA);
		TheB  = std::log(sig->TheB);
		The0  = std::log(sig->The0);
		Angle = sig->Angle;
		//CosAB = std::abs(std::cos(sig->PhiA - sig->PhiB));
		//T_AB0 = sig->T_0 - std::abs(sig->T_A - sig->T_B);
		//True =  sig->True;
		//W = sig->W;

		//if (method.find("Cuts") != std::string::npos)
		//	test = reader->EvaluateMVA(method.c_str(), eff);
		//else
		bool test = (reader ? reader->EvaluateMVA("Fisher method") > eff : true);

		if (test) {
			hCut->Fill(sig->True, sig->W);
			cut += sig->W;
		}
	}

	std::cout << " ****** SIGNAL RESULT IS = " << cut / tot << " ******\n";
	evts[0] = (tot > 0) ? cut / tot : 1.;
	//sig_evts->SetBinContent(m_bin, cut / tot);

	delete sig;

	std::cout << "Creating efficiency\n";

	for (int b = 1; b < hEff->GetNbinsX()+1; ++b) {
		double c = hEff->GetBinContent(b);
		if (c > 0.)
			hEff->SetBinContent(b, hCut->GetBinContent(b) / c);
		else
			hEff->SetBinContent(b, 1.);
	}
	hCut->Delete();


	hnl* bkg = new hnl(background);
	std::cout << "Loading background entries\n";
	tot = cut = 0;
	double cut_lnc = 0., cut_lnv = 0.;
	for (size_t ievt = 0; ievt < bkg->GetEntries(); ++ievt) {
		bkg->GetEntry(ievt);
		tot += bkg->W;

		if (chID && !bkg->ChID)
			continue;

		//E_A   = bkg->E_A;
		//E_B   = bkg->E_B;
		E_0   = bkg->E_0;
		M_0   = bkg->M_0;
		T_A   = bkg->T_A;
		T_B   = bkg->T_B;
		T_0   = bkg->T_0;
		TheA  = bkg->TheA;
		TheB  = bkg->TheB;
		The0  = bkg->The0;
		Angle = bkg->Angle;

		//if (method.find("Cuts") != std::string::npos)
		//	test = reader->EvaluateMVA(method.c_str(), eff);
		//else
		bool test = (reader ? reader->EvaluateMVA("Fisher method") > eff : true);

		if (test) {
			cut += bkg->W;
			if (bkg->ChA > 0)
				cut_lnc += bkg->W;
			else
				cut_lnv += bkg->W;
		}
	}
	std::cout << " ****** BACKGROUND RESULT IS = " << cut << "(" << cut_lnc << " + " << cut_lnv
		  << ") / " << tot << " ******\n";
	evts[1] = std::max(cut_lnc, tot / std::sqrt(12.) / bkg->GetEntries());
	evts[2] = std::max(cut_lnv, tot / std::sqrt(12.) / bkg->GetEntries());

	delete bkg;

	delete reader;

	outEff->cd();

	outEff->WriteObject(&evts, "events");
	hEff->Write();

	outEff->Close();
}

int main( int argc, char** argv )
{
	bool chID = false;

	int iarg;
	while((iarg = getopt(argc, argv, ":c")) != -1)
	{
		switch(iarg)
		{
			case 'c':
				chID =  true;
				break;
			default:
				break;
		}
	}

	std::cout << "charge ID is " << std::boolalpha << chID << "\n";
	CardDealer cd(argv[optind]);

	std::string sig_file(argv[optind+1]);
	sig_file.erase(0, sig_file.find_last_of('/')+1);
	sig_file.erase(0, sig_file.find_first_of('_')+1);
	sig_file.erase(sig_file.find(".root"));
	//chan_E_M_T_m_mass

	std::string channel = sig_file.substr(0, sig_file.find_first_of('_'));
	std::string m_str = sig_file.substr(sig_file.find_last_of('_')+1); 
	double mass = std::strtod(m_str.c_str(), NULL) / 1000.;

	std::cout << "Training channel " << channel << " with mass " << mass
		  << " (" << sig_file << ")\n";


	std::string bkg_file;
	if (!cd.Get("back_" + channel, bkg_file)) {
		std::cerr << "No background file for " << channel << "\n";
		return 1;
	}

	double min_eff, max_back;
	if (!cd.Get("efficiency", min_eff))
		min_eff = 1.;
	if (!cd.Get("background", max_back))
		max_back = 10.;


	std::string outname;
	if (!cd.Get("output", outname))
		outname = "efficiency.root";
	if (outname.find(".root") == std::string::npos)
		outname += ".root";
	outname.insert(outname.find(".root"), "_" + sig_file);

	if (chID)
		outname.insert(outname.find(".root"), "_ch");



	// energy direction bins
	std::vector<double> e_bins;
	if (!cd.Get("energy", e_bins)) {
		e_bins.reserve(126);
		std::generate_n(std::back_inserter(e_bins), 126,
				[n=0]() mutable { return 0.2 * n++; } );
	}

	bool train, evaluate;
	if (!cd.Get("train", train))
		train = false;
	if (!cd.Get("evaluate", evaluate))
		evaluate = false;


	//TMVA::Tools::Instance();

	TFile *inSignal = TFile::Open(argv[optind+1]);
	std::cout << "signal file " << inSignal->GetName() << "\n";
	TFile *inBackgr = TFile::Open(bkg_file.c_str());
	std::cout << "background file " << inBackgr->GetName() << "\n";

	TTree *signaltree = static_cast<TTree*>(inSignal->Get("hnl"));
	TTree *background = static_cast<TTree*>(inBackgr->Get("hnl"));

	// name of weight file, needed for reader as well
	std::string tmvaname = "cuts_" + channel + "_" + m_str + (chID ? "_ch" : "");

	// output name where training info is
	std::string tmvaout = "efftmva/TMVA.root";
	tmvaout.insert(tmvaout.find(".root"), "_" + sig_file);
	if (chID)
		tmvaout.insert(tmvaout.find(".root"), "_ch");

	if (train)	// can be easily skipper if already trained TMVA
		trainTMVA(tmvaout, tmvaname, signaltree, background, chID);
	if (evaluate)	// can be easily skipper if already trained TMVA
		evalTMVA(tmvaout, tmvaname, outname, e_bins,
				signaltree, background,
				chID, min_eff, max_back);

	inSignal->Close();
	inBackgr->Close();

	return 0;
}
