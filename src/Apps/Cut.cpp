// Do make to compile 
// Do AnalyseGENIE -f ghep_file.root -n events nbr to analyse

//____________________________________________________________________________
/*!

\program gtestEventLoop

\brief   Example event loop. Use that as a template for your analysis code.

\author  Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
		 STFC, Rutherford Appleton Laboratory

\created May 4, 2004

\cpright Copyright (c) 2003-2013, GENIE Neutrino MC Generator Collaboration
		 For the full text of the license visit http://copyright.genie-mc.org
		 or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#include <string>

#include "TSystem.h"
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TIterator.h"
#include "TH1.h"
#include "TGraph.h"

#include "EVGCore/EventRecord.h"
#include "GHEP/GHepParticle.h"
#include "Ntuple/NtpMCFormat.h"
#include "Ntuple/NtpMCTreeHeader.h"
#include "Ntuple/NtpMCEventRecord.h"
#include "Messenger/Messenger.h"
#include "PDG/PDGCodes.h"
#include "Utils/CmdLnArgParser.h"

using namespace genie;

void GetCommandLineArgs (int argc, char ** argv);
int gOptNEvt;
std::string gOptInpFilename;
std::string gOptOutFilename;
std::string gOptFluxFile;
std::string gOptSplineFile;
//___________________________________________________________________
int main(int argc, char ** argv)
{
	GetCommandLineArgs (argc, argv);
	// Open the ROOT file (XXX.ghep.root) (-i argument)
	// and get the TTree & its header
	TTree *tree = 0;
	genie::NtpMCTreeHeader *thdr = 0;
 	TFile file(gOptInpFilename.c_str(),"READ");
 	tree = dynamic_cast<TTree*>           ( file.Get("gtree")  );
 	thdr = dynamic_cast<NtpMCTreeHeader*> ( file.Get("header") );
 	if(!tree) return 1;
 	NtpMCEventRecord * mcrec = 0;
 	tree->SetBranchAddress("gmcrec", &mcrec);
 	// Get the nbr of evts to analyse (-n argument)
 	int nev = (gOptNEvt > 0) ? TMath::Min(gOptNEvt, (int)tree->GetEntries()) : (int) tree->GetEntries();
 	// Create a ROOT file (-o argument) with all the relevent histos
 	TFile outfile(gOptOutFilename.c_str(),"RECREATE");

//  	TH1D * h_pion_nrj    = new TH1D("pion_nrj",  "pion_nrj",    500,0,5);   // Pion energy
// 	TH1D * h_pion_mom    = new TH1D("pion_mom",  "pion_mom",    500,0,5);   // Pion mnomentum
//  	TH1D * h_pion_ang    = new TH1D("pion_ang",  "pion_ang",    180,0,180); // Pion angle
	unsigned int Index, nProton;
	double InvMass, NAngle, Angle, NEnergy, muTheta, piTheta;
	TTree * Trace = new TTree("trace", "Signal trace");
	Trace->Branch("Index", &Index, "Index/i");
	Trace->Branch("NEnergy", &NEnergy, "NEnergy/D");
	Trace->Branch("InvMass", &InvMass, "M2/D");
	Trace->Branch("NAngle", &NAngle, "NAngle/D");
	Trace->Branch("Angleup", &Angle, "Angle/D");
	Trace->Branch("muTheta", &muTheta, "muTheta/D");
	Trace->Branch("piTheta", &piTheta, "piTheta/D");
	Trace->Branch("Proton", &nProton, "sProton/i");

  	TH1D * hinProbe = new TH1D("Probe", "Probe energy", 1000,0,10);
  	TH1D * hNEnergy = new TH1D("NEnergy", "Sterile energy", 1000,0,10);
  	TH1D * hInvMass = new TH1D("InvMass", "Sterile invariant mass", 1000,0,10);
  	TH1D * hnnAngle = new TH1D("nnAngle", "Sterile entrance angle", 1000,0,TMath::Pi());
  	TH1D * hAngleup = new TH1D("Angleup", "Separation angle", 1000,0,TMath::Pi());
	
	// Get the splines used (converted in root using gspl2root)
	// This is needed to normalize your distribution
	TGraph* XS_graph = 0;
	//TFile splinefile(gOptSplineFile.c_str(),"READ");
	//gDirectory->GetObject("nu_mu_n/qel_cc_n",XS_graph);

	int cMu, cPi;
	int cMuF, cPiF;
	int total = 0;

	// Loop over all events
	for(int i = 0; i < nev; ++i)
       	{
		// get next tree entry
		tree->GetEntry(i);
		// get the GENIE event
		EventRecord & event = *(mcrec->event);
		// Print out all events infos
		//LOG("myAnalysis", pNOTICE) << event;

		// Initialize the particles info for each event
		TLorentzVector mu4, pi4; 	//four momentum of particles
		cMu = 0;
		cPi = 0;
		cMuF = 0;
		cPiF = 0;
		nProton = 0;

		//Any process is analysed
		GHepParticle * neu = event.Probe();
		TLorentzVector & neu4vec = *(neu->P4());
		double total_xs = (XS_graph == 0) ? 1 : XS_graph->Eval(neu->E());
		hinProbe->Fill(neu4vec.E());

		GHepParticle * p = 0;
		TIter event_iter(&event);

		while((p=dynamic_cast<GHepParticle *>(event_iter.Next())))
		{
			if (p->Pdg() == 13 ) //Muon
			{
				++cMu;
				if (p->Status() == 1)
				{
					++cMuF;
					mu4 = *(p->P4());
				}
			}
			if (p->Pdg() == 211) //Pion
			{
				++cPi;
				if (p->Status() == 1)
				{
					++cPiF;
					pi4 = *(p->P4());
				}
			}
			if (p->Pdg() == 2212) //Proton
				++nProton;
//			if (cMu*cPi == 1)
//				break;
		}

		std::cout << "All\tMuons " << cMu << "\tPions " << cPi << "\tBoth " << cMu*cPi << std::endl;
		std::cout << "Fin\tMuons " << cMuF << "\tPions " << cPiF << "\tBoth " << cMuF*cPiF << std::endl;
		if (cMuF*cPiF == 1)
		{
			++total;
			NEnergy = mu4.E()+pi4.E();

			muTheta = mu4.Angle(TVector3(0,0,1));	//Angle wrt to z
			piTheta = pi4.Angle(TVector3(0,0,1));  	//Angle wrt to z
			Angle = pi4.Angle(mu4.Vect());		//Separation between mu and pi

			TLorentzVector Tot4(mu4+pi4);
			InvMass = Tot4.M();
			NAngle = Tot4.Angle(TVector3(0,0,1));

			hNEnergy->Fill(NEnergy, total_xs);
			hInvMass->Fill(InvMass, total_xs);
			hnnAngle->Fill(NAngle, total_xs);
			hAngleup->Fill(Angle, total_xs);

			Index = i;
			Trace->Fill();
		}

		mcrec->Clear(); // clear current mc event record
	}//end loop over events

	// Scale to take binning and nbr generated events into account
//	double rate = 1;
//	double rate = 10 *  1 /nbrcoh;
//	h_pion_nrj->Scale(rate/0.01);
//	h_pion_mom->Scale(rate/0.01);
//	h_pion_ang->Scale(rate);
//	h_Q2->Scale(rate/0.01);
	// close input GHEP event file
	outfile.Write();
	outfile.Close();

	file.Close();
	LOG("myAnalysis", pNOTICE)  << "Done!";
	std::cout << "Total " << total << std::endl;
	return 0;
}

//___________________________________________________________________
void GetCommandLineArgs(int argc, char ** argv)
{
	LOG("myAnalysis", pINFO) << "Parsing commad line arguments";
	genie::CmdLnArgParser parser(argc,argv);

  	// get GENIE event sample
	if( parser.OptionExists('i') ) {
		LOG("myAnalysis", pINFO) << "Reading event sample filename";
		gOptInpFilename = parser.ArgAsString('i');
	} else {
		LOG("myAnalysis", pFATAL) << "Unspecified input filename - Exiting";
		exit(1);
	}

  	// number of events to analyse
	if( parser.OptionExists('n') ) {
		LOG("myAnalysis", pINFO) << "Reading number of events to analyze";
		gOptNEvt = parser.ArgAsInt('n');
	} else {
		LOG("myAnalysis", pINFO)<< "Unspecified number of events to analyze - Use all";
		gOptNEvt = -1;
	}

	// Output root file
	if( parser.OptionExists('o') ) {
		LOG("myAnalysis", pINFO) << "Creating output ROOT file";
		gOptOutFilename = parser.ArgAsString('o');
	} else {
		LOG("myAnalysis", pFATAL) << "Unspecified output filename";
	}

	// Flux file used in root file
	if( parser.OptionExists('f') ) {
		LOG("myAnalysis", pINFO) << "Reading flux file";
		gOptFluxFile = parser.ArgAsString('f');
	} else {
		LOG("myAnalysis", pINFO) << "Unspecified flux filename";
	}

	// Splines used in root file
	if( parser.OptionExists('s') ) {
		LOG("myAnalysis", pINFO) << "Reading splines file";
		gOptSplineFile = parser.ArgAsString('s');
	} else {
		LOG("myAnalysis", pINFO) << "Unspecified splines filename";
	}
}
//_________________________________________________________________________________
