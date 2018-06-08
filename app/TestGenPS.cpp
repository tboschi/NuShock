#include <iostream>

#include "TGenPhaseSpace.h"
#include "TLorentzVector.h"
#include "TH1.h"
#include "TFile.h"
#include "TTree.h"

int main()
{
	TFile *ffile = new TFile("saveps.root", "RECREATE");
	TTree *tev = new TTree("tev", "tev");

	double E0, E1, E2, T0, T1, T2, P0, P1, P2;
	double px0, py0, pz0, px1, py1, pz1, px2, py2, pz2;
	double sep, w;

	tev->Branch("E0", &E0, "fE0/D");
	tev->Branch("E1", &E1, "fE1/D");
	tev->Branch("E2", &E2, "fE1/D");
	tev->Branch("T0", &T0, "fT0/D");
	tev->Branch("T1", &T1, "fT1/D");
	tev->Branch("T2", &T2, "fT2/D");
	tev->Branch("P0", &P0, "fP0/D");
	tev->Branch("P1", &P1, "fP1/D");
	tev->Branch("P2", &P2, "fP2/D");
	tev->Branch("px0", &px0, "fpx0/D");
	tev->Branch("py0", &py0, "fpy0/D");
	tev->Branch("pz0", &pz0, "fpz0/D");
	tev->Branch("px1", &px1, "fpx1/D");
	tev->Branch("py1", &py1, "fpy1/D");
	tev->Branch("pz1", &pz1, "fpz1/D");
	tev->Branch("px2", &px2, "fpx2/D");
	tev->Branch("py2", &py2, "fpy2/D");
	tev->Branch("pz2", &pz2, "fpz2/D");
	tev->Branch("sep", &sep, "fsep/D");
	tev->Branch( "w",  &w,  "fw/D");

	TGenPhaseSpace event;
	TLorentzVector neutrino(0.0, 0.0, 0.0, 1.0);
	std::cout << "Beta " << neutrino.Beta() << std::endl;
	double masses[3] = {0.0, 0.05, 0.05};
	event.SetDecay(neutrino, 3, masses);

	for (unsigned int n = 0; n < 1e5; ++n)
	{
		w = event.Generate();
		TLorentzVector *neut  = event.GetDecay(0);
		TLorentzVector *muon1 = event.GetDecay(1);
		TLorentzVector *muon2 = event.GetDecay(2);

		E0 = neut->E();
		E1 = muon1->E();
		E2 = muon2->E();
		T1 = neut->Theta();
		T1 = muon1->Theta();
		T2 = muon2->Theta();
		P1 = neut->Phi();
		P1 = muon1->Phi();
		P2 = muon2->Phi();

		px0 = neut->Px();
		py0 = neut->Py();
		pz0 = neut->Pz();
		px1 = muon1->Px();
		py1 = muon1->Py();
		pz1 = muon1->Pz();
		px2 = muon2->Px();
		py2 = muon2->Py();
		pz2 = muon2->Pz();

		sep = muon1->Angle(muon2->Vect());

		tev->Fill();
	}

	tev->Write();
	ffile->Close();

	return 0;
}
