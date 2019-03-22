#include <iostream>

#include "TFile.h"
#include "TH1D.h"

int main(int argc, char** argv)
{
	TFile * OutF = new TFile("testroot.root", "RECREATE");

	TH1D *htest = new TH1D("htest", "htest", 100, -10, 10);

	htest->FillRandom("gaus");
	htest->Write();

	OutF->Close();

	return 0;
}
