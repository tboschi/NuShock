#include<iostream>
#include <fstream>

#include "Tools.h"
#include "Flux.h"
#include "TH1D.h"

int main(int argc, char** argv)
{
	std::string File(argv[1]);
	Flux *ff = new Flux(File);
	ff->Scale(1e20);

	std::ofstream Out(argv[2]);

	std::string Det(argv[3]);
	Detector * TheBox = new Detector(Det);
	ff->Scale(1.0/pow(TheBox->Get("Baseline"), 2));

	ff->Add();

	TH1D *hh;
	double Energy = ff->RangeStart();
	double End   = ff->RangeEnd();
	double Width = ff->BinWidth();
	double BinN  = ff->BinNumber();
	for (unsigned int i = 0; i < BinN; ++i)
	{
		Out << Energy << "\t";
		if (hh = ff->Get(Flux::Total))
			Out << hh->GetBinContent(i) << "\t";
		if (hh = ff->Get(Flux::Pion))
		 	Out << hh->GetBinContent(i) << "\t";
		if (hh = ff->Get(Flux::PPion))
		 	Out << hh->GetBinContent(i) << "\t";
		if (hh = ff->Get(Flux::Kaon))
		 	Out << hh->GetBinContent(i) << "\t";
		if (hh = ff->Get(Flux::Kaon0))
		 	Out << hh->GetBinContent(i) << "\t";
		if (hh = ff->Get(Flux::Charm))
		 	Out << hh->GetBinContent(i) << "\t";
		if (hh = ff->Get(Flux::Muon))
		 	Out << hh->GetBinContent(i) << "\t";
		if (hh = ff->Get(Flux::TauE))
		 	Out << hh->GetBinContent(i) << "\t";
		if (hh = ff->Get(Flux::TauM))
			Out << hh->GetBinContent(i) << "\t";

		Out << std::endl;
		Energy += Width;
	}
		

	delete ff;
	delete TheBox;
	Out.close();

	return 0;
}
