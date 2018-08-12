#include<iostream>
#include <fstream>

#include "Tools.h"
#include "Flux.h"
#include "TH1D.h"

int main(int argc, char** argv)
{
	std::string File(argv[1]);
	Flux *ff = new Flux(File);

	std::ofstream Out(argv[2]);

	ff->Scale(1e20);
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
		else if (hh = ff->Get(Flux::Pion))
		 	Out << hh->GetBinContent(i) << "\t";
		else if (hh = ff->Get(Flux::PPion))
		 	Out << hh->GetBinContent(i) << "\t";
		else if (hh = ff->Get(Flux::Kaon))
		 	Out << hh->GetBinContent(i) << "\t";
		else if (hh = ff->Get(Flux::Kaon0))
		 	Out << hh->GetBinContent(i) << "\t";
		else if (hh = ff->Get(Flux::Charm))
		 	Out << hh->GetBinContent(i) << "\t";
		else if (hh = ff->Get(Flux::Muon))
		 	Out << hh->GetBinContent(i) << "\t";
		else if (hh = ff->Get(Flux::TauE))
		 	Out << hh->GetBinContent(i) << "\t";
		else if (hh = ff->Get(Flux::TauM))
			Out << hh->GetBinContent(i) << "\t";

		Energy += Width;
	}
		

	delete ff;
	Out.close();

	return 0;
}
