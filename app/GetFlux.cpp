#include <iostream>
#include <fstream>

#include "TFile.h"
#include "TH1.h"

#include "tools.h"
#include "flux.h"


int main(int argc, char** argv)
{
	std::string file(argv[1]);
	Flux fx(file);
	fx.Scale(1.0e20 / 574.0 / 574.0);
	fx.Add();

	std::vector<Flux::Hist> fluxType;

	fluxType.push_back(Flux::Total);
	fluxType.push_back(Flux::TauE);
	fluxType.push_back(Flux::TauM);
	fluxType.push_back(Flux::Pion);
	fluxType.push_back(Flux::Kaon);
	fluxType.push_back(Flux::Kaon0);
	fluxType.push_back(Flux::Muon);
	fluxType.push_back(Flux::Charm);
	fluxType.push_back(Flux::PPion);

	std::ofstream out(argv[2]);
	for (int i = 1; i < fx.BinNumber()+1; ++i)
	{
		out << fx.Get(Flux::Total)->GetBinLowEdge(i) << "\t";
		for (int name = 0; name < fluxType.size(); ++name)
			if (fx.Get(fluxType[name]))
				out << fx.Get(fluxType[name])->GetBinContent(i) << "\t";
		out << std::endl;
	}
	out.close();

	return 0;
}
