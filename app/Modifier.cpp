#include <iostream>
#include <cstring>
#include <fstream>

#include "detector/Flux.h"
#include "detector/Driver.h"

int main(int argc, char **argv) {

	using modifier = std::vector<std::array<double, 5> >;
	std::map<Flux::Parent, modifier> modifs;
	for (int f = 1; f < argc; ++f) {
		std::cout << "file " << argv[f] << "\n";
		std::string file0(argv[f]);
		std::string fileB = file0;
		fileB.replace(file0.find("T0"), 2, "TB");

		auto c0 = Flux::fromROOT(file0);
		auto cB = Flux::fromROOT(fileB);
		c0.insert(std::make_move_iterator(cB.begin()),
			  std::make_move_iterator(cB.end()));

		file0.erase(file0.find(".root"));
		file0.erase(0, file0.find("T0_")+3);

		double mass = std::strtod(file0.c_str(), NULL) * 1.e-3;
		std::cout << file0 << "\tMass " << mass << "\n";

		for (const auto &c : c0) {
			std::array<double, 5> mod;
			mod[0] = mass;
			if (c.second->FindFirstBinAbove() >= 0.)
				mod[1] = c.second->GetBinLowEdge(c.second->FindFirstBinAbove());
			else
				mod[1] = 0.;
			if (c.second->FindLastBinAbove() >= 0.)
				mod[2] = c.second->GetBinLowEdge(c.second->FindLastBinAbove()+1);
			else
				mod[2] = 0.;
			mod[3] = c.second->GetMaximum();
			mod[4] = c.second->Integral("WIDTH");

			modifs[c.first].push_back(std::move(mod));
		}
	}

	for (const auto &mod : modifs) {
		std::string name = "data/raw_" + Flux::toString(mod.first) + ".dat";
		std::ofstream out(name.c_str());
		for (auto &m : mod.second)
			out << m[0] << "\t" << m[1] << "\t" << m[2]
			    << "\t" << m[3] << "\t" << m[4] << "\n";
		out.close();
	}

	return 0;
}
