#include<iostream>
#include <fstream>

#include "Tools.h"
#include "Physics.h"

int main()
{
	PhaseSpace *ThePS = new PhaseSpace();
	double mix[3] = {1.0, 1.0, 1.0};
	double mass = 0.4;
	ThePS->SetNeutrino(mass, mix, 1, 1, -1);
	ThePS->SetMass(mass);

	std::ofstream Out("nll.dat");

	double x = 0;
	double y = pow(Const::fMMuon/mass, 2);
	double z = pow(Const::fMElectron/mass, 2);

	for (double s = 0; s < 1; s += 0.01)
		for (double cos0 = -1; cos0 < 1; cos0 += 0.01)
			Out << s << "\t" << cos0 << "\t" << ThePS->M2_WW(s, cos0, x, y, z) << "\t" << ThePS->M2_WZ(s, cos0, x, y, z)  << std::endl;

	return 0;
}
