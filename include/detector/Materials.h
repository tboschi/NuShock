#ifndef MATERIAL_H
#define MATERIAL_H

#include <string>
#include <cmath>
#include "physics/Const.h"

struct Material {

	enum Name {
		undefined, 

		Air,

		LAr,
		LiquidArgon = LAr,

		GasAr,
		GasseousArgon = GasAr,

		Fe,
		Iron = Fe,

		Pb,
		Lead = Pb,
	};

	static Material::Name fromString(std::string name) {
		if (name == "Air")
			return Name::Air;
		else if (name == "LAr" || name == "LiquidArgon")
			return Name::LAr;
		else if (name == "GasAr" || name == "GasseousArgon")
			return Name::GasAr;
		else if (name == "Fe" || name == "Iron")
			return Name::Fe;
		else if (name == "Pb" || name == "Lead")
			return Name::Pb;
		else
			return Name::undefined;
	}

	static std::string toString(Material::Name name) {
		switch (name) {
			case Name::Air:
				return "Air";
			case Name::LAr:
				return "LAr";
			case Name::GasAr:
				return "GasAr";
			case Name::Fe:
				return "Fe";
			case Name::Pb:
				return "Pd";
			default:
				return "undefined";
		}
	}

	static double Density(Material::Name name)
	{
		switch (name)
		{
			case Air:
				return 1.2e-3; // g / cm^3
			case LAr:
				return 1.3945; // g / cm^3
			case GasAr:
				return 0.1020; // g / cm^3
			case Fe:
				return 7.874; // g / cm^3
			case Pb:
				return 11.34; // g / cm^3
			default:
				return 0.;
		}
	}

	// atomic mass
	static double A(Material::Name name)
	{
		switch (name)
		{
			case Name::Air:
				return 28.965;
			case Name::LAr:
			case Name::GasAr:
				return 40.;
			case Name::Fe:
				return 56.;
			case Name::Pb:
				return 207.;
			default:
				return 0;
		}
	}

	// atomic number
	static double Z(Material::Name name)
	{
		switch (name)
		{
			case Name::Air:
				return 14.453;
			case Name::LAr:
			case Name::GasAr:
				return 18.;
			case Name::Fe:
				return 26.;
			case Name::Pb:
				return 82.;
			default:
				return 0.;
		}
	}

	static double IonizationEnergy(Material::Name name)
	{
		switch (name)
		{
			case Name::Air:
				return 85.7;	// eV
			case Name::LAr:
			case Name::GasAr:
				return 188.;	// eV
			case Name::Fe:
				return 286.;	// eV
			case Name::Pb:
				return 823.;	// eV
			default:
				return 0.;	// eV
		}
	}

	//assuming same for positron and electron
	static double CriticalEnergy(Material::Name name)
	{
		switch (name)
		{
			case Name::Air:
				return 0.08792;	//GeV
			case Name::LAr:
				return 0.03284;	//GeV
			case Name::GasAr:
				return 0.03803;	//GeV
			case Name::Fe:
				return 0.02168;	//GeV
			case Name::Pb:
				return 0.00743;	//GeV
			default:
				return 0.; // GeV
		}
	}

	// in cm
	static double RadiationLength(Material::Name name)
	{
		double rad = 0.;
		switch (name)
		{
			case Name::Air:
				rad = 36.62; // g cm-2
				break;
			case LAr:
				rad = 19.55; // g cm-2
				break;
			case GasAr:
				rad = 19.55; // g cm-2
				break;
			case Fe:
				rad = 13.84; // g cm-2
				break;
			case Pb:
				rad = 6.37; // g cm-2
				break;
			default:
				return 0.; // g cm-2
		}

		return rad / Density(name);	// cm
	}

	// in cm
	static double NuclearRadiationLength(Material::Name name)
	{
		double rad = 0.;
		switch (name)
		{
			case Name::Air:
				rad = 90.1; // g cm-2
				break;
			case LAr:
				rad = 119.7; // g cm-2
				break;
			case GasAr:
				rad = 119.7; // g cm-2
				break;
			case Fe:
				rad = 132.1; // g cm-2
				break;
			case Pb:
				rad = 199.6; // g cm-2
				break;
			default:
				return 0.; // g cm-2
		}

		return rad / Density(name);	// cm
	}

	// gives Bethe energy loss dE/dx in GeV/m (0.1*)
	static double Bethe(Material::Name name, double mass, double beta, double gamma)
	{
		double e2m = Const::MElectron / mass;
		//electron mass in MeV -> *1000
		double Wmax = (2000. * Const::MElectron * std::pow(beta * gamma, 2))
			    / (1. + e2m * (2. * gamma + e2m));

		double logArg = 2000. * Const::MElectron * std::pow(beta * gamma, 2)
			      * Wmax / (1.e-12 * IonizationEnergy(name)); 	//Everything in MeV

		beta *= beta;
		// constant from PDG MeV mol-1 cm2
		return 0.0307075 * Density(name) * Z(name) / (A(name) * beta)
		       * (0.5 * std::log(logArg) - beta);
	}

};

#endif
