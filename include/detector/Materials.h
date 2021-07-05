/* Material specifications
 * add your own by:
 * 	creating new enum Name entry
 * 	extending fromString and toString
 * 	completing other properties accordingly
 *
 * All numbers are taken from https://pdg.lbl.gov/2011/AtomicNuclearProperties
 */

#ifndef MATERIAL_H
#define MATERIAL_H

#include <string>
#include <cmath>
#include <stdexcept>

#include "physics/Const.h"

namespace Material {

	enum class Name {
		LAr,
		argon_liquid = LAr,

		gasAr,
		argon_gas = gasAr,

		Fe,
		iron = Fe,

		Pb,
		lead = Pb,

		air,

		concrete,
	};

	Name fromString(std::string name);
	std::string toString(Material::Name name);

	double Density(Name name);

	// average atomic number over atmoic mass ratio
	double ZA(Name name);

	double IonizationEnergy(Name name);

	//assuming same for positron and electron
	double CriticalEnergy(Name name);

	// in cm
	double RadiationLength(Name name);

	// in cm
	double NuclearRadiationLength(Name name);

	// gives Bethe energy loss dE/dx in GeV/m (0.1*)
	double Bethe(Name name, double mass, double beta, double gamma);
}

inline std::ostream & operator<<(std::ostream &os, const Material::Name &mat) {
	return os << Material::toString(mat);
}

#endif
