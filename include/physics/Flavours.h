#ifndef FLAVOURS_H
#define FLAVOURS_H

#include <string>

struct Nu
{
	enum Flavour {
		E0 = 0,
		electron = E0,
		M0 = 1,
		muon = M0,
		T0 = 2,
		tau = T0,
		EB = 3,
		antielectron = EB,
		MB = 4,
		antimuon = MB,
		TB = 5,
		antitau = TB,

		_undefined = -1
	};

	static Flavour fromString(std::string flv) {
		if (flv == "nuEB")
			return Flavour::EB;
		else if (flv == "nuMB")
			return Flavour::MB;
		else if (flv == "nuTB")
			return Flavour::TB;
		else if (flv == "nuE0")
			return Flavour::E0;
		else if (flv == "nuM0")
			return Flavour::M0;
		else if (flv == "nuT0")
			return Flavour::T0;
		else
			throw std::invalid_argument("Nu: unknown flavour \"" + flv + "\"");
	}
	
	static Flavour fromPDG(int pdg) {
		switch (pdg) {
			case -12:
				return Flavour::EB;
			case -14:
				return Flavour::MB;
			case -16:
				return Flavour::TB;
			case 12:
				return Flavour::E0;
			case 14:
				return Flavour::M0;
			case 16:
				return Flavour::T0;
			default:
				return _undefined;
		}
	}

	static std::string toString(int pdg) {
		return toString(fromPDG(pdg));
	}

	static std::string toString(Flavour flv) {
		switch (flv) {
			case Flavour::EB:
				return "nuEB";
			case Flavour::MB:
				return "nuMB";
			case Flavour::TB:
				return "nuTB";
			case Flavour::E0:
				return "nuE0";
			case Flavour::M0:
				return "nuM0";
			case Flavour::T0:
				return "nuT0";
			default:
				return "";
		}
	}

	static bool isParticle(Flavour flv) {
		switch (flv) {
			case Flavour::E0:
			case Flavour::M0:
			case Flavour::T0:
				return true;
			case Flavour::EB:
			case Flavour::MB:
			case Flavour::TB:
				return false;
			default:
				return false;
		}
	}

	static bool isAntiparticle(Flavour flv) {
		return !isParticle(flv);
	}

	static std::vector<Nu::Flavour> Part() {
		return { Flavour::E0, Flavour::M0, Flavour::T0 };
	}

	static std::vector<Nu::Flavour> Anti() {
		return { Flavour::EB, Flavour::MB, Flavour::TB };
	}

	static std::vector<Nu::Flavour> All() {
		return { Flavour::E0, Flavour::M0, Flavour::T0,
			 Flavour::EB, Flavour::MB, Flavour::TB };
	}

	//static constexpr const Flavour Part[] = { Flavour::E0, Flavour::M0, Flavour::T0 };
	//static constexpr const Flavour Anti[] = { Flavour::EB, Flavour::MB, Flavour::TB };
	//static constexpr const Flavour All[]  = { Flavour::E0, Flavour::M0, Flavour::T0,
	//					  Flavour::EB, Flavour::MB, Flavour::TB };
};

inline std::ostream & operator<<(std::ostream &os, const Nu::Flavour &flv) {
	return os << Nu::toString(flv);
}


#endif
