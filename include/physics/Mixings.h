#ifndef MIXINGS_H
#define MIXINGS_H

#include <array>
#include <cmath>

#include "physics/Flavours.h"

struct Mixing {

	Mixing(double ue = 1., double um = 1., double ut = 1.) : _mix{{ue, um, ut}} {};

	/*
	Mixing(std::initializer_list<Nu::Flavour> flvs) : _mix{{0., 0., 0.}} {
		for (const Nu::Flavour &flv : flvs)
			switch (flv) {
				case Nu::Eb:
				case Nu::E_:
					_mix[0] = 1.;
					break;
				case Nu::Mb:
				case Nu::M_:
					_mix[1] = 1.;
					break;
				case Nu::Tb:
				case Nu::T_:
					_mix[2] = 1.;
					break;
				default:
					break;
			}
	}
	*/

	Mixing &operator*=(double x) {
		_mix[0] *= x;
		_mix[1] *= x;
		_mix[2] *= x;

		return *this;
	}

	Mixing operator*(double x) const {
		Mixing rhs = *this;
		return rhs *= x;
	}

	double operator()(Nu::Flavour flv, int e = 1) const {
		switch (flv) {
			case Nu::EB:
			case Nu::E0:
				return Ue(e);
			case Nu::MB:
			case Nu::M0:
				return Um(e);
			case Nu::TB:
			case Nu::T0:
				return Ut(e);
			default:
				return 0.;
		}
	}

	double Ue(int e = 1) const {
		return fast_pow(_mix[0], e);
	}

	double Um(int e = 1) const {
		return fast_pow(_mix[1], e);
	}

	double Ut(int e = 1) const {
		return fast_pow(_mix[2], e);
	}

	// set mixing value to 1 or 0
	void Flatten() {
		for (double &m : _mix)
			m = bool(m);
	}

	private:
		std::array<double, 3> _mix;

		inline double fast_pow(double val, int e) const {
			switch (e) {
				case 1:
					return val;
				case -1:
					return 1./val;
				case 2:
					return val * val;
				default:
					return std::pow(val, e);
			}
		}
};

inline Mixing operator*(double x, const Mixing &mix) {
	return mix * x;
}


inline std::ostream & operator<<(std::ostream &os, const Mixing &mix) {
	return os << "Mixing<Ue=" << mix.Ue() << ", Um=" << mix.Um()
		  << ", Ut=" << mix.Ut() << ">";
}

#endif
