/*
 * Neutrino class
 * Author: Tommaso Boschi
 */

#ifndef NEUTRINO_H
#define NEUTRINO_H

#include <iostream>
#include <string>

#include "physics/Const.h"
#include "physics/Particle.h"

class Neutrino : public Particle
{
	public:
		// options seen as binary
		//<majorana/dirac>  <anti/particle> <un/polarized> <left/right>
		//     <1/0>             <1/0>          <1/0>         <1/0>
		//       b3                b2             b1            b0
		// e.g. b0000 == particle | dirac | left
		// if bit b3 is 1, bit b2 is ignored
		// if bit b1 is 1, bit b0 is ignored
		enum Options
		{
			left         = 0,
			right        = 1,
			polarised    = 0,
			unpolarized  = 2,
			particle     = 0,
			antiparticle = 4,
			dirac	     = 0,
			majorana     = 8
		};

		// default pdg code is nu_e
		Neutrino(double mass = 0., size_t opts = 0, int pdg = 12);

		bool operator==(const Neutrino & rhs) const;
		bool operator!=(const Neutrino & rhs) const;

		size_t GetOptions() const;
		void SetOptions(size_t opts);
		void AddOptions(size_t opts);
		
		int Helicity() const;
		static int Helicity(size_t opts) {
			if (opts & Neutrino::unpolarized)
				return 0;
			return 2 * int(opts & Neutrino::right) - 1;
		}

		bool IsDirac() const;
		static bool IsDirac(size_t opts) {
			return !IsMajorana(opts);
		}

		bool IsMajorana() const;
		static bool IsMajorana(size_t opts) {
			return (opts & Neutrino::majorana);
		}

		// always true if majorana
		bool IsParticle() const;
		static bool IsParticle(size_t opts) {
			return IsMajorana(opts) || !IsAntiparticle(opts);
		}

		// always true if majorana
		bool IsAntiparticle() const;
		static bool IsAntiparticle(size_t opts) {
			return IsMajorana(opts) || (opts & Neutrino::antiparticle);
		}

		// tells if neutrino is particle or antiparticle
		// in dirac sense. this can be set by last flag
		bool IsDiracParticle() const;
		static bool IsDiracParticle(size_t opts) {
			return !IsDiracAntiparticle(opts);
		}

		// always true if majorana
		bool IsDiracAntiparticle() const;
		static bool IsDiracAntiparticle(size_t opts) {
			return (opts & Neutrino::antiparticle);
		}

	private:
		size_t _opts;

	friend std::ostream & operator<<(std::ostream &os, const Neutrino &N);
};

#endif
