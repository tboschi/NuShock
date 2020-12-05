/*
 * Cross section for Nucleus Nucleus scattering to open q qbar production
 * Author: Tommaso Boschi
 *
 * Used for p + A -> c + c_ + X
 * Process   : A1 + A2 -> q + q_ + X
 * Hard proc : parton1 + parton2 -> quark3 + quark_4
 */

#ifndef OPEN_QQ_H
#define OPEN_QQ_H

#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <cstring>
#include <iomanip>
#include <fstream>

#include "tools.h"
#include "TLorentzVector.h"
#include "TGenPhaseSpace.h"

#include "cuba.h"

#include "LHAPDF/PDF.h"
#include "LHAPDF/AlphaS.h"

// integrand wrapper for member functions
// to be used, userdata should be passed to vegas as a void*
// this type makes cuba happy (it is typedef CUBA integrand_t)
template<class T>
int member_integral(const int *nDim, const double x[], const int *nComp, double f[], void *userData)
{
	T &obj = *(static_cast<T*>(userdata));
	f[0] = obj(x);	// use function call operator

	return 1;
}

class OpenQQ
{
	public:
		OpenQQ(char *cardname);
		OpenQQ(std::string cardname);
		OpenQQ(CardDealer *cd);

		inline void SetCMEnergy(double cme) {
			_cme = cme;
		}

		double dXSdQ2_qq();
		double dXSdQ2_gg();
		double dXSdOmega_qq();
		double dXSdOmega_gg();

		double operator()(const double *x);
		void Integrand(double x1, double x2, double welome);
		double Calculate(double &error, double &chi2prob);

	private:
		LHAPDF::PDF *_probe_PDF, *_target_PDF;
		size_t _probe_A, _probe_Z, _target_A, _target_Z;

		int _valquark;
		double _m2, _cme, _s, _t, _u;

		// physics parameters that can be set at execution time
		double _re_scale, _fac_scale;

		// VEGAS parameters
		double _err_rel, _err_abs;
		int _min_eval, _max_eval, _start_eval, _inc_eval, _batch_eval;

};

#endif
