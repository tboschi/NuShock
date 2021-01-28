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

#include "tools/CardDealer.h"
#include "physics/Const.h"

#include "TLorentzVector.h"
#include "TGenPhaseSpace.h"

#include "cuba.h"

#include "LHAPDF/PDF.h"
#include "LHAPDF/AlphaS.h"

// integrand wrapper for member functions
// to be used, userdata should be passed to vegas as a void*
// this type makes cuba happy (it is typedef CUBA integrand_t)
template <class T>
int member_integral(const int *ndim, const double x[], const int *ncomp, double f[], void *userdata)
{
	T &obj = *(static_cast<T*>(userdata));
	f[0] = obj(x);	// use function call operator

	return 1;
}

class OpenQQ
{
	public:
		OpenQQ(char *cardname);
		OpenQQ(const std::string &cardname);
		OpenQQ(CardDealer *cd);
		~OpenQQ();

		void Init(const CardDealer &cd);

		inline void SetCMEnergy(double cme) {
			_cme = std::pow(cme, 2);
		}

		//double dXSdQ2_qq();
		//double dXSdQ2_gg();
		double dXSdOmega_qq();
		double dXSdOmega_gg();

		double operator()(const double *x);
		double Integrand(double x1, double x2, double omega);
		double Integrate(double &error, double &chi2prob);

	private:
		LHAPDF::PDF *_probe_PDF, *_target_PDF;
		size_t _probe_A, _probe_Z, _target_A, _target_Z;

		int _valquark;
		double _aS, _mf2, _m2, _cme, _s, _t, _u;

		// physics parameters that can be set at execution time
		double _re_scale, _fac_scale, _ps_cut;

		// VEGAS parameters
		double _err_rel, _err_abs;
		int _min_evals, _max_evals, _start_evals, _inc_evals, _batch_evals;

};

#endif
