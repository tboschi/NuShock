#include "physics/OpenQQ.h"

OpenQQ::OpenQQ(char * cardname)
{
	CardDealer cd(cardname);
	Init(cd);
}

OpenQQ::OpenQQ(const std::string &cardname)
{
	CardDealer cd(cardname);
	Init(cd);
}

OpenQQ::OpenQQ(CardDelaer *cd)
{
	Init(*cd);
}

void OpenQQ::Init(const CardDealer &cd)
{
	std::string probe_pdf_name;
	if (!cd.Get("probe_pdf", probe_pdf_name))
		throw std::invalid_argument("No PDF for probe particle set\n");
	if (!cd.Get("probe_mass", _probe_A))
		_probe_A = 1;
	if (!cd.Get("probe_atomic", _probe_Z))
		_probe_Z = 1;

	std::string target_pdf_name;
	if (!cd.Get("target_pdf", target_pdf_name))
		throw std::invalid_argument("No PDF for target particle set\n");
	if (!cd.Get("target_mass", _target_A))
		_target_A = 1;
	if (!cd.Get("target_atomic", _target_Z))
		_target_Z = 1;

	_probe_PDF  = LHAPDF::mkPDF(probe_pdf_name);
	_target_PDF = LHAPDF::mkPDF(target_pdf_name);

	//_partonPS = new TGenPhaseSpace;

	// care mostly about charm
	char quark;
	if (!cd.Get("valence_quark", quark))
		quark = 'c';

	double valmass;
	switch(quark)
	{
		case 'u':
			_valquark = 1;
			valmass = Const::MQuarkU;
			break;
		case 'd':
			_valquark = 2;
			valmass = Const::MQuarkD;
			break;
		case 's':
			_valquark = 3;
			valmass = Const::MQuarkS;
			break;
		case 'c':
			_valquark = 4;
			valmass = Const::MQuarkC;
			break;
		case 'b':
			_valquark = 5;
			valmass = Const::MQuarkB;
			break;
		case 't':
			_valquark = 6;
			valmass = Const::MQuarkT;
			break;
		default:
			throw std::invalid_argument("Unknown quark type\n");
	}

	// useful variable to have
	_m2 = std::pow(valmass, 2);
	_cme = -1;

	// QCD scale
	if (!cd.("renormalization_scale", _re_scale))
		_re_scale = 1.6;
	if (!cd.("factorization_scale", _fac_scale))
		_fac_scale = 2.1;

	// VEGAS parameters
	if (!cd.("relative_error"), _err_rel)
		_err_rel = 1.0e-6;	//relative error for each component
	if (!cd.("absolute_error"), _abs_rel)
		_err_abs = 1.0e-9;	//absolute error
	if (!cd.("min_evaluations"), _min_evals)
		_min_evals = 1e3;	//minimum number of evaluation
	if (!cd.("max_evaluations"), _max_evals)
		_max_eval = 1e6;		//maximum number of evaluation
	if (!cd.("start_evaluations"), _start_evals)
		_start_evals = 10;
	if (!cd.("increment_evaluations"), _inc_evals)
		_inc_evals = 10;
	if (!cd.("batch_evaluations"), _batch_evals)
		_batch_evals = 1000;
}

/*
// set partons such that CME is sqrt(s)
void OpenQQ::SetPartons(double s)	//probe on target, fixed by experiment and constant
{
	_parton1 = TLorentzVector(0, 0,  std::sqrt(s) / 2., std::sqrt(s) / 2.);
	_parton2 = TLorentzVector(0, 0, -std::sqrt(s) / 2., std::sqrt(s) / 2.);

	double mass[2] = {_m2, _m2};
	TLorentzVector S = _parton1 + _parton2;

	if (_partonPS->SetDecay(S, 2, Mass)) {
		SetFromPS();
	}
}

void OpenQQ::SetFromPS()
{
	_partonPS->Generate();

	_quark3.SetVect(_partonPS->GetDecay(0)->Vect());
	_quark3.SetE(_partonPS->GetDecay(0)->E());

	_quark4.SetVect(_partonPS->GetDecay(1)->Vect());
	_quark4.SetE(_partonPS->GetDecay(1)->E());

	TLorentzVector S1 = *Part1 + *Part2;
	TLorentzVector S2 = *Quark3 + *Quark4;
	std::cout << "S " << S1.M2() << "\t" << S2.M2() << std::endl;
	SetCMEnergy(S2.M2());

	TLorentzVector T1 = *Part1 - *Quark3;
	TLorentzVector T2 = *Part2 - *Quark4;
	std::cout << "T " << T1.M2() << "\t" << T2.M2() << std::endl;
	SetQ2(-T1.M2());
}
*/


//XSections

double OpenQQ::dXSdQ2_qq()	//differential cross section (dXS/dQ2) for qq_ into qq_
{
	double aS = _probe_PFD->alphaQ2(_m2 * std::pow(_re_scale, 2));
	return  4. * aS * Const::pi / (9. * std::pow(_s, 4)) *
		(_t * _t + _u * _u + 2. * _m2 * _s) ;
}

double OpenQQ::dXSdQ2_gg()	//differential cross section (dXS/dQ2) for gg into qq_
{
	double aS = _probe_PFD->alphaQ2(_m2 * std::pow(_re_scale, 2));
	return  aS * Const::pi / (_s * _s) * (4./3. - 3. _t * _u / (_s * _s)) / 8. *
		(_t / _u + _u / _t + 4. * _m2 _s / (_t * _u) * (1. - _m2 * _s / (_t * _u) ) );
}

double OpenQQ::dXSdOmega_qq()	//differential cross section (dXS/dOmega) for qq_ into qq_
{
	double aS = _probe_PFD->alphaQ2(_m2 * std::pow(_re_scale, 2));
	return aS / (9. * std::pow(_s, 3)) * std::sqrt(1. - 4 * _m2 / _s) *
		(_t * _t + _u * _u + 2 * _m2 * _s);
	}
}

double OpenQQ::dXSdOmega_gg()	//differential cross section (dXS/dOmega) for gg into qq_
{
	double aS = _probe_PFD->alphaQ2(_m2 * std::pow(_re_scale, 2));
	return aS / (32. * _s) * std::sqrt(1. - 4. * _m2 / _s) *
		(6. * _t * _u / _s - _m2 *  (_s - 4 * _m2)  / (3. * _t * _u) + 
		 4. * (_t * _u - 2 * _m2 * (2. * _m2 + _t)) / (3. * _t * _t) +
		 4. * (_t * _u - 2 * _m2 * (2. * _m2 + _u)) / (3. * _u * _u) -
		 3. * (_t * _u - _m2 * (_u - _t)) / (_s * _t) -
		 3. * (_t * _u - _m2 * (_t - _u)) / (_s * _u) );
}

//for vegas integeration
//double void OpenQQ::Integrand(const int *nDim, const double *x, const int *nComp, double f[])
double OpenQQ::operator()(const double *input)
{
	// total integration
	double tau0 = 4 * _m2 / _cme;	//\hat{s} / s

	// normalized variables in thei integration range
	// because VEGAS inputs are [0:1]
	double x1 = (1 - tau0) * input[0] + tau0;
	double x2 = (1 - tau0 / x1) * input[1] + tau0 / x1;
	double omega = 2 * input[2] - 1;

	// set mandelstam variables use in xsec formulae
	_s = _cme * x1 * x2;
	_t = _s / 2. * (1. - omega * std::sqrt(1. - 4. * _m2 / _s)); // it is actually t - m2
	_u = _s - _t;	// it is actually u - m2

	//jacobian	pdf integration		* phasespace
	double jac = (1 - tau0) * (1 - tau0/x1) * 4.0*Const::pi;
	return jac * Integrand(x1, x2, omega);
}

double OpenQQ::Integrand(double x1, double x2, double omega)
{
	if (std::abs(omega) > _ps_cut)
		return 0.0;

	// check what 2.1 is! and get from card file
	double mf2 = _m2 * std::pow(_fac_scale, 2);
	if (_probe_PDF->inRangeQ2(mf2) && _target_PDF->inRangeQ2(mf2) && 
	    _probe_PDF->inRangeX(x1) && _target_PDF->inRangeX(x1) && 
	    _probe_PDF->inRangeX(x2) && _target_PDF->inRangeX(x2) ) {
		//gluon component
		double pdf_gg = _probe_PDF->xfxQ2(0, x1, mf2) * _target_PDF->xfxQ2(0, x2, mf2)
			      + _probe_PDF->xfxQ2(0, x2, mf2) * _target_PDF->xfxQ2(0, x1, mf2);

		double pdf_qq = 0.0;
		for (int i = 1; i < _valquark; ++i) {
			//q/p and q_/A
			pdf_qq += _probe_PDF->xfxQ2( i, x1, mf2) * _target_PDF->xfxQ2(-i, x2, mf2)
				+ _probe_PDF->xfxQ2( i, x2, mf2) * _target_PDF->xfxQ2(-i, x1, mf2);

			//q_/p and q/A
			pdf_qq += _probe_PDF->xfxQ2(-i, x1, mf2) * _target_PDF->xfxQ2( i, x2, mf2);
				+ _probe_PDF->xfxQ2(-i, x2, mf2) * _target_PDF->xfxQ2( i, x1, mf2);
		}

		return (pdf_gg * dXSdOmega_gg() + pdf_qq * dXSdOmega_qq()) / (x1 * x2);
	}

	return 0.0;
}

// total integration of opencc xsec using VEGAS
double OpenQQ::Integrate(double &error, double &chi2prob)
{

	char *state = NULL;
	void *spin = NULL;

	//cuba wants an integrand_t type passed to Vegas
	//can pass this pointer as user data and recast it in member integral
	integrand_t intcast = static_cast<integrand_t>(&member_integral<OpenQQ>); 
	void *ud = static_cast<void*>(this);

	//output
	int trial, fail;
	double result;

	Vegas(3, 1, intcast, ud, 1, 	//ndim, ncomp, integrand_t, userdata, nvec
	      _err_rel, _err_abs, 0, 0, 			//epsrel, epsabs, verbosity, seed
	      _min_eval, _max_eval, _start_eval, _inc_eval, _batch_eval,
	      0, 0, 0,				//gridno, statefile, spin
	      &trial, &fail, &result, &error, &chi2prob);

	return Const::GeV2ub * _target_A * result;
}
