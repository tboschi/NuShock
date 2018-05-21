/*
 * Tools/Kine
 * namespace Kine for kinematic functions
 * 
 * Author: Tommaso Boschi
 */

#ifndef TOOLS_H
#define TOOLS_H

#include "cuba.h"

namespace Inte
{
	//integration
	template<class TempClass>
	int Integrand(const int *nDim, const double x[], const int *nComp, double f[], void *UserData)
	{
		TempClass *TempObject = static_cast<TempClass*>(UserData);

		TempObject->Integrand(nDim, x, nComp, f);
	}

	template<class TempClass>
	double VegasIntegration(TempClass *TempObject, int nDim, int &Trial, int &Fail, double &Error, double &Chi2Prob)
	{
		//input
		double EpsRel = 1.0e-8;		//relative error for each component
		double EpsAbs = 1.0e-12;	//absolute error
		int MinEval = 1e5;		//minimum number of evaluation
		int MaxEval = 1e6;		//maximum number of evaluation
		int nStart = 1000;
		int nIncrease = 500;
		int nBatch = 1000;
		void *UserData = TempObject;
		char *state = NULL;
		void *spin = NULL;

		integrand_t IntCast = reinterpret_cast<integrand_t>(&Integrand<TempClass>); 
	
		//output
		double Integral;

		Vegas(nDim, 1, IntCast, UserData, 1, 	//ndim, ncomp, integrand_t, userdata, nvec
		      EpsRel, EpsAbs, 0, 0, 		//epsrel, epsabs, verbosity, seed
		      MinEval, MaxEval, nStart, nIncrease, nBatch,
		      0, state, spin,			//gridno, statefile, spin
		      &Trial, &Fail, &Integral, &Error, &Chi2Prob);
	
		return Integral;
	}

	template<class TempClass>
	double BooleIntegration(TempClass *TempObject)
	{
		double a = 0, b = 0;
		double h = 1.0/100.0;
		double Integral = 0;
		for (a = 0; b + 1e-12 < 1.0; a = b)
		{
			b = a + h;

			Integral += 7  * TempObject->Integrand(a);
			Integral += 32 * TempObject->Integrand((3*a +b) / 4.0);
			Integral += 12 * TempObject->Integrand((a+b)/2.0 );
			Integral += 32 * TempObject->Integrand((a+3*b)/4.0 );
			Integral += 7  * TempObject->Integrand(b);
		}	
	
		return Integral * h/90.0;
	}	
}

#endif
