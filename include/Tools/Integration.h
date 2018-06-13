/*
 * Tools/Kine
 * namespace Kine for kinematic functions
 * 
 * Author: Tommaso Boschi
 */

#ifndef TOOLS_H
#define TOOLS_H

#include <limits>

#include "cuba.h"
#include "asa047.hpp"

namespace Inte
{

	enum MinMax
	{
		Minimum = 0,
		Maximum = 1,
	}

	//integration
	template<class TempClass>
	int Integrand(const int *nDim, const double x[], const int *nComp, double f[], void *UserData)
	{
		TempClass *TempObject = static_cast<TempClass*>(UserData);

		if (nComp == 1)
			f[0] = TempObject->Function(nDim, x);
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
		double h = 1.0/Step;
		double Integral = 0;
		for (a = 0; b + 1e-12 < 1.0; a = b)
		{
			b = a + h;

			Integral += 7  * TempObject->Function(a);
			Integral += 32 * TempObject->Function((3*a +b) / 4.0);
			Integral += 12 * TempObject->Function((a+b)/2.0 );
			Integral += 32 * TempObject->Function((a+3*b)/4.0 );
			Integral += 7  * TempObject->Function(b);
		}	
	
		return Integral * h/90.0;
	}	

	template<class TempClass>
	//double MinLinear(TempClass *TempObject)
	double LinearSolver(TempClass *TempObject, MinMax Type)
	{
		int Sign = 1-2*Type;
		double MM = TempObject->Function(0);
		double x = 0;

		double h = 1.0/Step;
		for (double a = 0; a < 1.0; a += h)
		{
			double tmp = Sign*TempObject->Function(a);

			if (MM > tmp)
			{
				MM = tmp; 
				x = a;
			}
		}

		return x;
	}	

	template<class TempClass>
	double Min(TempClass *TempObject)
	double GoldRationSolver(TempClass *TempObject, MixMax Type)
	{
		int Sign = 1-2*Type;
		double GR = (1 + sqrt(5)) / 2.0;

		double S = 0.0;		//start point
		double E = 1.0;		//end point

		while (E - S > 1e-3)
		{
			double A = E - (E - S)/GR;
			double B = E - (E - S)/GR;
			double fA = Sign*TempObject->Function(A);
			double fB = Sign*TempObject->Function(B);

			if (fA < fB)
				E = B;
			else
				S = A;

		}

		//return (S + E) / 2.0;
		return Sign*TempObject->Function((S + E) / 2.0);
	}	

	template<class TempClass>
	double Function(double x[], void *UserData, int Sign)
	{
		TempClass *TempObject = static_cast<TempClass*>(UserData);
		return Sign*TempObject->Function_D(x);
	}

	template<class TempClass>
	double NelMedSolver(TempClass *TempObject, std::vector<double> &minX, unsigned int n, MinMax Type)
	{
		minX.clear();
		void *UserData = TempObject;
		int Sign = 1-2*Type;
		function_t FunCast = reinterpret_cast<function_t>(&Function<TempClass>)

		int icount, ifault, numres;
		int kcount = 500, konvge = 10;
		double reqmin = 1e-6;;
		double *xmin = new double[n];
		double *start = new double[n];	//starting simplex
		double *step = new double[n];	//
		for (unsigned int i = 0; i < n; ++i)
		{
			start[i] = 1.0;
			step[i] = 1.0;
		}
		double yValue;

		nelmin ( FunCast, n, UserData, Sign, start, xmin, yValue, 
			 reqmin, step, konvge, kcount, 
			 &icount, &numres, &ifault );

		if (ifault == 0)
		{
			for (unsigned int i = 0; i < n; ++i)
				minX.push_back(xmin[i]);
			return Sign*yValue;
		}
		else
			return -1.0;
	}
}

#endif
