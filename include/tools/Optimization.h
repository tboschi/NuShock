/*
 * Author: Tommaso Boschi
 */

#ifndef OPTIMIZATION_H
#define OPTIMIZATION_H

// nelmer mead implementation
#include "tools/asa047.hpp"

// functions return x value at min/max
namespace Optimization
{
	// 1D func, specific blueprint of function required
	template <typename T, typename Func>
	T GoldenRatio(Func fn, T err = 1.e-6)
	{
		T gr = (1 + std::sqrt(5)) / 2.0;

		T s = 0.0;		//start point
		T e = 1.0;		//end point

		while (e - s > err)
		{
			T a = e - (e - s)/gr;
			T b = s + (e - s)/gr;

			if (fn(a) < fn(b))
				e = b;
			else
				s = a;
		}

		return e < s ? s : e;
	}

	// ndim function
	template <typename T, typename Func>
	std::vector<T> NelderMead(Func fn, int ndim, T err = 1.e-6)
	{
		int icount, ifault, numres;
		int kcount = 100, konvge = 10;

		T *xmin = new T[ndim]; // output
		T *start = new T[ndim];
		T *step = new T[ndim];
		for (int i = 0; i < ndim; ++i) {
			start[i] = 0.5;
			step[i] = 0.1;
		}

		// result
		T yval;

		// minimizing function with Nelder Mead 
		nelmin (fn, ndim, start, xmin, &yval, 
			err, step, konvge, kcount, 
			&icount, &numres, &ifault);

		// to maximize is defined negative to find minimum
		// maximum is given by inverting minimum
		std::vector<T> ret(xmin, xmin + ndim);

		delete xmin;
		delete start;
		delete step;

		return ret;
	}
}

#endif

