/*
 * Author: Tommaso Boschi
 */

#ifndef OPTIMIZATION_H
#define OPTIMIZATION_H

#include <vector>
#include <functional>

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
	std::vector<T> NelderMead(Func fn, size_t ndim, T err = 1.e-4,
			std::initializer_list<T> lh = {}, std::initializer_list<T> ls = {},
			bool verb = false)
	{
		int icount, ifault, numres;
		int kcount = 1000000, konvge = 100;
		std::vector<T> hint(lh), side(ls);

		T *xmin = new T[ndim]; // output
		T *start = new T[ndim];
		T *step = new T[ndim];
		for (size_t i = 0; i < ndim; ++i) {
			start[i] = hint.size() > i ? hint[i] : 0.5;
			step[i]  = side.size() > i ? side[i] : 0.5;
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

		if (verb) {
			std::cout << "NM <" << icount << ", " << ifault << ", " << numres << ">";
			if (kcount > 1e6)
				std::cout << " no convergence;";
			else
				std::cout << " convergence;";
			std::cout << " optimal value at:";
			for (auto d : ret)
				std::cout << "  " << d;
			std::cout << " = " << fn(xmin) << "\n";
		}


		delete xmin;
		delete start;
		delete step;

		return ret;
	}
}

#endif

