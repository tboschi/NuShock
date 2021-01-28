/*
 * Author: Tommaso Boschi
 */

#ifndef INTEGRATION_H
#define INTEGRATION_H

#include <cmath>
#include <functional>


namespace Integration
{
	template<size_t N, typename T, typename Func>
	struct _integrand {
		T v;
		const Func &fn;

		template<typename... Args>
		T operator()(Args... rest) const {
			return fn(v, rest...);
		}
	};

	template <size_t N, typename T, typename Func,
		 typename std::enable_if<(N > 1), bool>::type = true>
	T Boole(Func fn, size_t steps = 1000) {
		auto subfn = [&](T v) -> T {
			// if c++11
			auto subint = _integrand<N, T, Func> {v, fn};
			// if c++14
			//auto subint = [v, &fn](auto... rest) {
				//return fn(v, rest...);
			//};

			return Boole<N-1, T>(subint, steps);
		};

		return Boole<1, T>(subfn, steps);
	}

	template <size_t N, typename T, typename Func,
		 typename std::enable_if<(N == 1), bool>::type = true>
	T Boole(Func fn, size_t steps = 1000) {
		T h = 1./steps;
		T res = 0.;

		T left = 7. * fn(0.);
		for (size_t s = 0; s < steps; ++s) {
			res += left;

			res += 32. * fn((s + .25) * h);
			res += 12. * fn((s + 0.5) * h);
			res += 32. * fn((s + .75) * h);

			left = 7. * fn((s + 1.) * h);

			res += left;
		}	

		return res * h / 90.;
	}


	/*
	// fn is a 1D function
	template <typename Functor>
	double Boole1D(Functor fn, size_t steps = 1000)
	{
		double h = 1./steps;
		double res = 0.;

		double left = 7. * fn(0.);
		for (size_t s = 0; s < steps; ++s) {
			res += left;

			res += 32. * fn((s + .25) * h);
			res += 12. * fn((s + 0.5) * h);
			res += 32. * fn((s + .75) * h);

			left = 7. * fn((s + 1.) * h);

			res += left;
		}	

		return res * h / 90.;
	}

	template <typename Functor>
	double Boole2D(Functor fn, size_t steps = 1000)
	{
		auto subfn = [&](double v) ->double {
				auto subint = std::bind(fn, v, std::placeholders::_1);
				return Boole1D(subint, steps);
			 };

		return Boole1D(subfn, steps);
	}

		double h = 1./steps;
		double res = 0.;

		double left = 7. * subfn(0.);
		for (size_t s = 0; s < steps; ++s) {
			res += left;

			res += 32. * subfn((s + .25) * h);
			res += 12. * subfn((s + 0.5) * h);
			res += 32. * subfn((s + .75) * h);

			left = 7. * subfn((s + 1.) * h);

			res += left;
		}	

		return res * h / 90.;
	}
	*/
}

#endif
