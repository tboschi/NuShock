/*
 * Author: Tommaso Boschi
 */

#ifndef ROOTFINDING_H
#define ROOTFINDING_H

#include <cmath>
#include <iomanip>
#include <functional>


namespace RootFinding
{
	// fn is a 1D function
	// between start and end there is at least one zero!
	template <typename Functor, typename T>
	T BinarySearch(Functor fn, T start = 0., T end = 1., T err = 1e-6)
	{
		// resolution too small
		if (std::abs(start - end) <= err)
			return (start + end) / 2.;

		if (start > end) // basic check
			std::swap(start, end);

		T fS = fn(start);
		T fE = fn(end);

		if (fS == 0)
			return start;

		if (fE == 0)
			return end;

		if (fS * fE > 0) 
			return std::nan("");

		// else fS * fE < 0
		while (std::abs(start - end) > err) {
			T fM = fn((start + end) / 2.);

			if (fS * fM < 0) {
				end = (start + end) / 2.;
				fE = fM;
			}
			else if (fE * fM < 0) {
				start = (start + end) / 2.;
				fS = fM;
			}
			else // (fM == 0)
				break;
		}

		return (start + end) / 2.;
	}

	// determine if there is a solution in the interval
	// return end if interval contains at least one solution
	// return value between (start, end) if interval contains at least two solutions
	// return start if interval is not good
	template<typename Functor, typename T>
	T CheckInterval(Functor fn, T start = 0., T end = 1., T err = 1e-6, size_t max_depth = 8)
	{
		// too small
		if (std::abs(start - end) <= err) {
			//std::cout << "end\n";
			return start;
		}

		T fS = fn(start);
		T fE = fn(end);

		if (fS * fE < 0)
			return end;

		// divide interval
		for (std::size_t depth = 0; depth < max_depth; ++depth) {
			T split = (end - start) / (1UL << depth);
			for (std::size_t n = 1UL; n < (1UL << depth); n += 2UL) {
				T fM = fn(start + n * split);
				if (fS * fM < 0)
					return start + n * split;
			}
		}

		return start;
	}
}

#endif

/*
			//return Bisect(fn, start, (start+end) / 2., err);
		// not sure if between start and end there is solution, so split
		cc += 2;
		if (fn(start) * fn(end) > 0.) {
			//std::cout << "split\n";
			auto low = Zeros(fn, cc, start, (start + end)/2., err);
			auto hig = Zeros(fn, cc, (start + end)/2., end, err);
			res.insert(res.end(), std::make_move_iterator(low.begin()),
					      std::make_move_iterator(low.end()));
			res.insert(res.end(), std::make_move_iterator(hig.begin()),
					      std::make_move_iterator(hig.end()));
			return res;
		}

		double mid = Bisect(fn, cc, start, end, err);
		// solution is sufficiently far from start and end
		if (std::abs(mid - start) > err && std::abs(mid - end) > err) {
			res.push_back(mid);
			//return res;
		}

		if (std::abs(start - mid) > err) {
			std::vector<double> low = Zeros(fn, cc, start, mid-err, err);
			res.insert(res.end(), std::make_move_iterator(low.begin()),
					std::make_move_iterator(low.end()));
		}
		if (std::abs(end - mid) > err) {
			std::vector<double> hig = Zeros(fn, cc, mid+err, end, err);
			res.insert(res.end(), std::make_move_iterator(hig.begin()),
					      std::make_move_iterator(hig.end()));
		}

		return res;
*/
