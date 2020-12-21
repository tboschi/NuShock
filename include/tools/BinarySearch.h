#ifndef BinarySearch_H
#define BinarySearch_H

#include <iostream>

// solve func(x) = c for x
// returns all solutions

namespace BinarySearch {

	// return best in the direction of start
	template<typename F, typename T>
	T search(F func, T c, T start = 0., T end = 1., T err = 1e-6)
	{
		if (start > end) // basic check
			std::swap(start, end);

		T fS = func(start) - c;
		T fE = func(end) - c;

		if (fS * fE > 0)
			return search(func, c, start, (start+end) / 2., err);

		if (fS == 0)
			return start;

		if (fE == 0)
			return end;

		// else fS * fE < 0
		while (std::abs(start - end) > err) {
			T fM = func((start + end) / 2.) - c;

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

		return (start + end) / 2.0;
	}

	template<typename F, typename T>
	std::vector<T> solve(F func, T c, T start = 0., T end = 1., T err = 1e-6)
	{
		std::vector<T> res;
		if (std::abs(start - end) > err) {
			T mid = search(func, c, start, end, err);
			std::vector<T> low = solve(func, c, start, mid, err);
			std::vector<T> hig = solve(func, c, mid, end, err);

			res.insert(res.end(), std::make_move_iterator(low.begin()),
					      std::make_move_iterator(low.end()));
			res.insert(res.end(), std::make_move_iterator(hig.begin()),
					      std::make_move_iterator(hig.end()));
		}
		return res;
	}
}

#endif
