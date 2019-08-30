#ifndef BinarySearch_H
#define BinarySearch_H

#include <iostream>
#include <cmath>
#include <algorithm>

template<class T>
double bisect(const T& obj, double start = 0, double end = 1, double err = 1e-6)
{
	if (start > end)
	{
		double tmp = end;
		end = start;
		start = tmp;
	}

	double fS = obj(start);
	double fE = obj(end);

	if (fS * fE > 0)
	{
		//std::cout << "No zero point between " << start << " and " << end << std::endl;
		if (fS > 0 && fE > 0)
			return end;
		else 
			return start;
	}
	else if (fS == 0)
		return start;
	else if (fE == 0)
		return end;

	while (std::abs(start - end) > err)
	{
		//std::cout << "Ranging " << S << " (" << fS << ") - "
		//			<< E << " (" << fE << ")" << std::endl;
		double mid = (start + end) / 2.0;
		double fM = obj(mid);

		if (fM == 0)
		{
			start = mid;
			end = mid;
		}
		else if (fS * fM < 0)
		{
			end = mid;
			fE = fM;
		}
		else if (fE * fM < 0)
		{
			start = mid;
			fS = fM;
		}
	}

	return (start + end) / 2.0;
}

template<class T>
double findInterval(const T& obj, double start = 0, double end = 1, int depth = 8)
{
	double fS = obj(start);
	double fE = obj(end);

	for (int i = 0; i < depth; ++i)	//max 2^8 = 256 intervals
	{
		double split = (end - start) / pow (2, i+1);
		for (int n = 1; n < pow(2, i+1); n += 2)
		{
			double fM = obj(start + n*split);

			if (fS * fM < 0)
				return start + n * split;
		}
	}

	return start + end;
}

#endif
