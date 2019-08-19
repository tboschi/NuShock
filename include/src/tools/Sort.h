#ifndef Sort_H
#define Sort_H

#include <iostream>
#include <fstream>
#include <vector>

class Sort
{
	public:
		enum Type
		{
			descending = 0,
			ascending
		};

		Sort(const std::vector<double> &vv, Type ascent = ascending) :
			vX(vv),
			kAscending(ascent)
		{
		}

		bool operator()(int i, int j) const
		{
			return (2*kAscending - 1) * ( vX[i] - vX[j] ) < 0;
		}

	private:
		std::vector<double> vX;
		Type kAscending;
};

#endif
