#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <cstring>
#include <getopt.h>

#include "tools.h"
//#include "analysis.h"

// n observed events
// s is signal
// b is background
double llratio(double n, double s, double b) {
	// correct also if any term is zero
	return 2 * (s - std::max(n - b, 0.))
	     - 2 * (n > 0 ? n * (log(s + b) - log(std::max(n, b))) : 0);
}
// function above has minimum value for n = s + b

double poisson(double n, double s, double b) {
	s += b;
	if (n > 0) {
		if (true) {	//Poisson
			double p = std::exp(- s / n) * s;
			long double ret = p / n;
			for (--n; n > 0; --n) {
				ret *= p / n;
			}

			return ret;
		}
		else	//normal appoximation
			return exp(- pow(n-s, 2) / (2 * s)) /
				std::sqrt(2 * Const::pi * s);
	}

	return exp(-s);
}

/* poisson is exp(-s) * pow(s, n) / n!
 * CL is reached by sum_n exp(-s) * pow(s, n) / n!
 * exp(-s) * sum[ pow(s, n) / n! ]
 * and n contiguos in a range around n0
 * exp(-s) * sum[ pow(s, n0 + n - n0) / n0! * n0! / n! ]
 * exp(-s) * pow(s, n0) / n0! * sum[ pow(s, n - n0) * n0! / n! ]
*/

// function computes a factorized term in the poissonian sum
double partial(double n0, double n, double s, double b) {
	// n > n0 correct order
	s += b;
	long double ret = 1.;
	for (double k = std::max(n, n0); k > std::min(n, n0); --k)
		ret *= s / k;

	return n >= n0 ? ret : 1./ret;
}


double belt(double bak, double CL, int &nA, int &nB)
{
	double sig = std::sqrt(bak);
	double nL = bak, nR = bak;
	while (nL < bak + 1) {
		sig += 0.001;

		double num = std::round(sig + bak);
		double best = poisson(num, sig, bak);
		double sum = 1.;	 // only sum after the peak was reached
		double cor = best;
		// one direction, likelihood should be unimodal

		nL = num;
		nR = num;
		while (sum < CL / best) {

			// left point is better
			if (nL > 0 && llratio(nL-1, sig, bak) <= llratio(nR+1, sig, bak)) {
				--nL;
				sum += partial(num, nL, sig, bak);
			}
			else {	// right point is better
				++nR;
				sum += partial(num, nR, sig, bak);
			}

			if (nL < bak + 1)
				break;
		}

		//std::cout << "signal " << sig << " and belt " << nA << " - " << nB << "\n" << std::endl;
	}

	nA = nL;
	nB = nR;

	return sig;
}

// expand belt function in multiple dimensions
std::vector<double> confidence_region(double CL, const std::vector<double> &bak)
{

	double sig = std::sqrt(bak);
	double nL = bak, nR = bak;
	while (nL < bak + 1) {
		sig += 0.001;

		double num = std::round(sig + bak);
		double best = poisson(num, sig, bak);
		double sum = 1.;	 // only sum after the peak was reached
		double cor = best;
		// one direction, likelihood should be unimodal

		nL = num;
		nR = num;
		while (sum < CL / best) {

			// left point is better
			if (nL > 0 && llratio(nL-1, sig, bak) <= llratio(nR+1, sig, bak)) {
				--nL;
				sum += partial(num, nL, sig, bak);
			}
			else {	// right point is better
				++nR;
				sum += partial(num, nR, sig, bak);
			}

			if (nL < bak + 1)
				break;
		}

		//std::cout << "signal " << sig << " and belt " << nA << " - " << nB << "\n" << std::endl;
	}

	nA = nL;
	nB = nR;

	return sig;

int main(int argc, char** argv)
{
	const struct option longopts[] = 
	{
		{"background", 	required_argument, 	0, 'b'},
		{"confidence", 	required_argument, 	0, 'C'},
		{"help", 	no_argument,	 	0, 'h'},
		{0,	0, 	0,	0},
	};

	int index; 
	int iarg = 0;
	opterr = 1;

	int b = 0;
	double CL = 0.90;	//90% C.L.

	while((iarg = getopt_long(argc,argv, "b:C:h", longopts, &index)) != -1)
	{
		switch(iarg)
		{
			case 'b':
				b = strtol(optarg, NULL, 10);
				break;
			case 'C':
				CL = strtod(optarg, NULL);
				CL /= 100.0;
				break;
			case 'h':
				return 1;
			default:
				break;
		}
	}

	int nA, nB;
	double s = belt(b, CL, nA, nB);
	
	std::cout << "For " << CL*100.0 << " %, mean signal is " << s << " from " << nA << " to " << nB << std::endl;

	return 0;
}
