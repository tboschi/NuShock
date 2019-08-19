#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>

#include "src/tools/Sort.h"
#include <omp.h>


double Poisson(int n, double s)
{
	if (n > 0)
	{
		double p = s/n;
		double ret = 1;
		for (; n > 0; --n)
			ret *= exp(-p) * s / n;

		return ret;
	}
	else
		return exp(-s) * pow(s, n);
}

double Ratio(int n, double s, double b)
{
	return exp(b-s) * pow(s / b, n);
}

void feldcous(double nFHC, double nRHC, int np, int nm, double &p0, double &rB)
{
	p0 = Poisson(np, nFHC) * Poisson(nm, nRHC);
	rB = Ratio(np, nFHC, std::max(np, 0))
		* Ratio(nm, nRHC, std::max(nm, 0));
}

double Belt(int b, double CL, int &nA, int &nB)
{
	double s = sqrt(2.5*b);
	nA = 0;
	while (nA < b+1)
	{
		std::vector<double> vp0, vrB;
		std::vector<int> vN, vI;

		s += 0.005;

		double lim = 70*log(b + 270);	//from plot
		for (int n = std::max(0, b-1), j = 0; n < b + lim; ++n, ++j)
		//for (unsigned int n = 140; n < 160; ++n, ++j)
		{
			double p0 = Poisson(n, s+b);
			double rB = Ratio(n, s+b, std::max(0, int(n-b)) + b);

			vp0.push_back(p0);
			vrB.push_back(rB);
			vN.push_back(n);
			vI.push_back(j);

			//std::cout << "n " << n << "\ts " << s << "\tu " << std::max(0,int(n-b));
			//std::cout << "\tp0 " << p0 << "\tPA " << PA;
			//std::cout << "\tR " << Ratio << std::endl;
		}
		std::sort(vI.begin(), vI.end(), Sort(vrB, Sort::descending));

		double sum = vp0[vI.front()];
		nA = nB = vN[vI.front()];
		for (int i = 1; i < vI.size(); ++i)
		{
			if (sum <= CL)
			{
				int j = vI[i];
				sum += vp0[j];

				if (nA > vN[j])
					nA = vN[j];
				if (nB < vN[j])
					nB = vN[j];
			}
			else
				break;
		}
		//std::cout << "signal " << s << "\tsum " << sum << "\tA " << nA << "\tB " << nB << "\tfrom " << std::max(0,b-5) << "\tto " << b+lim+20 << std::endl;
		//std::cout << s << "\tA " << nA << "\tB " << nB << std::endl;
	}

	return s;
}


int Limits(double x, int &min_x, int &max_x)
{
	double a = 10.0;
	double b = 1 + 2.0 / a;
	min_x = std::max(0.0, -a + x / b);
	max_x = std::max(0.0,  a + x * b);

	return max_x - min_x;
}


bool IsSensitiveToLNV(double nFHC, double nRHC, double nLNV, double CL)
{
//	std::cout << "Checking..." << std::endl;
	int min_LNC_np, max_LNC_np, min_LNC_nm, max_LNC_nm;
	int min_LNV_np, max_LNV_np, min_LNV_nm, max_LNV_nm;

	int size_LNC_np = Limits(nFHC, min_LNC_np, max_LNC_np);
	int size_LNC_nm = Limits(nFHC, min_LNC_nm, max_LNC_nm);
	min_LNC_nm = 0;
	int size_LNV_np = Limits(nLNV, min_LNV_np, max_LNV_np);
	int size_LNV_nm = Limits(nLNV, min_LNV_nm, max_LNV_nm);

	std::vector<int> idx_LNC(size_LNC_np*size_LNC_nm), idx_LNV(size_LNV_np*size_LNV_nm);
	std::vector<int> vPos_LNC(size_LNC_np*size_LNC_nm), vNeg_LNC(size_LNC_np*size_LNC_nm);
	std::vector<int> vPos_LNV(size_LNV_np*size_LNV_nm), vNeg_LNV(size_LNV_np*size_LNV_nm);
	std::vector<double> vProb0_LNC(size_LNC_np*size_LNC_nm), vRatio_LNC(size_LNC_np*size_LNC_nm);
	std::vector<double> vProb0_LNV(size_LNV_np*size_LNV_nm), vRatio_LNV(size_LNV_np*size_LNV_nm);

#pragma omp parallel //for default(none) shared(std::cout, size, nFHC, nRHC, nLNV, idx_LNC, idx_LNV, vPos_LNC, vPos_LNV, vNeg_LNC, vNeg_LNV, vProb0_LNC, vProb0_LNV, vRatio_LNV, vRatio_LNC)
	//for (int np = 0; np < size; ++np)
	//{
	//	for (int nm = 0; nm < size; ++nm)
	//	{
	//		double p0, rB;
	//		int t = np*size + nm;
	//		//std::cout << "t " << omp_get_thread_num() << " : " << t << std::endl;

	//		///////// LNC
	//		feldcous(nFHC, nRHC, np, nm, p0, rB);

	//		vPos_LNC[t]   = np;
	//		vNeg_LNC[t]   = nm;
	//		vProb0_LNC[t] = p0;
	//		vRatio_LNC[t] = rB;
	//		idx_LNC[t]    = t;

	//		///////// LNV
	//		feldcous(nLNV, nLNV, np, nm, p0, rB);

	//		vPos_LNV[t]   = np;
	//		vNeg_LNV[t]   = nm;
	//		vProb0_LNV[t] = p0;
	//		vRatio_LNV[t] = rB;
	//		idx_LNV[t]    = t;
	//	}
	//}
	{
#pragma omp for
		for (int t = 0; t < size_LNC_np * size_LNC_nm; ++t)
			//for (int np = min_LNC_np; np < max_LNC_np; ++np)
			//	for (int nm = min_LNC_nm; nm < max_LNC_nm; ++nm)
		{
			int np = min_LNC_np + t / size_LNC_np;
			int nm = min_LNC_nm + t % size_LNC_nm;
			//std::cout << "LNC : (" << np << " , " << nm << ")" << std::endl;

			double p0 = Poisson(np, nFHC) * Poisson(nm, nRHC);
			double rB = Ratio(np, nFHC, std::max(np, 0))
				  * Ratio(nm, nRHC, std::max(nm, 0));

			vPos_LNC[t]   = np;
			vNeg_LNC[t]   = nm;
			vProb0_LNC[t] = p0;
			vRatio_LNC[t] = rB;
			idx_LNC[t]    = t;
		}
#pragma omp for
		for (int t = 0; t < size_LNV_np * size_LNV_nm; ++t)
			//for (int np = min_LNV_np; np < max_LNV_np; ++np)
			//	for (int nm = min_LNV_nm; nm < max_LNV_nmize; ++nm)
		{
			int np = min_LNV_np + t / size_LNV_np;
			int nm = min_LNV_nm + t % size_LNV_nm;
			//std::cout << "LNV : (" << np << " , " << nm << ")" << std::endl;

			double p0 = Poisson(np, nLNV) * Poisson(nm, nLNV);
			double rB = Ratio(np, nLNV, std::max(np, 0))
				  * Ratio(nm, nLNV, std::max(nm, 0));

			vPos_LNV[t]   = np;
			vNeg_LNV[t]   = nm;
			vProb0_LNV[t] = p0;
			vRatio_LNV[t] = rB;
			idx_LNV[t]    = t;
		}
	}

#pragma omp parallel sections
	{
	#pragma omp section
		std::sort(idx_LNC.begin(), idx_LNC.end(), Sort(vRatio_LNC, Sort::descending));

	#pragma omp section
		std::sort(idx_LNV.begin(), idx_LNV.end(), Sort(vRatio_LNV, Sort::descending));
	}

	double sum_LNC = 0;
	int imax = 0;
	for (int i = 0; i < idx_LNC.size(); ++i)
	{
		int ii = idx_LNC[i];

		if (sum_LNC < CL)
		{
			sum_LNC += vProb0_LNC[ii];
			imax = i+1;
		}
		else
			break;
	}

	double sum_LNV = 0;
	int jmax = 0;
	for (int j = 0; j < idx_LNV.size(); ++j)
	{
		int jj = idx_LNV[j];

		if (sum_LNV < CL)
		{
			sum_LNV += vProb0_LNV[jj];
			jmax = j+1;
		}
		else
			break;
	}

	//std::cout << "CL : " << sum_LNC << "\t" << sum_LNV << std::endl;
	if (sum_LNC < CL || sum_LNV < CL)
		return false;

	bool overlap = false;
	for (int i = 0; i < imax; ++i)
	{
		int ii = idx_LNC[i];
		for (int j = 0; j < jmax; ++j)
		{
			int jj = idx_LNV[j];

			if (vPos_LNC[ii] == vPos_LNV[jj] &&
			    vNeg_LNC[ii] == vNeg_LNV[jj])
			{
				//std::cout << "overlap at (" << vPos_LNC[ii] << ", " << vNeg_LNC[ii] << ") [" << i << "]" << std::endl;
				overlap = true;
				break;
			}
		}

		if (overlap)
			break;
	}

	return !overlap;
	//returns true  if there is no overlap -> yes sensitivity
	//returns false if there is overlap    -> no  sensitivity
}
