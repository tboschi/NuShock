#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>

#include <omp.h>

class decrease
{
	private:
		std::vector<double> vInd;
	public:
		decrease(std::vector<double>& vect) : vInd(vect) {}
		bool operator()(int i, int j) const { return vInd.at(i) > vInd.at(j); }
};

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
		std::sort(idx_LNC.begin(), idx_LNC.end(), decrease(vRatio_LNC));

	#pragma omp section
		std::sort(idx_LNV.begin(), idx_LNV.end(), decrease(vRatio_LNV));
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
