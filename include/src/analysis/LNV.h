#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <cmath>

#include "src/tools/Sort.h"
#include <omp.h>

double Poisson(int n, double s)
{
	if (n > 0)
	{
		if (s < 300)	//Poisson
		{
			double p = s/n;
			double ret = 1;
			for (; n > 0; --n)
				ret *= exp(-p) * s / n;

			return ret;
		}
		else	//normal appoximation
			return exp(- pow(n-s, 2) / (2 * s)) /
				sqrt(2 * Const::pi * s);
	}
	else
		return exp(-s) * pow(s, n);
}

double Ratio(int n, double s, double b)
{
	if (s < 300) //poisson
		return exp(b-s) * pow(s / b, n);
	else	//normal appoximation
		return sqrt(b / s) *
			exp(-0.5 * ( pow(n - s, 2) / s - pow(n - b, 2) / b) );
}

double Belt(int b, double CL, int &nA, int &nB)
{
	double s = sqrt(2.5*b);
	nA = 0;
	while (nA < b+1)
	{
		std::vector<double> vp0, vrB;
		std::vector<int> vN, vI;

		s += 0.01;

		int lim = b + 2*s + 1;
		//std::cout << "lim " << std::max(lim, 20) << std::endl;
		for (int n = b, j = 0; n < std::max(lim, 20); ++n, ++j)
		{
			double p0 = Poisson(n, s+b);
			double rB = Ratio(n, s+b, std::max(0, int(n-b)) + b);

			vp0.push_back(p0);
			vrB.push_back(rB);
			vN.push_back(n);
			vI.push_back(j);

			//std::cout << "n " << n << "\ts " << s << "\tu " << std::max(0,int(n-b))
			//	  << "\tp0 " << p0 << "\tR " << rB << std::endl;
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
				//std::cout << vN[j] << "\t";
			}
			else
				break;
		}
		//std::cout << std::endl;
		//std::cout << "signal " << s << "\tsum " << sum << "\tA " << nA << "\tB " << nB << "\tfrom " << 0 << "\tto " << b+lim << std::endl;
		//std::cout << s << "\tA " << nA << "\tB " << nB << std::endl;
	}

	return s;
}

double Belt(int b, double CL)
{
	int nA, nB;
	return Belt(b, CL, nA, nB);
}

int Limits(double x, int &min_x, int &max_x)
{
	min_x = std::max( 0.0, std::floor(x - 3*sqrt(x) - 15));
	max_x = std::max(50.0, std::ceil (x + 3*sqrt(x) + 15));

	return max_x - min_x;
}

int Limits(double x)
{
	int min_x, max_x;
	return Limits(x, min_x, max_x);
}

double ConfidenceRegion_sloppy(double CL,
		      int &imax, int &size_np, int &size_nm, 
		      double sig_np, double sig_nm, 
		      double bak_np = 0, double bak_nm = 0) 
{
	int min_np, max_np;
	int min_nm, max_nm;

	size_np = Limits(sig_np, min_np, max_np);
	size_nm = Limits(sig_np, min_nm, max_nm);

	std::vector<int> idx(size_np * size_nm);
	std::vector<int> vNp(size_np * size_nm), vNm(size_np * size_nm);
	std::vector<double> vProb0(size_np * size_nm), vRatio(size_np * size_nm);

#pragma omp parallel for
		for (int t = 0; t < size_np * size_nm; ++t)
		{
			int np = min_np + t / size_nm;
			int nm = min_nm + t % size_nm;
			//std::cout << "LNC : (" << np << " , " << nm << ")" << std::endl;

			if (pow(double(np - sig_np) / size_np, 2) +
			    pow(double(nm - sig_nm) / size_nm, 2) > 0.25)
				continue;

			double p0 = Poisson(np, sig_np + bak_np) * Poisson(nm, sig_nm + bak_nm);
			double rB = Ratio(np, sig_np + bak_np,
						std::max(int(np - bak_np), 0) + bak_np)
				  * Ratio(nm, sig_nm + bak_nm,
						std::max(int(nm - bak_nm), 0) + bak_nm);

			vNp[t]	  = np;
			vNm[t]    = nm;
			vProb0[t] = p0;
			vRatio[t] = rB;
			idx[t]    = t;
		}

	std::sort(idx.begin(), idx.end(), Sort(vRatio, Sort::descending));

	double sum = 0;

	double rad_min_np = vNp[idx[0]];
	double rad_max_np = vNp[idx[0]];
	double rad_min_nm = vNm[idx[0]];
	double rad_max_nm = vNm[idx[0]];

	imax = 1;
	for ( ; imax < idx.size(); ++imax)
	{
		int ii = idx[imax];

		if (vNp[ii] < rad_min_np)
			rad_min_np = vNp[ii];
		if (vNp[ii] > rad_max_np)
			rad_max_np = vNp[ii];

		if (vNm[ii] < rad_min_nm)
			rad_min_nm = vNm[ii];
		if (vNm[ii] > rad_max_nm)
			rad_max_nm = vNm[ii];


		if (sum < CL)
			sum += vProb0[ii];
		else
			break;
	}

	++imax;

	size_np = rad_max_np - rad_min_np;
	size_nm = rad_max_nm - rad_min_nm;

	return sum;
}

double ConfidenceRegion(double CL,
			std::vector<int> &vNp_ord, std::vector<int> &vNm_ord,
			double sig_np, double sig_nm, 
			double bak_np = 0, double bak_nm = 0) 
{
	int min_np, max_np;
	int min_nm, max_nm;

	int size_np = Limits(sig_np + bak_np, min_np, max_np);
	int size_nm = Limits(sig_nm + bak_nm, min_nm, max_nm);

	std::vector<int> idx, vNp, vNm;
	std::vector<double> vProb0, vRatio;
	idx.reserve(size_np * size_nm);
	vNp.reserve(size_np * size_nm);
	vNm.reserve(size_np * size_nm);
	vProb0.reserve(size_np * size_nm);
	vRatio.reserve(size_np * size_nm);

	//std::cout << "\tcalculating probabilities" << std::endl;
//#pragma omp parallel for
	for (int t = 0; t < size_np * size_nm; ++t)
	{
		int np = min_np + t / size_nm;
		int nm = min_nm + t % size_nm;

		if (pow(double(np - sig_np - bak_np) / size_np, 2) +
		    pow(double(nm - sig_nm - bak_nm) / size_nm, 2) > 0.25)
			continue;

		double p0 = Poisson(np, sig_np + bak_np) *
			Poisson(nm, sig_nm + bak_nm);
		double rB = Ratio(np, sig_np + bak_np,
				std::max(np - bak_np, 0.0) + bak_np) *
			Ratio(nm, sig_nm + bak_nm,
					std::max(nm - bak_nm, 0.0) + bak_nm);

		idx.push_back(idx.size());
		vNp.push_back(np);
		vNm.push_back(nm);
		vProb0.push_back(p0);
		vRatio.push_back(rB);
	}

	//std::cout << "\tsorting" << std::endl;
	//std::sort(idx.begin(), idx.end(), Sort(vRatio, Sort::descending));

	double sum = 0;
	vNp_ord.clear();
	vNm_ord.clear();
	//for (int i; i < idx.size(); ++i)
	//{
	//	int ii = idx[i];

	//	vNp_ord.push_back(vNp[ii]);
	//	vNm_ord.push_back(vNm[ii]);

	//	if (sum < CL)
	//		sum += vProb0[ii];
	//	else
	//		break;
	//}

	while (sum < CL)
	{
		std::vector<int>::iterator ii = std::max_element(idx.begin(), idx.end(),
								 Sort(vRatio, Sort::ascending));
		//std::vector<double>::iterator ir = std::max_element(vRatio.begin(), vRatio.end());
		if (ii != idx.end())
		{
			int n = *ii;

			sum += vProb0[n];
			vNp_ord.push_back(vNp[n]);
			vNm_ord.push_back(vNm[n]);

			idx.erase(ii);
		}
		else
			break;
	}

	//std::cout << "\tsizes " << vNp_ord.size() << ", " << vNm_ord.size() << std::endl;

	if (sum < CL)
		std::cerr << "WARNING: CL reached is " << sum << std::endl;

	return sum;
}

bool IsSensitiveToLNV_sloppy(double CL,
		      int &size_LNC_np, int &size_LNC_nm, int &size_LNV_np, int &size_LNV_nm,
		      double s_FHC, double s_RHC, double r_FHC, double r_RHC,
		      double sbFHC = 0, double sbRHC = 0, double rbFHC = 0, double rbRHC = 0)
{
	//std::ofstream all("all.dat", std::ios_base::app);
	std::ofstream pla("pla.dat", std::ios_base::app);
	std::ofstream radlnc("radlnc.dat", std::ios_base::app);
	std::ofstream radlnv("radlnv.dat", std::ios_base::app);
//	std::cout << "Checking..." << std::endl;
	int min_LNC_np, max_LNC_np, min_LNC_nm, max_LNC_nm;
	int min_LNV_np, max_LNV_np, min_LNV_nm, max_LNV_nm;

	size_LNC_np = Limits(s_FHC, min_LNC_np, max_LNC_np);
	size_LNC_nm = Limits(s_RHC, min_LNC_nm, max_LNC_nm);
	size_LNV_np = Limits(r_FHC, min_LNV_np, max_LNV_np);
	size_LNV_nm = Limits(r_RHC, min_LNV_nm, max_LNV_nm);


	std::cout << s_FHC << "\tLNC : " <<  min_LNC_np << "\t" << max_LNC_np << "\t" << min_LNC_nm << "\t" << max_LNC_nm << std::endl;
	std::cout << r_FHC << "\tLNV : " <<  min_LNV_np << "\t" << max_LNV_np << "\t" << min_LNV_nm << "\t" << max_LNV_nm << std::endl;

	std::vector<int> idx_LNC(size_LNC_np*size_LNC_nm), idx_LNV(size_LNV_np*size_LNV_nm);
	std::vector<int> vPos_LNC(size_LNC_np*size_LNC_nm), vNeg_LNC(size_LNC_np*size_LNC_nm);
	std::vector<int> vPos_LNV(size_LNV_np*size_LNV_nm), vNeg_LNV(size_LNV_np*size_LNV_nm);
	std::vector<double> vProb0_LNC(size_LNC_np*size_LNC_nm), vRatio_LNC(size_LNC_np*size_LNC_nm);
	std::vector<double> vProb0_LNV(size_LNV_np*size_LNV_nm), vRatio_LNV(size_LNV_np*size_LNV_nm);

	int idP_LNC = -1, idP_LNV = -1;
	int idM_LNC = -1, idM_LNV = -1;
	double pMax_LNC = -1, pMax_LNV = -1;

#pragma omp parallel //for default(none) shared(std::cout, size, nFHC, nRHC, nLNV, idx_LNC, idx_LNV, vPos_LNC, vPos_LNV, vNeg_LNC, vNeg_LNV, vProb0_LNC, vProb0_LNV, vRatio_LNV, vRatio_LNC)
	{
#pragma omp for
		for (int t = 0; t < size_LNC_np * size_LNC_nm; ++t)
			//for (int np = min_LNC_np; np < max_LNC_np; ++np)
			//	for (int nm = min_LNC_nm; nm < max_LNC_nm; ++nm)
		{
			int np = min_LNC_np + t / size_LNC_nm;
			int nm = min_LNC_nm + t % size_LNC_nm;
			//std::cout << "LNC : (" << np << " , " << nm << ")" << std::endl;

			double dist = pow(double(np - s_FHC) / size_LNC_np, 2) +
				      pow(double(nm - s_RHC) / size_LNC_nm, 2); 
			if (dist > 0.25)
				continue;

			double p0 = Poisson(np, s_FHC + sbFHC) * Poisson(nm, s_RHC + sbRHC);
			double rB = Ratio(np, s_FHC + sbFHC,
						std::max(int(np - sbFHC), 0) + sbFHC)
				  * Ratio(nm, s_RHC + sbRHC,
						std::max(int(nm - sbRHC), 0) + sbRHC);

			vPos_LNC[t]   = np;
			vNeg_LNC[t]   = nm;
			vProb0_LNC[t] = p0;
			vRatio_LNC[t] = rB;
			idx_LNC[t]    = t;

			//all << nFHC << "\t" << nRHC << "\t" << r_FHC << "\t" << r_RHC << "\t" << np << "\t" << nm << "\t" << p0 << "\t" << rB << std::endl;
			if (pMax_LNC < p0)
			{
				pMax_LNC = p0;
				idP_LNC = np;
				idM_LNC = nm;
			}
		}
		//all << "\n\n" << std::endl;
#pragma omp for
		for (int t = 0; t < size_LNV_np * size_LNV_nm; ++t)
			//for (int np = min_LNV_np; np < max_LNV_np; ++np)
			//	for (int nm = min_LNV_nm; nm < max_LNV_nmize; ++nm)
		{
			int np = min_LNV_np + t / size_LNV_nm;
			int nm = min_LNV_nm + t % size_LNV_nm;
			//std::cout << "LNV : (" << np << " , " << nm << ")" << std::endl;

			double dist = pow(double(np - r_FHC) / size_LNV_np, 2) +
				      pow(double(nm - r_RHC) / size_LNV_nm, 2); 
			if (dist > 0.25)
				continue;

			double p0 = Poisson(np, r_FHC + rbFHC) * Poisson(nm, r_RHC + rbRHC);
			double rB = Ratio(np, r_FHC + rbFHC,
					      std::max(int(np - rbFHC), 0) + rbFHC)
				  * Ratio(nm, r_RHC + rbRHC,
					      std::max(int(nm - rbRHC), 0) + rbRHC);

			vPos_LNV[t]   = np;
			vNeg_LNV[t]   = nm;
			vProb0_LNV[t] = p0;
			vRatio_LNV[t] = rB;
			idx_LNV[t]    = t;

			//all << nFHC << "\t" << nRHC << "\t" << vFHC << "\t" << vRHC << "\t" << np << "\t" << nm << "\t" << p0 << "\t" << rB << std::endl;
			if (pMax_LNV < p0)
			{
				pMax_LNV = p0;
				idP_LNV = np;
				idM_LNV = nm;
			}
		}
		//all << "\n\n" << std::endl;
	}

#pragma omp parallel sections
	{
	#pragma omp section
		std::sort(idx_LNC.begin(), idx_LNC.end(), Sort(vRatio_LNC, Sort::descending));

	#pragma omp section
		std::sort(idx_LNV.begin(), idx_LNV.end(), Sort(vRatio_LNV, Sort::descending));
	}

	//std::cout << "idx LNC : " << vPos_LNC[idx_LNC.front()] << " == "
	//			  << vPos_LNC[idP_LNC] << ", "
	//			  << vNeg_LNC[idx_LNC.front()] << " == "
	//			  << vNeg_LNC[idM_LNC]<< std::endl;
	//std::cout << "idx LNV : " << vPos_LNV[idx_LNV.front()] << " == "
	//			  << vPos_LNV[idP_LNV] << ", "
	//			  << vNeg_LNV[idx_LNV.front()] << " == "
	//			  << vNeg_LNV[idM_LNV]<< std::endl;

	//std::cout << "size LNC: " << size_LNC_np * size_LNC_nm * pMax_LNC << " > "
	//			  << CL * 18.0 / Const::pi << " ? " << std::boolalpha
	//			  << (size_LNC_np * size_LNC_nm * pMax_LNC > CL * 18.0 / Const::pi) << std::endl;
	//std::cout << "size LNV: " << size_LNV_np * size_LNV_nm * pMax_LNV << " > "
	//			  << CL * 18.0 / Const::pi << " ? " << std::boolalpha
	//			  << (size_LNV_np * size_LNV_nm * pMax_LNV > CL * 18.0 / Const::pi) << std::endl;

	//std::cout << "#dirac area\n";
	double sum_LNC = 0;
	int imax = 0;

	double rad_min_Pos_LNC = 1e5;
	double rad_max_Pos_LNC = -1;
	double rad_min_Neg_LNC = 1e5;
	double rad_max_Neg_LNC = -1;

	for (int i = 0; i < idx_LNC.size(); ++i)
	{
		int ii = idx_LNC[i];

		pla << s_FHC << "\t" << s_RHC << "\t" << r_FHC << "\t" << r_RHC << "\t" << vPos_LNC[ii] << "\t" << vNeg_LNC[ii] << "\t" << sum_LNC << std::endl;

		if (vPos_LNC[ii] < rad_min_Pos_LNC)
			rad_min_Pos_LNC = vPos_LNC[ii];
		if (vPos_LNC[ii] > rad_max_Pos_LNC)
			rad_max_Pos_LNC = vPos_LNC[ii];

		if (vNeg_LNC[ii] < rad_min_Neg_LNC)
			rad_min_Neg_LNC = vNeg_LNC[ii];
		if (vNeg_LNC[ii] > rad_max_Neg_LNC)
			rad_max_Neg_LNC = vNeg_LNC[ii];


		if (sum_LNC < CL)
		{
			sum_LNC += vProb0_LNC[ii];
			imax = i+1;
		}
		else
			break;
	}
	radlnc << s_FHC << "\t" << s_RHC << "\t" << rad_min_Pos_LNC << "\t" << rad_max_Pos_LNC << "\t" << rad_min_Neg_LNC << "\t" << rad_max_Neg_LNC << std::endl;

	pla << "\n\n" << std::endl;
	//std::cout << "#major area\n";
	double sum_LNV = 0;
	int jmax = 0;

	double rad_min_Pos_LNV = 1e5;
	double rad_max_Pos_LNV = -1;
	double rad_min_Neg_LNV = 1e5;
	double rad_max_Neg_LNV = -1;

	for (int j = 0; j < idx_LNV.size(); ++j)
	{
		int jj = idx_LNV[j];

		pla << s_FHC << "\t" << s_RHC << "\t" << r_FHC << "\t" << r_RHC << "\t" << vPos_LNV[jj] << "\t" << vNeg_LNV[jj] << "\t" << sum_LNV << std::endl;

		if (vPos_LNV[jj] < rad_min_Pos_LNV)
			rad_min_Pos_LNV = vPos_LNV[jj];
		if (vPos_LNV[jj] > rad_max_Pos_LNV)
			rad_max_Pos_LNV = vPos_LNV[jj];

		if (vNeg_LNV[jj] < rad_min_Neg_LNV)
			rad_min_Neg_LNV = vNeg_LNV[jj];
		if (vNeg_LNV[jj] > rad_max_Neg_LNV)
			rad_max_Neg_LNV = vNeg_LNV[jj];
		if (sum_LNV < CL)
		{
			sum_LNV += vProb0_LNV[jj];
			jmax = j+1;
		}
		else
			break;
	}
	radlnv << r_FHC << "\t" << r_RHC << "\t" << rad_min_Pos_LNV << "\t" << rad_max_Pos_LNV << "\t" << rad_min_Neg_LNV << "\t" << rad_max_Neg_LNV << std::endl;

	pla << "\n\n" << std::endl;

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

bool IsSensitiveToLNV(double CL,
		      std::vector<int> &dp_dir, std::vector<int> &dm_dir,
		      std::vector<int> &dp_maj, std::vector<int> &dm_maj,
		      double sig_np_dir, double sig_nm_dir,
		      double sig_np_maj, double sig_nm_maj,
		      double bak_np_dir = 0, double bak_nm_dir = 0,
		      double bak_np_maj = 0, double bak_nm_maj = 0)
{
	int size_np_dir = Limits(sig_np_dir + bak_np_dir);
	int size_nm_dir = Limits(sig_nm_dir + bak_nm_dir);
	int size_np_maj = Limits(sig_np_maj + bak_np_maj);
	int size_nm_maj = Limits(sig_nm_maj + bak_nm_maj);

	//std::cout << "distance " << pow(sig_np_dir + bak_np_dir - sig_np_maj - bak_np_maj, 2) + pow(sig_nm_dir + bak_nm_dir - sig_nm_maj - bak_nm_maj, 2)
	//	  << " vs " << pow(std::max(size_np_dir, size_nm_dir), 2) / 4.0 + pow(std::max(size_np_maj, size_nm_maj), 2) / 4.0
	//	  << std::endl;

	//regions won't overlap if true
	//because distance is greater than covered regions
	if ( pow(sig_np_dir + bak_np_dir - sig_np_maj - bak_np_maj, 2) +
	     pow(sig_nm_dir + bak_nm_dir - sig_nm_maj - bak_nm_maj, 2) >
	     pow(std::max(size_np_dir, size_nm_dir), 2) / 4.0 +
	     pow(std::max(size_np_maj, size_nm_maj), 2) / 4.0 )
		return true;

	std::vector<int> vNp_dir, vNm_dir;
	std::vector<int> vNp_maj, vNm_maj;
	double sum_dir, sum_maj;
	//std::cout << "defining regions" << std::endl;
//#pragma omp parallel sections //shared(vN
	{
//#pragma omp section
	sum_dir = ConfidenceRegion(CL, vNp_dir, vNm_dir,
			sig_np_dir, sig_nm_dir,
			bak_np_dir, bak_nm_dir);
//#pragma omp section
	sum_maj = ConfidenceRegion(CL, vNp_maj, vNm_maj,
			sig_np_maj, sig_nm_maj,
			bak_np_maj, bak_nm_maj);
	}

	//std::cout << "sizes " << vNp_dir.size() << ", " << vNm_dir.size() << std::endl;
	//std::cout << "sizes " << vNp_maj.size() << ", " << vNm_maj.size() << std::endl;

	//std::cout << "comparing sums" << std::endl;
	if (sum_dir < CL || sum_maj < CL)
		return false;

	//std::cout << "comparing regions" << std::endl;
	bool overlap = false;
	for (int i = 0; i < vNp_dir.size(); ++i)
	{
		for (int j = 0; j < vNp_maj.size(); ++j)
		{
			if (vNp_dir[i] == vNp_maj[j] &&
			    vNm_dir[i] == vNm_maj[j])
			{
				overlap = true;
				break;
			}
		}

		if (overlap)
			break;
	}

	dp_dir.push_back(*std::min_element(vNp_dir.begin(), vNp_dir.end()));
	dp_dir.push_back(*std::max_element(vNp_dir.begin(), vNp_dir.end()));
	dm_dir.push_back(*std::min_element(vNm_dir.begin(), vNm_dir.end()));
	dm_dir.push_back(*std::max_element(vNm_dir.begin(), vNm_dir.end()));

	dp_maj.push_back(*std::min_element(vNp_maj.begin(), vNp_maj.end()));
	dp_maj.push_back(*std::max_element(vNp_maj.begin(), vNp_maj.end()));
	dm_maj.push_back(*std::min_element(vNm_maj.begin(), vNm_maj.end()));
	dm_maj.push_back(*std::max_element(vNm_maj.begin(), vNm_maj.end()));

	return !overlap;
	//returns true  if there is no overlap -> yes sensitivity
	//returns false if there is overlap    -> no  sensitivity
}

bool IsSensitiveToLNV(double CL,
		      double sig_np_dir, double sig_nm_dir,
		      double sig_np_maj, double sig_nm_maj,
		      double bak_np_dir = 0, double bak_nm_dir = 0,
		      double bak_np_maj = 0, double bak_nm_maj = 0)
{
	std::vector<int> dp_dir, dm_dir, dp_maj, dm_maj;
	return IsSensitiveToLNV(CL, dp_dir, dm_dir, dp_maj, dm_maj,
		      sig_np_dir, sig_nm_dir, sig_np_maj, sig_nm_maj,
		      bak_np_dir, bak_nm_dir, bak_np_maj, bak_nm_maj);
}
