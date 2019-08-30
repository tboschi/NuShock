/*
 * Cross section for Nucleus Nucleus scattering to open q qbar production
 * Used for p + A -> c + c_ + X
 * Author: Tommaso Boschi
 */

#ifndef HADRONSCATTER_H
#define HADRONSCATTER_H

#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <cstring>
#include <iomanip>
#include <fstream>

#include "tools.h"
#include "TLorentzVector.h"
#include "TGenPhaseSpace.h"

#include "cuba.h"

#include "LHAPDF/PDF.h"
#include "LHAPDF/AlphaS.h"

class Hadron
{
	public:
		Hadron(std::string ProbeName, std::string TargetName, int OutQuark);

		void SetProb(double Px, double Py, double Pz, double E);
		void SetTarg(double Px, double Py, double Pz, double E);
		void SetParton(double s);
		void SetFromPS();
		double Pt(int e = 1);

		double M(int e = 1);
		double s(int e = 1);
		double s_0(int e = 1);

		double s_(int e = 1);
		double t_(int e = 1);
		double u_(int e = 1);
		double Q2_(int e = 1);
		void SetCMEnergy(double X);
		double CosT_(int e = 1);
		void SetOmega(double X);
		void SetQ2(double X);
		double Q2Lim(double &Min, double &Max);

		double aS(int e = 1.0);

		double dXSdQ2_qq();
		double dXSdQ2_gg();
		double dXSdOmega_qq();
		double dXSdOmega_gg();

		double GetPDFSet(double Q2, double x1, double x2);
		double operator()(const double *x, int nComp);
		//void Integrand(const int *nDim, const double *x, const int *nComp, double f[]);
		double Variable(double *x);

		double Total(double &Error, double &Chi2Prob);
		double Spectrum(double &Error, double &Chi2Prob);
		double XSec(double &Error, double &Chi2Prob);

		bool PSCut();

	private:
		LHAPDF::PDF *ProbPdf, *TargPdf;

		int integrandType;
		double fs_, fQ2_, fcosT_;

		TLorentzVector *Prob, *Targ;
		TLorentzVector *Part1, *Part2, *Quark3, *Quark4;

		TGenPhaseSpace *PartonPS;

		int ValQuark;
		double ValMass;
		std::ofstream *fout;
};

#endif
