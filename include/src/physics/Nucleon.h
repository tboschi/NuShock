/*
 * Cross section for nu - N scattering
 * Author: Tommaso Boschi
 */

#ifndef NUCLEON_H
#define NUCLEON_H

#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <cstring>
#include <iomanip>

#include "Tools.h"
#include "TLorentzVector.h"
#include "TGenPhaseSpace.h"
#include "cuba.h"

class Nucleon
{
	public:
		Nucleon(bool Neutrino, bool Nucleon);

		void SetPS(double M_Sterile);
		double GeneratePS();

		double dSigmadQ2(int Neut = 0);
		double dSigmadOmega(int Neut = 0);
		double SigmaTot(int Neut = 0);
		double Amp2(int Neut = 0);
		
		double Q2Lim(double &Min, double &Max);
		double Variable(double dt);
		double Integrate(double (Nucleon::*FF)(int), double A, double B, int Neut = 0);
		void Integrand(const int *nDim, const double *x, const int *nComp, double f[]);
		//static int Integrand(const int *nDim, const double x[], const int *nComp, double f[], void* UserData);
		//double VegasIntegration(double &Error, double &Chi2Prob);

		double A();
		double B();
		double C();

		double F1(int e = 1);
		double F2(int e = 1);
		double GA(int e = 1);
		double F1EM(bool N, int e = 1);
		double F2EM(bool N, int e = 1);
		double GAEM(int e = 1);
		double F1S(int e = 1);
		double F2S(int e = 1);
		double GAS(int e = 1);
		double GEsachs(bool N, int e = 1);
		double GMsachs(bool N, int e = 1);
		double GDipole();

		double s();
		double s_();
		double t();
		double t_();
		double u();
		double u_();
		double Q2();

		int Tau3();
		int Sign();
		double GetHeavyE();

		void SetQ2(double X);

		void SetNeutrino(bool B);
		void SetNeutrino(bool B, TLorentzVector &nu);
		void SetNucleon(bool B);

		void SetProbe(TLorentzVector &Probe);
		void SetProbeM(double dM);
		void SetProbeE(double dE);
		void SetTarget(TLorentzVector &Probe);
		void SetTargetM(double dM);
		void SetTargetE(double dE);
		void SetSterile(TLorentzVector &Probe);
		void SetSterileM(double dM);
		void SetSterileE(double dE);
		void SetRecoil(TLorentzVector &Probe);
		void SetRecoilM(double dM);
		void SetRecoilE(double dE);

		double Acos(double X);

		bool IsEnergyConserved();
		bool IsChanged();
		bool IsAllowed();
		void ResetFormFactors();

	private:
		TLorentzVector *Pl, *Pi, *Ph, *Pf;	//Pl + Pi = Ph + Pf
		TGenPhaseSpace *Event;

		double El, Mn, Ml, Mh;
		double fA, fB, fC;
		double fF1, fF2, fGA, fF1EM, fF2EM, fGAEM, fF1S, fF2S, fGAS, fGEsachs, fGMsachs, fGDipole;

		double fSigmaTot, fQ2Min, fQ2Max;
		double Mh_prev, El_prev;

		const double M_Proton;
		const double M_Neutron;

		int t3, Nu;
};

#endif
