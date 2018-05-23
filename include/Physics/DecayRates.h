/*
 * Decay rate calculator for a sterile neutrino below 2GeV in a minimal model
 * Author: Tommaso Boschi
 */

#ifndef DECAYRATES_H
#define DECAYRATES_H

#include <iostream>
#include <fstream>
#include <vector>
#include <cstring>

//ROOT include
#include "TLorentzVector.h"
#include "TGenPhaseSpace.h"

#include "Tools/Tools.h"
#include "Decay/ThreeBody.h"


class DecayRates
{
	public:
		enum Channel
		{
			_ALL,
			//unclassified
			_nnn,		//3 body	N -> 3 nu
			_nGAMMA,	//2 body	N -> nu photon
			//pure leptonic
			_nEE,		//3 body	N -> nu e e
			_nEMU,		//3 body	N -> nu e mu (via U_e)
			_nMUE,		//3 body	N -> nu mu e (via U_m)
			_nMUMU,		//3 body	N -> nu mu mu
			_nET,		//3 body	N -> nu e tau (via U_e)
			_nTE,		//3 body	N -> nu tau e (via U_t)
			_nMUT,		//3 body	N -> nu mu tau (via U_m)
			_nTMU,		//3 body	N -> nu tau mu (via U_t)
			//pion
			_nPI0,		//2 body	N -> nu pi0
			_EPI,		//2 body	N -> e pi
			_MUPI,		//2 body	N -> mu pi
			_TPI,		//2 body	N -> tau pi
			//kaon
			_EKA,		//2 body	N -> e K
			_MUKA,		//2 body	N -> mu K
			//rho decay 100% in pions
			_nRHO0,		//2 body	N -> rho0
			_ERHO,		//2 body	N -> e rho
			_MURHO,		//2 body	M -> mu rho
			//kaon*
			_EKAx,		//2 body	N -> e K*
			_MUKAx,		//2 body	N -> mu K*
			//other (eta, phi, omega.. )
			_nOMEGA,	//2 body	N -> nu w
			_nETA,		//2 body	N -> nu eta
			_nETAi,		//2 body	N -> nu eta'
			_nPHI,		//2 body	N -> nu phi
			//charm
			_ECHARM,	//2 body	N -> e D+
			//Channels for experimental comparison (EPI, MUPI, nEE, nMUE, nMUMU)
			_ExpALL,	//
		};

		DecayRates(double MSterile = 0.0, double Ue = 0.0, double Um = 0.0, double Ut = 0.0);

		double GetMass();
		double GetUe();
		double GetUm();
		double GetUt();

		void SetMass(double X);
		void SetUe(double X);
		void SetUm(double X);
		void SetUt(double X);

	protected:
		std::vector<double> I_var;

		std::vector<double> vMass;

		//Masses
		const double M_Neutrino;
		const double M_Photon;
		const double M_Electron;
		const double M_Muon;
		const double M_Pion;
		const double M_Pion0;
		const double M_Kaon;
		const double M_Kaon0;
		const double M_Eta;
		const double M_Rho;
		const double M_Rho0;
		const double M_Omega;
		const double M_Kaonx;
		const double M_Kaon0x;
		const double M_Etai;
		const double M_Phi;
		const double M_Tau;
		const double M_Charm;

	private:
		double M_Sterile, M_Sterile_prev;
		double Ue, Um, Ut;

};

#endif
