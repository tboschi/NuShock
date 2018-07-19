/*
 * Decay rate calculator for a sterile neutrino below 2GeV in a minimal model
 * Author: Tommaso Boschi
 */

#ifndef AMPLITUDE_H
#define AMPLITUDE_H

#include <iostream>
#include <fstream>
#include <vector>
#include <cstring>

#include "Tools.h"
#include "cuba.h"
#include "Tools/asa047.hpp"

class Amplitude
{
	public:
		enum Process
		{
			DecayRates,
			Production,
			Undefined,
		};

		enum Channel
		{
			_undefined = 0,

			_ALL,

			//decay modes
			//unclassified
			_nnn,		//3 body	N -> 3 nu
			_nGAMMA,	//2 body	N -> nu photon
			//pure leptonic
			_nEE,		//3 body	N -> nu e e
			_nEM,		//3 body	N -> nu e mu (via U_e)
			_nME,		//3 body	N -> nu mu e (via U_m)
			_nMM,		//3 body	N -> nu mu mu
			_nET,		//3 body	N -> nu e tau (via U_e)
			_nTE,		//3 body	N -> nu tau e (via U_t)
			_nMT,		//3 body	N -> nu mu tau (via U_m)
			_nTM,		//3 body	N -> nu tau mu (via U_t)
			//pion
			_nPI0,		//2 body	N -> nu pi0
			_EPI,		//2 body	N -> e pi
			_MPI,		//2 body	N -> mu pi
			_TPI,		//2 body	N -> tau pi
			//kaon
			_EKA,		//2 body	N -> e K
			_MKA,		//2 body	N -> mu K
			//rho decay 100% in pions
			_nRHO0,		//2 body	N -> rho0
			_ERHO,		//2 body	N -> e rho
			_MRHO,		//2 body	M -> mu rho
			//kaon*
			_EKAx,		//2 body	N -> e K*
			_MKAx,		//2 body	N -> mu K*
			//other (eta, phi, omega.. )
			_nOMEGA,	//2 body	N -> nu w
			_nETA,		//2 body	N -> nu eta
			_nETAi,		//2 body	N -> nu eta'
			_nPHI,		//2 body	N -> nu phi
			//charm
			_ECHARM,	//2 body	N -> e D+
			//Channels for experimental comparison (EPI, MPI, nEE, nME, nMM)
			_ExpALL,	//

			//production modes
			//pure leptonic
			_MuonE,		//3 body	mu  -> nu N e	(via Ue)
			_MuonM,		//3 body	mu  -> nu N e	(via Um)
			_TauEE,		//3 body	tau -> nu N e	(via Ue)
			_TauET,		//3 body	tau -> nu N e	(via Ut)
			_TauMM,		//3 body	tau -> nu N mu	(via Um)
			_TauMT,		//3 body	tau -> nu N mu	(via Ut)
			_TauPI,		//2 body	tau -> N pi	(via Ut)
			_Tau2PI,	//3 body	tau -> N pi pi	(via Ut)
			//pseudomeson leptonic
			_PionE,		//2 body	pi -> N e	(via Ue)
			_PionM,		//2 body	pi -> N mu	(via Um)
			_KaonE,		//2 body	K  -> N e	(via Ue)
			_KaonM,		//2 body	K  -> N mu	(via Um)
			_CharmE,	//2 body	Ds -> N e	(via Ue)
			_CharmM,	//2 body	Ds -> N mu	(via Um)
			_CharmT,	//2 body	Ds -> N tau	(via Ut)
			//pseudomeson semileptonic
			_Kaon0E,	//3 body	K0 -> pi+ N e	(via Ue)
			_Kaon0M,	//3 body	K0 -> pi+ N mu	(via Um)
			_KaonCE,	//3 body	K+ -> pi0 N e	(via Ue)
			_KaonCM,	//3 body	K+ -> pi0 N mu	(via Um)
		};

		Amplitude();

		void LoadMap();
		std::string ShowChannel(Channel Name);

		Process LoadMass(Channel Name);

		double Kallen(double X, double Y, double Z);
		double SqrtKallen(double X, double Y, double Z);

		double dGammad5_3B(double M2);
		double dGammad2_3B(double M2);
		double dGammad2_2B(double M2, double x, double y);
		double dGammad0_2B(double M2, double x, double y);

		double Limit(double &s, double x, double y, double z);
		double Limit(double &s, double &t, double x, double y, double z);

		double M2_LeptonPseudoMeson(double cos0, double x, double y);
		double M2_NeutrinoPseudoMeson(double cos0, double x, double y);
		double M2_LeptonVectorMeson(double cos0, double x, double y);
		double M2_NeutrinoVectorMeson(double cos0, double x, double y);
		double M2_WW(double s, double cos0, double x, double y, double z);
		double M2_WZ(double s, double t, double cos0s, double cos0t, double x, double y, double z);
		double M2_WZ(double u, double cos0u, double x, double y, double z);

		double M2_LeptonNeutrino(double s, double x, double y, double z);
		double M2_AntiLeptonNeutrino(double s, double x, double y, double z);
		double M2_LeptonTwo(double x, double y);
		double M2_LeptonThree(double x, double y, double z);
		double M2_MesonTwo(double x, double y);
		double M2_MesonThree(double s, double t, double x, double y, double z, double L_, double L0);

		bool IsChanged();
		virtual void Reset() { ; }

		double Mass(int E = 1.0);
		double MassN(int E = 1.0);
		double Ue(int E = 1.0);
		double Um(int E = 1.0);
		double Ut(int E = 1.0);
		int Helicity();
		bool GetFermion();
		bool GetParticle();

		void SetMass(double Mass);
		void SetMassN(double Mass);
		void SetUe(double Ue);
		void SetUm(double Um);
		void SetUt(double Ut);
		void SetFermion(bool Fermion);
		void SetParticle(bool Particle);
		void SetHelicity(int Helicity);
		void SetNeutrino(double Mass, double* Mixings, bool Fermion, bool Particle, int Helix);

		double Function(double x);
		double Function_D(double *x);

		unsigned int CC;

	protected:
		Channel Channel_prev;

		std::vector<double> vMass;
		std::vector<int> vPdg;
		std::vector<double> F_var;

		std::map<Amplitude::Channel, std::string> chMap;

		void SetFunction(double (Amplitude::*FF)(double));
		void SetFunction_D(double (Amplitude::*FF)(double *));

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
		const double M_Etai;
		const double M_Phi;
		const double M_Tau;
		const double M_Charm;
		const double M_CharmS;

	private:
		double M_Sterile, M_Sterile_prev;
		double M_Parent;
		double fUe, fUm, fUt;

		bool bFermion, bParticle;
		int iHel, iHel_prev;

		double (Amplitude::*fFunction)(double);
		double (Amplitude::*fFunction_D)(double*);
};

#endif
