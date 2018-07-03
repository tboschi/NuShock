/*
 * PhaseSpace generator for simulation of sterile neutrino decays (it is not needed for production)
 * It also returns correct 
 *
 * Author: Tommaso Boschi
 */

#ifndef PHASESPACE_H
#define PHASESPACE_H

#include <iostream>
#include <fstream>
#include <vector>
#include <cstring>

//ROOT include
#include "TLorentzVector.h"
#include "TGenPhaseSpace.h"

#include "Physics/Amplitude.h"

class PhaseSpace : public Amplitude
{
	public:
		enum Reference
		{
			RestFrame = 0,
			LabFrame = 1,
		};

		PhaseSpace();

		bool SetDecay(Channel Name);
		bool Generate(Channel Name);
		double Ratio(Channel Name);

		double nEE_ratio();
		double nMM_ratio();
		double NeutrinoLeptonAA_ratio(double &maxName_u, double &maxName_a, double (PhaseSpace::*Uu)(int));

		double nEM_ratio();
		double nME_ratio();
		double nET_ratio();
		double nTE_ratio();
		double nMT_ratio();
		double nTM_ratio();
		double NeutrinoLeptonAB_ratio(double &maxName);

		double nPI0_ratio();
		double nETA_ratio();
		double nETAi_ratio();
		double NeutrinoPseudoMeson_ratio(double &maxName);

		double EPI_ratio();
		double MPI_ratio();
		double TPI_ratio();
		double EKA_ratio();
		double MKA_ratio();
		double ECHARM_ratio();
		double LeptonPseudoMeson_ratio(double &maxName);

		double nRHO0_ratio();
		double nOMEGA_ratio();
		double nPHI_ratio();
		double NeutrinoVectorMeson_ratio(double &manName);

		double ERHO_ratio();
		double MRHO_ratio();
		double EKAx_ratio();
		double MKAx_ratio();
		double LeptonVectorMeson_ratio(double &manName);

		double MuonE_ratio();
		double MuonM_ratio();
		double TauEE_ratio();
		double TauET_ratio();
		double TauMM_ratio();
		double TauMT_ratio();
		double Kaon0E_ratio();
		double Kaon0M_ratio();
		double KaonCE_ratio();
		double KaonCM_ratio();

		//DECAY RATES

		double NeutrinoLeptonAA(double &d_Ul, double &d_Un, double M_Neut, double M_Lepton);
		double NeutrinoLeptonAB(double M_Neut, double M_LeptonA, double M_LeptonB);
		double Max_NeutrinoLeptonAA(double &max_Ul, double &max_Ua, double M_Neut, double M_Lepton);
		double Max_NeutrinoLeptonAB(double M_Neut, double M_LeptonA, double M_LeptonB);
		double NeutrinoLeptonLepton(double s, double t, double cos0, double cos1, double x, double y, double z, double gL, double gR);
		double max_NeutrinoLeptonLepton(double x, double y, double z, double gL, double gR);
		double max_NeutrinoLeptonLepton_D(double *p);

		double NeutrinoPseudoMeson(double M_Neut, double M_Meson);
		double Max_LeptonPseudoMeson(double M_Neut, double M_Meson);
		double LeptonPseudoMeson(double M_Neut, double M_Meson);
		double Max_NeutrinoPseudoMeson(double M_Neut, double M_Meson);
		double max_ToPseudoMeson(double x, double y);
		double max_ToPseudoMeson_cos0(double cos0);
		double ToPseudoMeson(double cos0, double x, double y);

		double NeutrinoVectorMeson(double M_Neut, double M_Meson);
		double Max_NeutrinoVectorMeson(double M_Neut, double M_Meson);
		double LeptonVectorMeson(double M_Neut, double M_Meson);
		double Max_LeptonVectorMeson(double M_Neut, double M_Meson);
		double max_ToVectorMeson(double x, double y);
		double max_ToVectorMeson_cos0(double cos0);
		double ToVectorMeson(double cos0, double x, double y);

		//PRODUCTION

		double LeptonNeutrino(double M_Lepton0, double M_Lepton, double M_Neut);
		double Max_LeptonNeutrino(double M_Lepton0, double M_Lepton, double M_Neut);
		double max_LeptonNeutrino(double x, double y, double z);
		double max_LeptonNeutrino_u(double u);
		double AntiLeptonNeutrino(double M_Lepton0, double M_Lepton, double M_Neut);
		double Max_AntiLeptonNeutrino(double M_Lepton0, double M_Lepton, double M_Neut);
		double max_AntiLeptonNeutrino(double x, double y, double z);
		double max_AntiLeptonNeutrino_s(double s);

		double MesonThree(double M_Meson0, double M_Meson, double M_Lepton, double L_, double L0);
		double Max_MesonThree(double M_Meson0, double M_Meson, double M_Lepton, double L_, double L0);
		double max_MesonThree(double x, double y, double z, double L_, double L0);
		double max_MesonThree_D(double *p);

		//KINEMATICS

		void Kinematic_2B(double &cos0);
		void Kinematic_3B(double &s, double &t, double &cos0, double &cos1);

		unsigned int Daughter();
		TLorentzVector* DaughterVector(unsigned int i, Reference Frame = RestFrame);
		Particle* Daughter(unsigned int i, Reference Frame = RestFrame);
		TLorentzVector *LabF();
		TLorentzVector *Rest();
		TLorentzVector *Parent(Reference Frame = RestFrame);

		void SetLabf(TLorentzVector &Vec);
		void SetRest(double Mass);

		void Reset();
		void SetFunction(double (PhaseSpace::*FF)(double));
		void SetFunction_D(double (PhaseSpace::*FF)(double*));

	private:
		TGenPhaseSpace *Event;
		TRandom3 *GenMT;
		TLorentzVector *P_labf, *P_rest;
		unsigned int nDaughter;

		double (Amplitude::*M2_Function)(double, double, double);

		double maxnnn,
		       maxnGAMMA,
		       maxnEE_e,	
		       maxnEE_a,	
		       maxnEM,	
		       maxnME,	
		       maxnMM_m,	
		       maxnMM_a,	
		       maxnET,	
		       maxnTE,	
		       maxnMT,	
		       maxnTM,	
		       maxnPI0,	
		       maxEPI,	
		       maxMPI,	
		       maxTPI,	
		       maxEKA,	
		       maxMKA,	
		       maxnRHO0,	
		       maxERHO,	
		       maxMRHO,	
		       maxEKAx,	
		       maxMKAx,	
		       maxnETA,	
		       maxnETAi,	
		       maxnOMEGA,
		       maxnPHI,	
		       maxECHARM,
		       maxMuonE,
		       maxMuonM,
		       maxTauEE,
		       maxTauET,
		       maxTauMM,
		       maxTauMT,
		       maxTauPI,
		       maxTau2PI,
		       maxPionE,
		       maxPionM,
		       maxKaonE,
		       maxKaonM,
		       maxCharmE,
		       maxCharmM,
		       maxCharmT,
		       maxKaon0E,
		       maxKaon0M,
		       maxKaonCE,
		       maxKaonCM;

	protected:
};

#endif
