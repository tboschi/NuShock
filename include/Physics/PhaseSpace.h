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
		double nEM_ratio();
		double nME_ratio();
		double nET_ratio();
		double nTE_ratio();
		double nMT_ratio();
		double nTM_ratio();
		double nPI0_ratio();
		double nETA_ratio();
		double nETAi_ratio();
		double EPI_ratio();
		double MPI_ratio();
		double TPI_ratio();
		double EKA_ratio();
		double MKA_ratio();
		double ECHARM_ratio();
		double nRHO0_ratio();
		double nOMEGA_ratio();
		double nPHI_ratio();
		double ERHO_ratio();
		double MRHO_ratio();
		double EKAx_ratio();
		double MKAx_ratio();

		double NeutrinoLeptonAA(double &d_Ul, double &d_Un, double M_Neut, double M_Lepton);
		double NeutrinoLeptonAB(double M_Neut, double M_LeptonA, double M_LeptonB);
		double Max_NeutrinoLeptonAA(double &max_Ul, double &max_Ua, double M_Neut, double M_Lepton);
		double Max_NeutrinoLeptonAB(double M_Neut, double M_LeptonA, double M_LeptonB);
		double NeutrinoLeptonLepton(double x, double y, double z, double s, double t, double cos0, double cos1, double gL, double gR);
		double max_NeutrinoLeptonLepton(double x, double y, double z, double gL, double gR);
		double max_NeutrinoLeptonLepton_D(const double *p);

		double NeutrinoPseudoMeson(double M_Neut, double M_Meson);
		double Max_LeptonPseudoMeson(double M_Neut, double M_Meson);
		double LeptonPseudoMeson(double M_Neut, double M_Meson);
		double Max_NeutrinoPseudoMeson(double M_Neut, double M_Meson);
		double max_ToPseudoMeson(double x, double y);
		double max_ToPseudoMeson_cos0(const double cos0);
		double ToPseudoMeson(double x, double y, double cos0);

		double NeutrinoVectorMeson(double M_Neut, double M_Meson);
		double Max_NeutrinoVectorMeson(double M_Neut, double M_Meson);
		double LeptonVectorMeson(double M_Neut, double M_Meson);
		double Max_LeptonVectorMeson(double M_Neut, double M_Meson);
		double max_ToVectorMeson(double x, double y);
		double max_ToVectorMeson_cos0(const double cos0);
		double ToVectorMeson(double x, double y, double cos0);

		void Kinematic_2B(double &cos0);
		void Kinematic_3B(double &s, double &t, double &cos0, double &cos1);

		unsigned int Daughter();
		TLorentzVector* DaughterVector(unsigned int i, Reference Frame = RestFrame);
		Particle* Daughter(unsigned int i, Reference Frame = RestFrame);
		TLorentzVector *LabF();
		TLorentzVector *Rest();
		TLorentzVector *NuVector(Reference Frame = RestFrame);

		void SetNLabf(TLorentzVector &Vec);
		void SetNRest(double Mass);

		void Reset();
		void SetFunction(double (PhaseSpace::*FF)(const double));
		void SetFunction_D(double (PhaseSpace::*FF)(const double*));

	private:
		TGenPhaseSpace *Event;
		TRandom3 *GenMT;
		TLorentzVector *N_labf, *N_rest;
		unsigned int nDaughter;

		double (Amplitude::*M2_Function)(double, double, double);

		double fnnn,
                       fnGAMMA,
                       fnEE_e,	
                       fnEE_a,	
                       fnEM,	
                       fnME,	
                       fnMM_m,	
                       fnMM_a,	
                       fnET,	
                       fnTE,	
                       fnMT,	
                       fnTM,	
                       fnPI0,	
                       fEPI,	
                       fMPI,	
                       fTPI,	
                       fEKA,	
                       fMKA,	
                       fnRHO0,	
                       fERHO,	
                       fMRHO,	
                       fEKAx,	
                       fMKAx,	
                       fnETA,	
                       fnETAi,	
                       fnOMEGA,
                       fnPHI,	
                       fECHARM;

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
                       maxECHARM;

	protected:
};

#endif
