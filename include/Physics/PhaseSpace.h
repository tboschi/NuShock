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

class PhaseSpace : protected Amplitude
{
	public:
		PhaseSpace();

		bool SetDecay(Channel Name, TLorentzVector &Rest);
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
		double Max_NeutrinoLeptonAA(double &max_Ul, double &max_Ua, double M_Neut, double M_Lepton);
		double NeutrinoLeptonAB(double M_Neut, double M_LeptonA, double M_LeptonB);
		double Max_NeutrinoLeptonAB(double M_Neut, double M_LeptonA, double M_LeptonB);
		double max_NeutrinoLeptonLepton(double x, double y, double z, double gL, double gR);
		double max_NeutrinoLeptonLepton_D(double *x);
		double NeutrinoLeptonLepton(double x, double y, double z, double s, double t, double cos0, double cos1, double gL, double gR);

		double NeutrinoPseudoMeson(double M_Neut, double M_Meson);
		double Max_LeptonPseudoMeson(double M_Neut, double M_Meson);
		double LeptonPseudoMeson(double M_Neut, double M_Meson);
		double Max_NeutrinoPseudoMeson(double M_Neut, double M_Meson);
		double max_LeptonPseudo(double x, double y);
		double max_LeptonPseudo_cos0(double cos0);
		double LeptonPseudo(double x, double y, double cos0);

		double NeutrinoVectorMeson(double M_Neut, double M_Meson);
		double Max_NeutrinoVectorMeson(double M_Neut, double M_Meson);
		double LeptonVectorMeson(double M_Neut, double M_Meson);
		double Max_LeptonVectorMeson(double M_Neut, double M_Meson);
		double max_LeptonVector(double x, double y);
		double max_LeptonVector_cos0(double cos0);
		double LeptonVector(double x, double y, double cos0);

		void Kinematic_2B(double &cos0);
		void Kinematic_3B(double &s, double &t, double &cos0, double &cos1);
		void Boost();
		TLorentzVector *GetNvec();
		unsigned int GetDaughter();
		TLorentzVector GetDaughter(int i);
		void SetNLabf(TLorentzVector &Vec);
		void SetNRest(double Mass);
		void Reset();

	private:
		double fnnn,
                       fnGAMMA,
                       fnEE_e,	
                       fnEE_mt,	
                       fnEMU,	
                       fnMUE,	
                       fnMUMU_m,	
                       fnMUMU_et,	
                       fnET,	
                       fnTE,	
                       fnMUT,	
                       fnTMU,	
                       fnPI0,	
                       fEPI,	
                       fMUPI,	
                       fTPI,	
                       fEKA,	
                       fMUKA,	
                       fnRHO0,	
                       fERHO,	
                       fMURHO,	
                       fEKAx,	
                       fMUKAx,	
                       fnETA,	
                       fnETAi,	
                       fnOMEGA,
                       fnPHI,	
                       fECHARM;

		double maxnnn,
                       maxnGAMMA,
                       maxnEE_e,	
                       maxnEE_mt,	
                       maxnEMU,	
                       maxnMUE,	
                       maxnMUMU_m,	
                       maxnMUMU_et,	
                       maxnET,	
                       maxnTE,	
                       maxnMUT,	
                       maxnTMU,	
                       maxnPI0,	
                       maxEPI,	
                       maxMUPI,	
                       maxTPI,	
                       maxEKA,	
                       maxMUKA,	
                       maxnRHO0,	
                       maxERHO,	
                       maxMURHO,	
                       maxEKAx,	
                       maxMUKAx,	
                       maxnETA,	
                       maxnETAi,	
                       maxnOMEGA,
                       maxnPHI,	
                       maxECHARM;

	protected:
};

#endif
