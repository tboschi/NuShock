/*
 * Decay rate simplifier for a three body decay of muon and kaon(0)
 * Author: Tommaso Boschi
 */

#ifndef THREEBODY_H
#define THREEBODY_H

#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <cstring>

#include "Tools/Tools.h"

class ThreeBody
{
	public:
		ThreeBody(std::string Parent, double MSterile = 0.0, double Ue = 0.0, double Um = 0.0, double Ut = 0.0);

		void InitMap();
		void InitConst();

		double ddGamma();
		double dGamma();
		double Gamma();
		double ddPhaseSpace();
		double dPhaseSpace();
		double PhaseSpace();
		double M2();
		double M2IntY();
		double M2IntXY();
		double M2_Z();
		double M2_ZIntY();
		double M2_WZ();
		double M2_WZIntY();

		double M2Lept();
		double M2LeptIntY();
		double M2Kaon();
		double M2KaonIntY(double Y);
		double M2KaonIntY();
		double M2Kaon0();
		double M2Kaon0IntY(double Y);
		double M2Kaon0IntY();
		double M2nEE();
		double M2nEEIntY();
		double M2nMUMU();
		double M2nEMU();
		double M2nMUE();
		double MaxGamma();

		double yLim(double &Min, double &Max);
		double xLim(double &Min, double &Max);
		double Integrate(double (ThreeBody::*FF)(), double A, double B);
		bool InLimX();
		bool InLimY();

		bool IsEnergyConserved();

		double fPlus();
		double fMinus();
		
		double a(double p = 1.0);
		double b(double p = 1.0);
		double c(double p = 1.0);
		double x(double p = 1.0);
		double y(double p = 1.0);

		void ElectronChannel();
		void MuonChannel();
		void TauChannel();

		std::string GetParent();
		double GetEnergyX();
		double GetEnergyY();
		double GetParentMass();
		double GetSterileMass();
		double GetUe();
		double GetUm();
		double GetUt();
		double GetUu();
		double GetDecayConst();
		double GetLambda1();
		double GetLambda0();

		void SetParent(std::string Name);
		void SetEnergyX(double X);
		void SetEnergyY(double X);
		void SetX(double X);
		void SetY(double X);
		void SetSterileMass(double X);
		void SetUe(double X);
		void SetUm(double X);
		void SetUt(double X);
		bool IsChanged();

	private:
		double M_Sterile, M_Parent;
		double M_Sterile_prev, M_Parent_prev;
		double U_e, U_m, U_t;
		double U_e_prev, U_m_prev, U_t_prev;

		const double M_Neutrino;
		const double M_Photon;
		const double M_Electron;
		const double M_Muon;
		const double M_Tau;
		const double M_Pion;
		const double M_Pion0;
		const double M_Kaon;
		const double M_Kaon0;

		const double fKaon;
		const double fLambda1;
		const double fLambda0;

		double fA, fB, fC;
		double fEX, fEY;
		double fMax;

		std::map<std::string, ChannelName> mapParent;
		std::string sParent;

		bool IsElectron, IsMuon, IsTau;
};

#endif
