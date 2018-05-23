#include "PhaseSpace.h"

PhaseSpace::PhaseSpace(Neutrino *N) :
	Event(new TGenPhaseSpace)
{
	DecayRates();
	SetNeutrino(N);

	Event = new TGenPhaseSpace;
	N_vec = new TLorentzVector;
	N_rest = new TLorentzVector;
}

PhaseSpace::SetDecay(TLorentzVector &Rest, std::vector<double> &vMass)
{
	double MassArray[10];
	int nDaughter = vMass.size();
	for (unsigned int i = 0; i < nDaughter; ++i)
		MassArray[i] = Masses.at(i);

	int Products;

	switch(Name)
	{
		case _ALL:
		case _nnn:
			Mass[0] = M_Neutrino;
			Result = nnn();
			break;
		case _nGAMMA:
			Result = nGAMMA();
			break;
		case _nEE:
			Result = nEE();
			break;
		case _nEMU:
			Result = nEMU();
			break;
		case _nMUE:
			Result = nMUE();
			break;
		case _nMUMU:
			Result = nMUMU();
			break;
		case _nET:
			Result = nET();
			break;
		case _nTE:
			Result = nTE();
			break;
		case _nMUT:
			Result = nMUT();
			break;
		case _nTMU:
			Result = nTMU();
			break;
		case _nPI0:
			Result = nPI0();
			break;
		case _EPI:
			Result = EPI();
			break;
		case _MUPI:
			Result = MUPI();
			break;
		case _TPI:
			Result = TPI();
			break;
		case _EKA:
			Result = EKA();
			break;
		case _MUKA:
			Result = MUKA();
			break;
		case _EKAx:
			Result = EKAx();
			break;
		case _MUKAx:
			Result = MUKAx();
			break;
		case _nRHO0:
			Result = nRHO0();
			break;
		case _ERHO:
			Result = ERHO();
			break;
		case _MURHO:
			Result = MURHO();
			break;
		case _nETA:
			Result = nETA();
			break;
		case _nETAi:
			Result = nETAi();
			break;
		case _nOMEGA:
			Result = nOMEGA();
			break;
		case _nPHI:
			Result = nPHI();
			break;
		case _ECHARM:
			Result = ECHARM();
			break;
		case _ExpALL:
			Result = ExpALL();
			break;
		default:
			Result = 0.0;
			break;
	}
}
