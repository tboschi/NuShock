#include "PhaseSpace.h"

PhaseSpace::PhaseSpace(Neutrino *N) :
	Event(new TGenPhaseSpace),
	GenMT(new TRandom3(0))
{
	DecayRates();
	SetNeutrino(N);

	TheSpace = new TGenPhaseSpace;
	N_vec = new TLorentzVector;
	N_rest = new TLorentzVector;
}

double DecayRates::SetMass(Channel Name, std::vector<double> &vMass)	//return the maximum value of phase space
{
	double Max = 0.0;
	vMass.clear();

	switch(Name)
	{
		case _ALL:
		case _nnn:
			vMass.push_back(M_Neutrino);
			vMass.push_back(M_Neutrino);
			vMass.push_back(M_Neutrino);
			break;
		case _nGAMMA:
			vMass.push_back(M_Neutrino);
			vMass.push_back(M_Photon);
			break;
		case _nEE:
		case _ExpALL:
			vMass.push_back(M_Neutrino);
			vMass.push_back(M_Electron);
			vMass.push_back(M_Electron);

			Max = Max_nEEMax();
			break;
		case _nEMU:
			vMass.push_back(M_Neutrino);
			vMass.push_back(M_Electron);
			vMass.push_back(M_Muon);

			Max = Max_nEMU();
			break;
		case _nMUE:
			vMass.push_back(M_Neutrino);
			vMass.push_back(M_Electron);
			vMass.push_back(M_Muon);

			Max = Max_nMUE();
			break;
		case _nMUMU:
			vMass.push_back(M_Neutrino);
			vMass.push_back(M_Muon);
			vMass.push_back(M_Muon);

			Max = Max_nMUMU();
			break;
		case _nET:
			vMass.push_back(M_Neutrino);
			vMass.push_back(M_Electron);
			vMass.push_back(M_Tau);

			Max = Max_nET();
			break;
		case _nTE:
			vMass.push_back(M_Neutrino);
			vMass.push_back(M_Electron);
			vMass.push_back(M_Tau);

			Max = Max_nTE();
			break;
		case _nMUT:
			vMass.push_back(M_Neutrino);
			vMass.push_back(M_Muon);
			vMass.push_back(M_Tau);

			Max = Max_nMUT();
			break;
		case _nTMU:
			vMass.push_back(M_Neutrino);
			vMass.push_back(M_Muon);
			vMass.push_back(M_Tau);
			
			Max = Max_nTMU();
			break;
		case _nPI0:
			vMass.push_back(M_Neutrino);
			vMass.push_back(M_Pion0);
			break;
		case _EPI:
			vMass.push_back(M_Electron);
			vMass.push_back(M_Pion);
			break;
		case _MUPI:
			vMass.push_back(M_Muon);
			vMass.push_back(M_Pion);
			break;
		case _TPI:
			vMass.push_back(M_Pion);
			vMass.push_back(M_Tau);
			break;
		case _EKA:
			vMass.push_back(M_Electron);
			vMass.push_back(M_Kaon);
			break;
		case _MUKA:
			vMass.push_back(M_Muon);
			vMass.push_back(M_Kaon);
			break;
		case _EKAx:
			vMass.push_back(M_Electron);
			vMass.push_back(M_Kaonx);
			break;
		case _MUKAx:
			vMass.push_back(M_Muon);
			vMass.push_back(M_Kaonx);
			break;
		case _nRHO0:
			vMass.push_back(M_Neutrino);
			vMass.push_back(M_Rho0);
			break;
		case _ERHO:
			vMass.push_back(M_Electron);
			vMass.push_back(M_Rho);
			break;
		case _MURHO:
			vMass.push_back(M_Muon);
			vMass.push_back(M_Rho0);
			break;
		case _nETA:
			vMass.push_back(M_Neutrino);
			vMass.push_back(M_Eta);
			break;
		case _nETAi:
			vMass.push_back(M_Neutrino);
			vMass.push_back(M_Etai);
			break;
		case _nOMEGA:
			vMass.push_back(M_Neutrino);
			vMass.push_back(M_Omega);
			break;
		case _nPHI:
			vMass.push_back(M_Neutrino);
			vMass.push_back(M_Phi);
			break;
		case _ECHARM:
			vMass.push_back(M_Electron);
			vMass.push_back(M_Charm);
			break;
		default:
			break;
	}

	return Max;
}
bool PhaseSpace::SetDecay(TLorentzVector &Rest, std::vector<double> &vMass)
{
	double MassArray[10];
	nDaughter = vMass.size();
	for (unsigned int i = 0; i < nDaughter; ++i)
		MassArray[i] = Masses.at(i);

	return TheSpace->SetDecay(Rest, nDaughter, MassArray);
}

double PhaseSpace::LoadPhaseSpace(double Max)
{
	double / dMax;
	return TheSpace->Generate();
}

void PhaseSpace::Boost()
{
	for (unsigned int i = 0; i < nDaughter; ++i)
		TheSpace->GetDecayProduct
	N_vec->Boost();
}

double DecayRates::MaxThreeGamma()
{
	SetFunction(&I_WZ_s);
	Inte::MaxMin(this, Max, Min);
}


TLorentzVector *DecayRates::GetNvec()
{
	return N_vec;
}

TLorentzVector DecayRates::GetDecayProduct(int i, int &ID)
//Particle *DecayRates::GetDecayProduct(int i, int &ID)
{
	ID = PdgCode[i];
	TLorentzVector Daughter = *(Event->GetDecay(i));
	Daughter.Boost(N_vec->BoostVector());

	return Daughter;
}

void DecayRates::SetNvec(TLorentzVector &X)
{
	*N_vec = X;
	N_rest->SetPxPyPzE(0, 0, 0, N_vec->M());
}


double DecayRates::NeutrinoLeptonAB(double M_LeptonA, double M_LeptonB, double M_Neutrino)
{
	SetFunction(dNeutrinoLeptonAB);
}

double DecayRates::dNeutrinoLeptonAB(double M_LeptonA, double M_LeptonB, double M_Neutrino, double s, double cos0)
{
	double dMN2 = M_Neutrino*M_Neutrino/GetMass()/GetMass();
	double dMA2 = M_LeptonA*M_LeptonB/GetMass()/GetMass();
	double dMB2 = M_LeptonB*M_LeptonB/GetMass()/GetMass();

	return GetMult() * Const::fGF2 / (16.0 * Const::fPi3) * 
		pow(GetMass(), 5) * I_WW_s(dMN2, dMA2, dMB2);
}

double DecayRates::MaxNeutrinoLeptonAB(double M_LeptonA, double M_LeptonB, double M_Neutrino)
{
	double dMN2 = M_Neutrino*M_Neutrino/GetMass()/GetMass();
	double dMA2 = M_LeptonA*M_LeptonB/GetMass()/GetMass();
	double dMB2 = M_LeptonB*M_LeptonB/GetMass()/GetMass();

	return GetMult() * Const::fGF2 / (16.0 * Const::fPi3) * 
		pow(GetMass(), 5) * Max_WW(dMN2, dMA2, dMB2);
}

double DecayRates::Max_WW(double x, double y, double z)
{
	I_var.clear();

	I_var.push_back(x);	//0
	I_var.push_back(y);	//1
	I_var.push_back(z);	//2

	SetFunction(&I_WW_s);

	return Inte::Max(this);
}

double PhaseSpace::Max_WW_s(double s)
{
	I_var.push_back(s);	//3

	SetFunction(&I_WW_cos0);
	double Ret = Inte::Max(this);		//max in cos0 for given s

	I_var.pop_back();
	SetFunction(&I_WW_s);
	return Ret;
}

double PhaseSpace::Max_WW_cos0(double cos0)
{
	const double &x = I_var.at(0);
	const double &y = I_var.at(1);
	const double &z = I_var.at(2);
	const double &s = I_var.at(3);

	double sInf = x*x + y*y + 2*sqrt(x*y);
	double sSup = 1 - 2*sqrt(z) + z;

	double S = sInf + (sSup - sInf) * s;
	double Cos0 = -1 + 2*cos0;

	return I_WW(x, y, z, S, Cos0);
}

double DecayRates::I_WW(double x, double y, double z, double theta)
{
	I_var.clear();

	I_var.push_back(x);	//0
	I_var.push_back(y);	//1
	I_var.push_back(z);	//2
	I_var.push_back(cos0);	//3

	SetFunction(&I_WW_s);				//Function will fix the integration volume
	return Inte::BooleIntegration(this); 
}

double DecayRates::I_WW_s(double s)
{
	const double &x = I_var.at(0);
	const double &y = I_var.at(1);
	const double &z = I_var.at(2);
	const double &cos0 = I_var.at(3);

	return I_WW_s(s, cos0, x, y, z);
}

