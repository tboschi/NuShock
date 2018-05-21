#include "Neutrino.h"

Neutrino::Neutrino(double Mass, unsigned int Options) :
	fMass = Mass
{
	SetParticle(Options);
	SetFermion(Options);
	SetHelicity(Options);

	SetMixings(1.0, 1.0, 1.0);

	TheDecay = new DecayRates();
	//TheCross = new CrossSection();
}

Neutrino::~Neutrino()
{
}

Neutrino::DecayWidth(Channel Name)
{
	TheDecay->SetParent(GetMass(), GetMixings(), GetFermion(), GetHelicity();
	return TheDecay->Gamma(Channel);
}

Neutrino::DecayBranch(Channel Name)
{
	TheDecay->SetParent(GetMass(), GetMixings(), GetFermion(), GetHelicity();
	return TheDecay->Branch(Channel);
}

//setter
//
void Neutrino::SetMass(double Mass)
{
	fMass = Mass;
}

void Neutrino::SetMixings(double Ue, double Um, double Ut)
{
	fMixings[0] = Ue;
	fMixings[1] = Um;
	fMixings[2] = Ut;
}

void Neutrino::SetHelicity(unsigned int Options)
{
	switch (Options & 3)
	{
		case 0:
			iHel = -1;
			break;
		case 1:
			iHel = 1;
			break;
		case 2:
		case 3:
			iHel = 0;	//unpolarised
		default:
			std::cerr << "SetHelicity " << Options << "\t:\t";
			std::cerr << "Invalid options" << std::endl;
			break;
	}
}

void Neutrino::SetFermion(unsigned int Options)
{
	switch (Options & 4)
	{
		case 0:
			bFermion = true;	//Dirac
			break;
		case 4:
			bFermion = false;	//Majorana
			break;
		default:
			std::cerr << "SetFermion " << Options << "\t:\t";
			std::cerr << "Invalid options" << std::endl;
			break;
	}
}

void Neutrino::SetParticle(unsigned int Options)
{
	switch (Options & 8)
	{
		case 0:
			bParticle = true;	//Particle
			break;
		case 8:
			bParticle = false;	//Antiparticle
			break;
		default:
			std::cerr << "SetParticle " << Options << "\t:\t";
			std::cerr << "Invalid options" << std::endl;
			break;
	}
}


//getter
//
double Neutrino::GetMass()
{
	return fMass;
}

double* Neutrino::GetMixings()
{
	return fMixings;
}

int Neutrino::GetHelicity()
{
	return iHel;
}

bool Neutrino::IsDirac()
{
	return bFermion;
}

bool Neutrino::IsMajorana()
{
	return !bFermion;
}

bool Neutrino::IsParticle()
{
	return bParticle;
}

bool Neutrino::IsAntiparticle()
{
	return !bParticle;
}
