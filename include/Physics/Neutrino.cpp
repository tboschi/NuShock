#include "Neutrino.h"

//Majorana can be treated as a neutrino + antineutrino
//
Neutrino::Neutrino(double Mass, unsigned int Options) :
	fMass(Mass)
{
	SetParticle(Options);
	SetFermion(Options);
	SetHelicity(Options);

	SetMixings(1.0, 1.0, 1.0);

	TheDecayRates = new DecayRates();
	TheProduction = new Production();
	ThePhaseSpace = new PhaseSpace();

	//TheCross = new CrossSection();
}

Neutrino::~Neutrino()
{
	delete TheDecayRates;
	delete TheProduction;
	delete ThePhaseSpace;
}

void Neutrino::SetParent(Amplitude *Object)
{
	Object->SetNeutrino(GetMass(), GetMixings(), GetFermion(), GetParticle(), GetHelicity());
}

double Neutrino::DecayWidth(Channel Name)
{
	SetParent(TheDecayRates);
	return TheDecayRates->Gamma(Name);
}

double Neutrino::DecayBranch(Channel Name)
{
	SetParent(TheDecayRates);
	return TheDecayRates->Branch(Name);
}

double Neutrino::ProductionWidth(Channel Name)
{
	SetParent(TheProduction);
	return TheProduction->Gamma(Channel);
}

double Neutrino::ProductionScale(Channel Name)
{
	SetParent(TheProduction);
	return TheProduction->Scale(Channel);
}

std::vector<TLorentzVector> Neutrino::PhaseSpace(Channel Name)
{
	SetParent(ThePhaseSpace);

	std::vector<TLorentzVector> vDaughter;
	if (ThePhaseSpace->Generate(Name))
	{
		for (unsigned int i = 0; i < TheSpace->GetNdaughter(); ++i)
			vDaughter.push_back(TheSpace->GetDaughter(i));
	}

	return vDaughter;
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

void Neutrino::SetEnergy(double Energy)
{
	fEnergy = Energy;
}

void Neutrino::SetEnergyKin(double Energy)
{
	fEnergy = GetMass() + Energy;
}

void Neutrino::SetHelicity(unsigned int Options)		//Left for particle is -1
{								//Right for particle is 1
	switch (Options & 3)					//For antiparticle, reverse
	{
		case 0:
			iHel = 1-2*GetParticle();
			break;
		case 1:
			iHel = 2*GetParticle()-1;
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

double Neutrino::GetEnergy()
{
	return fEnergy;
}

double Neutrino::GetEnergyKin()
{
	return fEnergy-GetMass();
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
