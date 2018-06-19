#include "Neutrino.h"

//Majorana can be treated as a neutrino + antineutrino
//
Neutrino::Neutrino(double Mass, unsigned int Options) :
	fMass(Mass)
{
	SetParticle(Options);
	SetFermion(Options);
	SetHelicity(Options);

	fMixings = new double[3];
	SetMixings(1.0, 1.0, 1.0);

	TheDecayRates = new DecayRates();	//to compute heavy neutrino decays
	TheProduction = new Production();	//to compute massive neutrino production widths
	TheProdLightN = new Production();	//to compute massless neutrino production widths
	ThePhaseSpace = new PhaseSpace();	//to generate phasespace for neutrino decays

	//TheCross = new CrossSection();
	TheProdLightN->SetNeutrino(0, Mixings(), 1, 1, 0);	//SM neutrino loaded
}

Neutrino::~Neutrino()
{
	delete TheDecayRates;
	delete TheProduction;
	delete ThePhaseSpace;
}

void Neutrino::SetParent(Amplitude *Object)
{
	Object->SetNeutrino(Mass(), Mixings(), GetFermion(), GetParticle(), Helicity());
}

//////////
///DECAY///
//////////
void Neutrino::DecayChannel(std::vector<std::string> &vChan)
{
	vChan = TheDecayRates->ListChannel();
}

double Neutrino::DecayTotal()
{
	return DecayWidth(Amplitude::_ALL);
}

double Neutrino::DecayWidth()
{
	return DecayWidth(DecayChannel());
}

double Neutrino::DecayWidth(std::string Name)
{
	return DecayWidth(TheDecayRates->FindChannel(Name));
}

double Neutrino::DecayWidth(Amplitude::Channel Name)
{
	SetParent(TheDecayRates);
	return TheDecayRates->Gamma(Name);
}

double Neutrino::DecayBranch()
{
	return DecayBranch(DecayChannel());
}

double Neutrino::DecayBranch(std::string Name)
{
	return DecayBranch(TheDecayRates->FindChannel(Name));
}

double Neutrino::DecayBranch(Amplitude::Channel Name)
{
	SetParent(TheDecayRates);
	return TheDecayRates->Branch(Name);
}

////////////////
///PRODUCTION///
////////////////
//
void Neutrino::ProductionChannel(std::vector<std::string> &vChan)
{
	vChan = TheProduction->ListChannel();
}

double Neutrino::ProductionWidth()
{
	return ProductionWidth(ProductionChannel());
}

double Neutrino::ProductionWidth(std::string Name)
{
	return ProductionWidth(TheProduction->FindChannel(Name));
}

double Neutrino::ProductionWidth(Amplitude::Channel Name)
{
	SetParent(TheProduction);
	return TheProduction->Gamma(Name);
}

double Neutrino::ProductionScale()
{
	return ProductionScale(ProductionChannel());
}

double Neutrino::ProductionScale(std::string Name)
{
	return ProductionScale(TheProduction->FindChannel(Name));
}

double Neutrino::ProductionScale(Amplitude::Channel Name)
{
	SetParent(TheProduction);
	return TheProduction->Gamma(Name)/(2*TheProdLightN->Gamma(Name));
}

std::vector<TLorentzVector> Neutrino::DecayPS()
{
	return GeneratePS(DecayChannel());
}

std::vector<TLorentzVector> Neutrino::ProductionPS()
{
	return GeneratePS(ProductionChannel());
}

std::vector<TLorentzVector> Neutrino::GeneratePS(Amplitude::Channel Name)
{
	SetParent(ThePhaseSpace);

	std::vector<TLorentzVector> vDaughter;
	if (ThePhaseSpace->Generate(Name))
	{
		for (unsigned int i = 0; i < ThePhaseSpace->Daughter(); ++i)
			vDaughter.push_back(ThePhaseSpace->Daughter(i));
	}

	return vDaughter;
}

void Neutrino::SetDecayChannel(std::string Name)
{
	chDecay = TheDecayRates->FindChannel(Name);
}

void Neutrino::SetProductionChannel(std::string Name)
{
	chProduction = TheProduction->FindChannel(Name);
}

Amplitude::Channel Neutrino::DecayChannel()
{
	return chDecay;
}

Amplitude::Channel Neutrino::ProductionChannel()
{
	return chProduction;
}

std::string Neutrino::DecayChannelName()
{
	return TheDecayRates->ShowChannel(DecayChannel());
}

std::string Neutrino::ProductionChannelName()
{
	return TheProduction->ShowChannel(ProductionChannel());
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
	fEnergy = Mass() + Energy;
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
double Neutrino::Mass()
{
	return fMass;
}

double* Neutrino::Mixings()
{
	return fMixings;
}

double Neutrino::Ue(int E)
{
	return pow(fMixings[0], E);
}

double Neutrino::Um(int E)
{
	return pow(fMixings[1], E);
}

double Neutrino::Ut(int E)
{
	return pow(fMixings[2], E);
}

double Neutrino::Energy()
{
	return fEnergy;
}

double Neutrino::EnergyKin()
{
	return fEnergy-Mass();
}

int Neutrino::Helicity()
{
	return iHel;
}

bool Neutrino::GetFermion()
{
	return bFermion;
}

bool Neutrino::IsDirac()
{
	return GetFermion();
}

bool Neutrino::IsMajorana()
{
	return !GetFermion();
}

bool Neutrino::GetParticle()
{
	return bParticle;
}

bool Neutrino::IsParticle()
{
	return GetParticle();
}

bool Neutrino::IsAntiparticle()
{
	return !GetParticle();
}
