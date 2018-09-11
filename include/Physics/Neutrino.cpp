#include "Physics/Neutrino.h"

//Majorana can be treated as a neutrino + antineutrino
//
Neutrino::Neutrino(double Mass, unsigned int Options) //:
	//fMass(Mass)
{
	SetMass(Mass);		//should initialise correctly
				//when SetEnergy is called, the neutrino is aligned on the z-axis (theta = 0, phi = 0)
	SetParticle(Options);
	SetFermion(Options);
	SetHelicity(Options);

	fMixings = new double[3];
	SetMixings(0.0, 0.0, 0.0);

	TheDecayRates = new DecayRates();	//to compute heavy neutrino decays, left helix
	TheProduction = new Production();	//to compute massive neutrino production widths
	TheProdLightN = new Production();	//to compute massless neutrino production widths
	ThePhaseSpace = new PhaseSpace();	//to generate phasespace for neutrino decays

	double MixOne[3] = {1.0, 1.0, 1.0};
	TheProdLightN->SetNeutrino(0, MixOne, 1, 1, -1);	//SM neutrino loaded

	//TheCross = new CrossSection();
	chDecay	     = Amplitude::_undefined;
	chProduction = Amplitude::_undefined;
}

Neutrino::~Neutrino()
{
	delete fMixings;
	fMixings = 0;

	delete TheDecayRates;
	TheDecayRates = 0;
	delete TheProduction;
	TheProduction = 0;
	delete TheProdLightN;
	TheProdLightN = 0;
	delete ThePhaseSpace;
	ThePhaseSpace = 0;
}

void Neutrino::SetParent(Amplitude *Object)
{
	Object->SetNeutrino(Mass(), Mixings(), GetFermion(), GetParticle(), Helicity());
}

//////////
///DECAY///
//////////
bool Neutrino::IsDecayAllowed()
{
	if (DecayChannel() == Amplitude::_undefined)
	{
		std::vector<Amplitude::Channel> vChan = TheDecayRates->ListChannels();
		bool Ret = false;
		for (unsigned int ch = 0; ch < vChan.size(); ++ch)
			Ret += IsDecayAllowed(vChan.at(ch));

		return Ret;
	}
	else
		return IsDecayAllowed(DecayChannel());
}

bool Neutrino::IsDecayAllowed(std::string Name)
{
	return IsDecayAllowed(TheDecayRates->FindChannel(Name));
}

bool Neutrino::IsDecayAllowed(Amplitude::Channel Name)
{
	SetParent(TheDecayRates);
	return TheDecayRates->IsAllowed(Name);
}

void Neutrino::DecayChannels(std::vector<std::string> &vChan)
{
	vChan.clear();
	std::vector<Amplitude::Channel> vAmpChan = TheDecayRates->ListChannels();
	for (unsigned int ch = 0; ch < vAmpChan.size(); ++ch)
		vChan.push_back(TheDecayRates->FindChannel(vAmpChan.at(ch)));
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
	return (IsMajorana() ? 2.0 : 1.0) * TheDecayRates->Gamma(Name);
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
bool Neutrino::IsProductionAllowed()
{
	if (ProductionChannel() == Amplitude::_undefined)
	{
		std::vector<Amplitude::Channel> vChan = TheProduction->ListChannels();
		bool Ret = false;
		for (unsigned int ch = 0; ch < vChan.size(); ++ch)
			Ret += IsProductionAllowed(vChan.at(ch));

		return Ret;
	}
	else
		return IsProductionAllowed(ProductionChannel());
}

bool Neutrino::IsProductionAllowed(std::string Name)
{
	return IsProductionAllowed(TheProduction->FindChannel(Name));
}

bool Neutrino::IsProductionAllowed(Amplitude::Channel Name)
{
	SetParent(TheProduction);
	return TheProduction->IsAllowed(Name);
}

void Neutrino::ProductionChannels(std::vector<std::string> &vChan)
{
	vChan.clear();
	std::vector<Amplitude::Channel> vAmpChan = TheProduction->ListChannels();
	for (unsigned int ch = 0; ch < vAmpChan.size(); ++ch)
		vChan.push_back(TheProduction->FindChannel(vAmpChan.at(ch)));
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
	//return TheProduction->Gamma(Name, true) / TheProdLightN->Gamma(Name) /
	//	(Helicity() ? 2.0 : 1.0);
	return TheProduction->Gamma(Name, true) / TheProdLightN->Gamma(Name);
}

std::vector<Particle*> Neutrino::DecayPS()	//neutrino is labframe
{
	return DecayPS(DecayChannel());
}

std::vector<Particle*> Neutrino::DecayPS(Amplitude::Channel Name)	//neutrino is labframe
{
	SetParent(ThePhaseSpace);
	TLorentzVector Vec = FourVector();
	ThePhaseSpace->SetLabf(Vec);

	std::vector<Particle*> vDaughter;
	if (ThePhaseSpace->Generate(Name))
		for (unsigned int i = 0; i < ThePhaseSpace->Daughters(); ++i)
			vDaughter.push_back(ThePhaseSpace->Daughter(i, PhaseSpace::LabFrame));

	return vDaughter;
}

std::vector<Particle*> Neutrino::ProductionPS(TLorentzVector &Vec)	//other particle is labframe
{
	return ProductionPS(ProductionChannel(), Vec);
}

std::vector<Particle*> Neutrino::ProductionPS(Amplitude::Channel Name, TLorentzVector &Vec)
{
	SetParent(ThePhaseSpace);
	ThePhaseSpace->SetLabf(Vec);

	std::vector<Particle*> vDaughter;
	if (ThePhaseSpace->Generate(Name))
		for (unsigned int i = 0; i < ThePhaseSpace->Daughters(); ++i)
			vDaughter.push_back(ThePhaseSpace->Daughter(i, PhaseSpace::LabFrame));

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
/*
void Neutrino::SetMass(double Mass)
{
	fMass = Mass;
}
*/

void Neutrino::SetMixings(double Ue, double Um, double Ut)
{
	fMixings[0] = Ue;
	fMixings[1] = Um;
	fMixings[2] = Ut;
}

/*
void Neutrino::SetEnergy(double Energy)
{
	fEnergy = Energy;
}

void Neutrino::SetEnergyKin(double Energy)
{
	fEnergy = Mass() + Energy;
}
*/

void Neutrino::SetHelicity(unsigned int Options)		//Left for particle is -1
{								//Right for particle is 1
	switch (Options & 3)
	{
		case 0:
			//iHel = 1-2*GetParticle();
			iHel = -1;
			break;
		case 1:
			//iHel = 2*GetParticle()-1;
			iHel =  1;
			break;
		case 2:
		case 3:
			iHel = 0;	//unpolarised
			break;
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

//this is just a flip of helicity actually
//if antiparticle -> flip helicity
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
/*
double Neutrino::Mass()
{
	return fMass;
}
*/

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

/*
double Neutrino::Energy()
{
	return fEnergy;
}

double Neutrino::EnergyKin()
{
	return fEnergy-Mass();
}
*/

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
