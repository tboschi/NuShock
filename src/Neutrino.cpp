#include "Neutrino.h"

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

	theDecayRates = new DecayRates();	//to compute heavy neutrino decays, left helix
	theProduction = new Production();	//to compute massive neutrino production widths
	theProdLightN = new Production();	//to compute massless neutrino production widths
	thePhaseSpace = new PhaseSpace();	//to generate phasespace for neutrino decays

	double MixOne[3] = {1.0, 1.0, 1.0};
	theProdLightN->SetNeutrino(0, MixOne, 1, 1, -1);	//SM neutrino loaded

	//theCross = new CrossSection();
	chDecay	     = Amplitude::_undefined;
	chProduction = Amplitude::_undefined;
}

Neutrino::Neutrino(const Neutrino &N)
{
	SetMass(N.Mass());
			
	bParticle = N.bParticle;
	bFermion = N.bFermion;
	iHel = N.iHel;

	fMixings = new double[3];
	SetMixings(N.Ue(), N.Um(), N.Ut());

	theDecayRates = new DecayRates();	//to compute heavy neutrino decays, left helix
	theProduction = new Production();	//to compute massive neutrino production widths
	theProdLightN = new Production();	//to compute massless neutrino production widths
	thePhaseSpace = new PhaseSpace();	//to generate phasespace for neutrino decays

	double MixOne[4] = {1.0, 1.0, 1.0};
	theProdLightN->SetNeutrino(0, MixOne, 1, 1, -1);	//SM neutrino loaded

	//theCross = new CrossSection();
	chDecay	     = N.chDecay;
	chProduction = N.chProduction;
}

Neutrino::~Neutrino()
{
	delete fMixings;

	delete theDecayRates;
	delete theProduction;
	delete theProdLightN;
	delete thePhaseSpace;
}

Neutrino & Neutrino::operator=(const Neutrino & N)
{
	if (this != &N)
	{
		SetMass(N.Mass());

		bParticle = N.bParticle;
		bFermion  = N.bFermion;
		iHel	  = N.iHel;

		chDecay	     = N.chDecay;
		chProduction = N.chProduction;

		delete theDecayRates;
		delete theProduction;
		delete theProdLightN;
		delete thePhaseSpace;
		theDecayRates = new DecayRates();
		theProduction = new Production();
		theProdLightN = new Production();
		thePhaseSpace = new PhaseSpace();

		double MixOne[3] = {1.0, 1.0, 1.0};
		theProdLightN->SetNeutrino(0, MixOne, 1, 1, -1);	//SM neutrino loaded

		delete fMixings;
		fMixings = new double[3];
		SetMixings(N.Ue(), N.Um(), N.Ut());
	}

	return *this;
}

void Neutrino::SetParent(Amplitude *Object)
{
	Object->SetNeutrino(Mass(), Mixings(), GetFermion(), GetParticle(), Helicity());
}

//////////
///DECAY///
//////////

double Neutrino::DecayThreshold(std::string name)
{
	if (name.empty())
		return DecayThreshold(DecayChannel());
	else
		return DecayThreshold(theDecayRates->FindChannel(name));
}

double Neutrino::DecayThreshold(Amplitude::Channel name)
{
	if (name == Amplitude::_undefined)
	{
		std::vector<Amplitude::Channel> vChan = theDecayRates->ListChannels();
		double Limit = Const::MZ;

		for (int ch = 0; ch < vChan.size(); ++ch)
		{
			double tmp = DecayThreshold(vChan[ch]);
			if (tmp < Limit)
				Limit = tmp;
		}

		return Limit;
	}
	else
	{
		SetParent(theDecayRates);
		return theDecayRates->MassThreshold(name);
	}
}

bool Neutrino::IsDecayAllowed(std::string name)
{
	if (name.empty())
		return IsDecayAllowed(DecayChannel());
	else
		return IsDecayAllowed(theDecayRates->FindChannel(name));
}

bool Neutrino::IsDecayAllowed(Amplitude::Channel name)
{
	if (name == Amplitude::_undefined)
	{
		std::vector<Amplitude::Channel> vChan = theDecayRates->ListChannels();
		bool Ret = false;
		for (int ch = 0; ch < vChan.size(); ++ch)
			Ret += IsDecayAllowed(vChan.at(ch));

		return Ret;
	}
	else
	{
		SetParent(theDecayRates);
		return theDecayRates->IsAllowed(name);
	}
}

void Neutrino::DecayChannels(std::vector<std::string> &vChan)
{
	vChan.clear();
	std::vector<Amplitude::Channel> vAmpChan = theDecayRates->ListChannels();
	for (int ch = 0; ch < vAmpChan.size(); ++ch)
		vChan.push_back(theDecayRates->FindChannel(vAmpChan.at(ch)));
}

double Neutrino::DecayTotal()
{
	return DecayWidth(Amplitude::_ALL);
}

double Neutrino::DecayWidth(std::string name)
{
	if (name.empty())
		return DecayWidth(DecayChannel());
	else
		return DecayWidth(theDecayRates->FindChannel(name));
}

double Neutrino::DecayWidth(Amplitude::Channel name)
{
	SetParent(theDecayRates);
	return theDecayRates->Gamma(name);
}

double Neutrino::DecayBranch(std::string name)
{
	if (name.empty())
		return DecayBranch(DecayChannel());
	else
		return DecayBranch(theDecayRates->FindChannel(name));
}

double Neutrino::DecayBranch(Amplitude::Channel name)
{
	SetParent(theDecayRates);
	return theDecayRates->Branch(name);
}

////////////////
///PRODUCTION///
////////////////
//

double Neutrino::ProductionThreshold(std::string name)
{
	if (name.empty())
		return ProductionThreshold(ProductionChannel());
	else
		return ProductionThreshold(theProduction->FindChannel(name));
}

double Neutrino::ProductionThreshold(Amplitude::Channel name)
{
	if (name == Amplitude::_undefined)
	{
		std::vector<Amplitude::Channel> vChan = theProduction->ListChannels();
		double Limit = 0.0;

		for (int ch = 0; ch < vChan.size(); ++ch)
		{
			double tmp = ProductionThreshold(vChan[ch]);
			if (tmp > Limit)
				Limit = tmp;
		}

		return Limit;
	}
	else
	{
		SetParent(theProduction);
		return theProduction->MassThreshold(name);
	}
}

bool Neutrino::IsProductionAllowed(std::string name)
{
	if (name.empty())
		return IsProductionAllowed(ProductionChannel());
	else
		return IsProductionAllowed(theProduction->FindChannel(name));
}

bool Neutrino::IsProductionAllowed(Amplitude::Channel name)
{
	if (name == Amplitude::_undefined)
	{
		std::vector<Amplitude::Channel> vChan = theProduction->ListChannels();
		bool Ret = false;
		for (int ch = 0; ch < vChan.size(); ++ch)
			Ret += IsProductionAllowed(vChan[ch]);

		return Ret;
	}
	else
	{
		SetParent(theProduction);
		return theProduction->IsAllowed(name);
	}
}

void Neutrino::ProductionChannels(std::vector<std::string> &vChan)
{
	vChan.clear();
	std::vector<Amplitude::Channel> vAmpChan = theProduction->ListChannels();
	for (int ch = 0; ch < vAmpChan.size(); ++ch)
		vChan.push_back(theProduction->FindChannel(vAmpChan[ch]));
}

double Neutrino::ProductionWidth(std::string name)
{
	if (name.empty())
		return ProductionWidth(ProductionChannel());
	else
		return ProductionWidth(theProduction->FindChannel(name));
}

double Neutrino::ProductionWidth(Amplitude::Channel name)
{
	SetParent(theProduction);
	return theProduction->Gamma(name);
}

double Neutrino::ProductionScale(std::string name)
{
	if (name.empty())
		return ProductionScale(ProductionChannel());
	else
		return ProductionScale(theProduction->FindChannel(name));
}

double Neutrino::ProductionScale(Amplitude::Channel name)
{
	SetParent(theProduction);
	//return theProduction->Gamma(name, true) / theProdLightN->Gamma(name) /
	//	(Helicity() ? 2.0 : 1.0);
	return theProduction->Gamma(name, true) / theProdLightN->Gamma(name);
}

std::vector<Particle> Neutrino::DecayPS(std::string name)	//neutrino is labframe
{
	if (name.empty())
		return DecayPS(DecayChannel());
	else
		return DecayPS(theDecayRates->FindChannel(name));
}

std::vector<Particle> Neutrino::DecayPS(Amplitude::Channel name)	//neutrino is labframe
{
	SetParent(thePhaseSpace);
	//TLorentzVector vec = FourVector();
	thePhaseSpace->SetLabFrame(FourVector());

	std::vector<Particle> daughters;
	double val;
	if (thePhaseSpace->Generate(name, val))
		for (int i = 0; i < thePhaseSpace->Daughters(); ++i)
			daughters.push_back(thePhaseSpace->Daughter(i, PhaseSpace::labFrame));
	else
		std::cout << "generation failed " << val << std::endl;

	return daughters;
}

std::vector<Particle> Neutrino::ProductionPS(const TLorentzVector &vec, std::string name)	//other particle is labframe
{
	if (name.empty())
		return ProductionPS(vec, DecayChannel());
	else
		return ProductionPS(vec, theProduction->FindChannel(name));
}

std::vector<Particle> Neutrino::ProductionPS(const TLorentzVector &vec, Amplitude::Channel name)
{
	SetParent(thePhaseSpace);
	thePhaseSpace->SetLabFrame(vec);

	std::vector<Particle> daughters;
	double val;
	if (thePhaseSpace->Generate(name, val))
		for (int i = 0; i < thePhaseSpace->Daughters(); ++i)
			daughters.push_back(thePhaseSpace->Daughter(i, PhaseSpace::labFrame));
	else
		std::cout << "generation failed " << val << std::endl;

	return daughters;
}

void Neutrino::SetDecayChannel(std::string name)
{
	chDecay = theDecayRates->FindChannel(name);
}

void Neutrino::SetProductionChannel(std::string name)
{
	chProduction = theProduction->FindChannel(name);
}

Amplitude::Channel Neutrino::DecayChannel() const
{
	return chDecay;
}

Amplitude::Channel Neutrino::ProductionChannel() const
{
	return chProduction;
}

std::string Neutrino::DecayChannelName() const
{
	return theDecayRates->ShowChannel(DecayChannel());
}

std::string Neutrino::ProductionChannelName() const
{
	return theProduction->ShowChannel(ProductionChannel());
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

double Neutrino::Ue(int E) const
{
	return pow(fMixings[0], E);
}

double Neutrino::Um(int E) const
{
	return pow(fMixings[1], E);
}

double Neutrino::Ut(int E) const
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
