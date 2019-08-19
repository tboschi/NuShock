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
double Neutrino::DecayThreshold()
{
	if (DecayChannel() == Amplitude::_undefined)
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
		return DecayThreshold(DecayChannel());
}

double Neutrino::DecayThreshold(std::string Name)
{
	return DecayThreshold(theDecayRates->FindChannel(Name));
}

double Neutrino::DecayThreshold(Amplitude::Channel Name)
{
	SetParent(theDecayRates);
	return theDecayRates->MassThreshold(Name);
}

bool Neutrino::IsDecayAllowed()
{
	if (DecayChannel() == Amplitude::_undefined)
	{
		std::vector<Amplitude::Channel> vChan = theDecayRates->ListChannels();
		bool Ret = false;
		for (int ch = 0; ch < vChan.size(); ++ch)
			Ret += IsDecayAllowed(vChan.at(ch));

		return Ret;
	}
	else
		return IsDecayAllowed(DecayChannel());
}

bool Neutrino::IsDecayAllowed(std::string Name)
{
	return IsDecayAllowed(theDecayRates->FindChannel(Name));
}

bool Neutrino::IsDecayAllowed(Amplitude::Channel Name)
{
	SetParent(theDecayRates);
	return theDecayRates->IsAllowed(Name);
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

double Neutrino::DecayWidth()
{
	return DecayWidth(DecayChannel());
}

double Neutrino::DecayWidth(std::string Name)
{
	return DecayWidth(theDecayRates->FindChannel(Name));
}

double Neutrino::DecayWidth(Amplitude::Channel Name)
{
	SetParent(theDecayRates);
	//return (IsMajorana() ? 2.0 : 1.0) * theDecayRates->Gamma(Name);
	return theDecayRates->Gamma(Name);
}

double Neutrino::DecayBranch()
{
	return DecayBranch(DecayChannel());
}

double Neutrino::DecayBranch(std::string Name)
{
	return DecayBranch(theDecayRates->FindChannel(Name));
}

double Neutrino::DecayBranch(Amplitude::Channel Name)
{
	SetParent(theDecayRates);
	return theDecayRates->Branch(Name);
}

////////////////
///PRODUCTION///
////////////////
//
double Neutrino::ProductionThreshold()
{
	if (ProductionChannel() == Amplitude::_undefined)
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
		return ProductionThreshold(ProductionChannel());
}

double Neutrino::ProductionThreshold(std::string Name)
{
	return ProductionThreshold(theProduction->FindChannel(Name));
}

double Neutrino::ProductionThreshold(Amplitude::Channel Name)
{
	SetParent(theProduction);
	return theProduction->MassThreshold(Name);
}

bool Neutrino::IsProductionAllowed()
{
	if (ProductionChannel() == Amplitude::_undefined)
	{
		std::vector<Amplitude::Channel> vChan = theProduction->ListChannels();
		bool Ret = false;
		for (int ch = 0; ch < vChan.size(); ++ch)
			Ret += IsProductionAllowed(vChan[ch]);

		return Ret;
	}
	else
		return IsProductionAllowed(ProductionChannel());
}

bool Neutrino::IsProductionAllowed(std::string Name)
{
	return IsProductionAllowed(theProduction->FindChannel(Name));
}

bool Neutrino::IsProductionAllowed(Amplitude::Channel Name)
{
	SetParent(theProduction);
	return theProduction->IsAllowed(Name);
}

void Neutrino::ProductionChannels(std::vector<std::string> &vChan)
{
	vChan.clear();
	std::vector<Amplitude::Channel> vAmpChan = theProduction->ListChannels();
	for (int ch = 0; ch < vAmpChan.size(); ++ch)
		vChan.push_back(theProduction->FindChannel(vAmpChan[ch]));
}

double Neutrino::ProductionWidth()
{
	return ProductionWidth(ProductionChannel());
}

double Neutrino::ProductionWidth(std::string Name)
{
	return ProductionWidth(theProduction->FindChannel(Name));
}

double Neutrino::ProductionWidth(Amplitude::Channel Name)
{
	SetParent(theProduction);
	return theProduction->Gamma(Name);
}

double Neutrino::ProductionScale()
{
	return ProductionScale(ProductionChannel());
}

double Neutrino::ProductionScale(std::string Name)
{
	return ProductionScale(theProduction->FindChannel(Name));
}

double Neutrino::ProductionScale(Amplitude::Channel Name)
{
	SetParent(theProduction);
	//return theProduction->Gamma(Name, true) / theProdLightN->Gamma(Name) /
	//	(Helicity() ? 2.0 : 1.0);
	return theProduction->Gamma(Name, true) / theProdLightN->Gamma(Name);
}

std::vector<Particle> Neutrino::DecayPS()	//neutrino is labframe
{
	return DecayPS(DecayChannel());
}

std::vector<Particle> Neutrino::DecayPS(Amplitude::Channel Name)	//neutrino is labframe
{
	SetParent(thePhaseSpace);
	//TLorentzVector vec = FourVector();
	thePhaseSpace->SetLabFrame(FourVector());

	std::vector<Particle> daughters;
	if (thePhaseSpace->Generate(Name))
		for (int i = 0; i < thePhaseSpace->Daughters(); ++i)
			daughters.push_back(thePhaseSpace->Daughter(i, PhaseSpace::labFrame));
	else
		std::cout << "generation failed" << std::endl;

	return daughters;
}

std::vector<Particle> Neutrino::ProductionPS(TLorentzVector &Vec)	//other particle is labframe
{
	return ProductionPS(ProductionChannel(), Vec);
}

std::vector<Particle> Neutrino::ProductionPS(Amplitude::Channel Name, TLorentzVector &vec)
{
	SetParent(thePhaseSpace);
	thePhaseSpace->SetLabFrame(vec);

	std::vector<Particle> daughters;
	if (thePhaseSpace->Generate(Name))
		for (int i = 0; i < thePhaseSpace->Daughters(); ++i)
			daughters.push_back(thePhaseSpace->Daughter(i, PhaseSpace::labFrame));
	else
		std::cout << "generation failed" << std::endl;

	return daughters;
}

void Neutrino::SetDecayChannel(std::string Name)
{
	chDecay = theDecayRates->FindChannel(Name);
}

void Neutrino::SetProductionChannel(std::string Name)
{
	chProduction = theProduction->FindChannel(Name);
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
	return theDecayRates->ShowChannel(DecayChannel());
}

std::string Neutrino::ProductionChannelName()
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
