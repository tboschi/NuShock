#include "Neutrino.h"

//Majorana can be treated as a neutrino + antineutrino
//
Neutrino::Neutrino(double mass, size_t opts)
{
	SetM(mass);
	SetOptions(opts);
}

std::ostream & operator<<(std::ostream &os, const Neutrino &N) {
	return os << "<Neutrino: mass " << N.M() << ", helicity = " << N.Helicity()
		  << ", " << (IsDirac() ? "dirac" : "majorana")
		  << ", vec (" << N.E() << ", " << N.Px() << ", " << N.Py() << ", " << N.Pz() << ")>";
}


// neutrinos are identical if masses and quantum numbers are
Neutrino & Neutrino::operator==(const Neutrino & rhs) const 
{
	return ( (M() == rhs.M()) && (_opts == rhs._opts) );
}

Neutrino & Neutrino::operator!=(const Neutrino & rhs) const 
{
	return !(*this == rhs);
}

int Neutrino::Helicity() const {
	return Neutrino::Helicity(_opts);
}

bool Neutrino::IsDirac() const {
	return Neutrino::IsDirac(_opts);
}

bool Neutrino::IsMajorana() const {
	return Neutrino::IsMajorana(_opts);
}

bool Neutrino::IsParticle() const {
	return Neutrino::IsParticle(_opts);
}

bool Neutrino::IsAntiparticle() const {
	return Neutrino::IsAntiParticle(_opts);
}

void Neutrino::SetOptions(const size_t &opts) {
	_opts = opts & 15;
}

void Neutrino::AddOptions(const size_t &opts) {
	_opts |= opts;
}



void Neutrino::SetParent(Amplitude *Object)
{
	Object->SetNeutrino(M(), _mixings, _opts);
}

/*
//////////
///DECAY///
//////////

Neutrino::Neutrino(double mass, unsigned int opts) :
	Particle(),
	_opts(opts),
	_decay(Channel::_undefined),
	_production(Channel::_undefined)
{
	SetM(mass);

	theDecayRates = new DecayRates();	//to compute heavy neutrino decays, left helix
	theProduction = new Production();	//to compute massive neutrino production widths
	theProdLightN = new Production();	//to compute massless neutrino production widths
	thePhaseSpace = new PhaseSpace();	//to generate phasespace for neutrino decays

	double MixOne[3] = {1.0, 1.0, 1.0};
	theProdLightN->SetNeutrino(0, MixOne, 1, 1, -1);	//SM neutrino loaded
}

Neutrino::Neutrino(const Neutrino &N)
	_opts(N._opts),
	_mixings(N._mixings)
	_decay(N._decay),
	_production(N._production)
{
	SetMass(N.M());

			
	SetMixings(N.Ue(), N.Um(), N.Ut());

	theDecayRates = new DecayRates();	//to compute heavy neutrino decays, left helix
	theProduction = new Production();	//to compute massive neutrino production widths
	theProdLightN = new Production();	//to compute massless neutrino production widths
	thePhaseSpace = new PhaseSpace();	//to generate phasespace for neutrino decays

	double MixOne[4] = {1.0, 1.0, 1.0};
	theProdLightN->SetNeutrino(0, MixOne, 1, 1, -1);	//SM neutrino loaded

	//theCross = new CrossSection();
	_decay	     = N._decay;
	_production = N._production;
}
double Neutrino::DecayThreshold() {
	return DecayThreshold(_decay);
}

double Neutrino::DecayThreshold(const std::string &name)
{
	if (name.empty())
		return DecayThreshold();
	else
		return DecayThreshold(Channel::fromString(name));
}

double Neutrino::DecayThreshold(Channel::Name chan)
{
	if (chan == Channel::_undefined) {
		double limit = Const::MZ;
		for (const auto chan : Channel::Decays ) {
			double tmp = DecayThreshold(chan);
			if (tmp < limit)
				limit = tmp;
		}

		return limit;
	}
	else {
		SetParent(theDecayRates);
		return theDecayRates->MassThreshold(chan);
	}
}

bool Neutrino::IsDecayAllowed() {
	return IsDecayAllowed(_decay);
}

bool Neutrino::IsDecayAllowed(const std::string &name)
{
	if (name.empty())
		return IsDecayAllowed();
	else
		return IsDecayAllowed(Channel::fromString(name));
}

bool Neutrino::IsDecayAllowed(Channel::Name chan)
{
	if (chan == Channel::_undefined) {
		return std::any_of(std::begin(Channel::Decays),
				   std::end(Channel::Decays),
				[](const Channel::Name &chan) { return IsDecayAllowed(chan); } );
	}
	else {
		SetParent(theDecayRates);
		return theDecayRates->IsAllowed(chan);
	}
}

std::vector<std::string> Neutrino::DecayChannels()
{
	std::vector<std::string> channels(std::size(Channel::Decays));
	std::transform(std::begin(Channel::Decays), std::end(Channel::Decays), channels.begin(),
			[](const Channel::Name &chan) { return Channel::toString(chan); });

	return channels;
}

double Neutrino::DecayTotal()
{
	return DecayWidth(Channel::_ALL);
}

double Neutrino::DecayWidth() {
	DecayWidth(_decay);
}

double Neutrino::DecayWidth(const std::string &name)
{
	if (name.empty())
		return DecayWidth();
	else
		return DecayWidth(Channel::fromString(name));
}

double Neutrino::DecayWidth(Channel::Channel chan)
{
	SetParent(theDecayRates);
	return theDecayRates->Gamma(chan);
}

double Neutrino::DecayBranch()
{
	return DecayBranch(_decay);
}

double Neutrino::DecayBranch(const std::string &name)
{
	if (name.empty())
		return DecayBranch();
	else
		return DecayBranch(Channel::fromString(name));
}

double Neutrino::DecayBranch(Channel::Channel name)
{
	SetParent(theDecayRates);
	return theDecayRates->Branch(name);
}

////////////////
///PRODUCTION///
////////////////
//

double Neutrino::ProductionThreshold() {
	return ProductionThreshold(_production);
}

double Neutrino::ProductionThreshold(const std::string &name)
{
	if (name.empty())
		return ProductionThreshold();
	else
		return ProductionThreshold(Channel::fromString(name));
}

double Neutrino::ProductionThreshold(Channel::Name chan)
{
	if (chan == Channel::_undefined) {
		double limit = 0.;
		for (const auto chan = Channel::Productions)
			double tmp = ProductionThreshold(chan);
			if (tmp > limit)
				limit = tmp;
		}

		return limit;
	}
	else {
		SetParent(theProduction);
		return theProduction->MassThreshold(name);
	}
}

bool Neutrino::IsProductionAllowed() {
	IsProductionAllowed(_production);
}

bool Neutrino::IsProductionAllowed(const std::string &name)
{
	if (name.empty())
		return IsProductionAllowed();
	else
		return IsProductionAllowed(Channel::fromString(name));
}

bool Neutrino::IsProductionAllowed(Channel::Name name)
{
	if (name == Channel::_undefined) {
		return std::any_of(std::begin(Channel::Productions),
				   std::end(Channel::Productions),
				   [](const Channel::Name &chan) {
				   	return IsProductionAllowed(chan); });
	}
	else {
		SetParent(theProduction);
		return theProduction->IsAllowed(name);
	}
}

std::vector<std::string> Neutrino::ProductionChannels()
{
	std::vector<std::string> channels(std::size(Channel::Productions));
	std::transform(std::begin(Channel::Productions), std::end(Productions), channels.begin(),
			[](const Channel::Name &chan) { return Channel::toString(chan); });

	return channels;
}

double Neutrino::ProductionWidth()
{
	return ProductionWidth(_production);
}

double Neutrino::ProductionWidth(const std::string &name)
{
	if (name.empty())
		return ProductionWidth();
	else
		return ProductionWidth(Channel::fromString(name));
}

double Neutrino::ProductionWidth(Channel::Name name)
{
	SetParent(theProduction);
	return theProduction->Gamma(name);
}

double Neutrino::ProductionScale() {
	return ProductionScale(_production);
}

double Neutrino::ProductionScale(const std::string &name)
{
	if (name.empty())
		return ProductionScale();
	else
		return ProductionScale(Channel::fromString(name));
}

double Neutrino::ProductionScale(Channel::Name chan)
{
	SetParent(theProduction);
	//return theProduction->Gamma(name, true) / theProdLightN->Gamma(name) /
	//	(Helicity() ? 2.0 : 1.0);
	return theProduction->Gamma(chan, true) / theProdLightN->Gamma(chan);
}

std::vector<Particle> Neutrino::DecayPS()	//neutrino is labframe
{
	return DecayPS(_decay);
}

std::vector<Particle> Neutrino::DecayPS(std::string name)	//neutrino is labframe
{
	if (name.empty())
		return DecayPS();
	else
		return DecayPS(Channel::toString(name));
}

std::vector<Particle> Neutrino::DecayPS(Channel::Name name)	//neutrino is labframe
{
	SetParent(thePhaseSpace);
	//TLorentzVector vec = FourVector();
	thePhaseSpace->SetLabFrame(FourVector());

	double val;
	if (thePhaseSpace->Generate(name, val)) {
		std::vector<Particle> daughters;
		daughters.reserve(thePhaseSpace->Daughters();
		for (size_t i = 0; i < thePhaseSpace->Daughters(); ++i)
			daughters.push_back(thePhaseSpace->Daughter(i, PhaseSpace::labFrame));
	}
	else
		std::cout << "generation failed " << val << std::endl;

	return std::vector<Particle>();
}

std::vector<Particle> Neutrino::ProductionPS(const TLorentzVector &vec)	//other particle is labframe
{
	return ProductionPS(vec, _production);
}


std::vector<Particle> Neutrino::ProductionPS(const TLorentzVector &vec, std::string name)	//other particle is labframe
{
	if (name.empty())
		return ProductionPS(vec);
	else
		return ProductionPS(vec, Channel::fromString(name));
}

std::vector<Particle> Neutrino::ProductionPS(const TLorentzVector &vec, Channel::Name name)
{
	SetParent(thePhaseSpace);
	thePhaseSpace->SetLabFrame(vec);

	double val;
	if (thePhaseSpace->Generate(name, val)) {
		std::vector<Particle> daughters;
		daughters.reserve(thePhaseSpace->Daughters();
		for (size_t i = 0; i < thePhaseSpace->Daughters(); ++i)
			daughters.push_back(thePhaseSpace->Daughter(i, PhaseSpace::labFrame));
	}
	else
		std::cout << "generation failed " << val << std::endl;

	return std::vector<Particle>();
}

void Neutrino::SetDecayChannel(const std::string &name)
{
	_decay = Channel::fromString(name);
}

void Neutrino::SetProductionChannel(const std::string &name)
{
	_production = Channel::fromString(name);
}

Channel::Name Neutrino::DecayChannel() const
{
	return _decay;
}

Channel::Name Neutrino::ProductionChannel() const
{
	return _production;
}

std::string Neutrino::DecayChannelName() const
{
	return Channel::toString(_decay);
}

std::string Neutrino::ProductionChannelName() const
{
	return Channel::toString(_production);
}

//setter
//
void Neutrino::SetMass(double Mass)
{
	fMass = Mass;
}

void Neutrino::SetMixings(double Ue, double Um, double Ut)
{
	_mixings[0] = Ue;
	_mixings[1] = Um;
	_mixings[2] = Ut;
}

void Neutrino::SetEnergy(double Energy)
{
	fEnergy = Energy;
}

void Neutrino::SetEnergyKin(double Energy)
{
	fEnergy = Mass() + Energy;
}

//this is just a flip of helicity actually
//if antiparticle -> flip helicity
void Neutrino::SetOptions(const size_t &opts)
{
	_opts = opts & 15;
}

void Neutrino::AddOptions(const size_t &opts)
{
	_opts |= opts;
}


double* Neutrino::Mixings()
{
	return _mixings.data();
}

double Neutrino::Ue(int e) const
{
	return std::pow(_mixings[0], e);
}

double Neutrino::Um(int E) const
{
	return std::pow(_mixings[1], e);
}

double Neutrino::Ut(int E) const
{
	return std::pow(_mixings[2], e);
}
*/
