#include "physics/Particle.h"

Particle::Particle(int pdgCode, double E, double px, double py, double pz) :
	TLorentzVector(px, py, pz, E),
	_pdg(pdgCode)
	//_charge(0)
{
	ChargeID();
}

Particle::Particle(int pdgCode, const TLorentzVector &fv) :
	TLorentzVector(fv),
	_pdg(pdgCode)
	//_charge(0)
{
	ChargeID();
}

std::ostream & operator<<(std::ostream &os, const Particle &p) {
	return os << "<pdg=" << p.Pdg() << ", M=" << p.M()
		  << ", Q=" << p.Q() << ", p=(" << p.E()
		  << ", " << p.Px() << ", " << p.Py()
		  << ", " << p.Pz() << ")>";
}


/*
//copy ctor
Particle::Particle(const Particle &p) :
	_pdg(p._pdg),
	TLorentzVector(p)
{
}
*/

// return apparent pdg (multplied by recorded charge)
int Particle::Pdg() const
{
	// opposite particle
	if (RealQ() * _charge < 0)
		return - _pdg;
	return _pdg;
}

// pdg assigned at creation of object
int Particle::RealPdg() const
{
	return _pdg;
}

// recorded charge
int Particle::Q() const
{
	return _charge;
}

// true charge associated to pdg code
int Particle::RealQ() const
{
	return Particle::Q(_pdg);
}

// return constant values for some particles
// in seconds
double Particle::LifeTime() const
{
	switch (std::abs(_pdg))
	{
		case 13:
			return 2.1969811e-6;	//muon
		case 211:
			return 2.6033e-8;	//pion
		case 111:
			return 8.52e-17;	//pion0
		default:
			return 1e26;		//lifetime of uni
	}
}

double Particle::LabLifeTime() const
{
	return LifeTime() * Gamma();
}

double Particle::LabSpace() const
{
	return LabLifeTime() * Beta();
}

double Particle::EKin() const
{
	return E() - M();
}

// overloading mass return function
double Particle::M() const
{
	return std::abs(TLorentzVector::M());
}

//////// non const functions
//
void Particle::SetPdg(int pdg)
{
	_pdg = pdg;
}

void Particle::ChargeID(int charge)
{
	_charge = charge;
}

void Particle::ChargeID()
{
	_charge = RealQ();
}

void Particle::ChargeMisID()
{
	_charge = -RealQ();
}

void Particle::SetFourVector(const TLorentzVector &vv)
{
	TLorentzVector::SetE(vv.E());
	TLorentzVector::SetPx(vv.Px());
	TLorentzVector::SetPy(vv.Py());
	TLorentzVector::SetPz(vv.Pz());
}

void Particle::SetFourVector(double E, double px, double py, double pz)
{
	TLorentzVector::SetE(E);
	TLorentzVector::SetPx(px);
	TLorentzVector::SetPy(py);
	TLorentzVector::SetPz(pz);
}

// set mass keeping same energy
void Particle::SetM(double mm)
{
	if (mm <= 0.) {	// set particle massless
		if (P() > 0.)	// if there is a momentum, keep it
			SetE(P());
		else	// redefine momentum
			TLorentzVector::SetPxPyPzE(0, 0, 1., 1.);
	}
	else if (E() <= mm)	// particle at rest
		TLorentzVector::SetE(std::sqrt(std::pow(P(), 2) + std::pow(mm, 2)));
	else	// maintain energy
		TLorentzVector::SetRho(std::sqrt(std::pow(E(), 2) - std::pow(mm, 2)));
}

// set energy keeping same mass
void Particle::SetE(double ee)
{
	if (ee <= M()) {
		TLorentzVector::SetE(M());
		if (P() > 0.)
			TLorentzVector::SetRho(0.);
	}
	else {
		if (P() > 0.)
			TLorentzVector::SetRho(std::sqrt(std::pow(ee, 2) - std::pow(M(), 2)));
		else
			TLorentzVector::SetPz(std::sqrt(std::pow(ee, 2) - std::pow(M(), 2)));
		TLorentzVector::SetE(ee);
	}
}

void Particle::SetEKin(double ee)
{
	TLorentzVector::SetE(ee + M());
}

void Particle::SetP(double pp)
{
	TLorentzVector::SetE(std::sqrt(std::pow(pp, 2) + M2()));
	TLorentzVector::SetRho(pp);
}
