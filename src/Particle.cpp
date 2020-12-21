#include "physics/Particle.h"

Particle::Particle(int pdgCode, double E, double px, double py, double pz) :
	_pdg(pdgCode),
	_charge(0),
	TLorentzVector(px, py, pz, E)
{
}

Particle::Particle(int pdgCode, const TLorentzVector &fv) :
	_pdg(pdgCode),
	_charge(0),
	TLorentzVector(fv)
{
}

std::ostream & operator<<(std::ostream &os, const Particle &p) {
	return os << "<pdg: " << p._pdg << ", mass " << p.M()
		  << ", charge " << p._charge << ", vec (" << p.E() << ", "
		  << p.Px() << ", " << p.Py() << ", " << p.Pz() << ")>";
}


/*
//copy ctor
Particle::Particle(const Particle &p) :
	_pdg(p._pdg),
	TLorentzVector(p)
{
}
*/


int Particle::Pdg() const
{
	return _pdg;
}

int Particle::Charge() const
{
	return _charge;
}

int Particle::RealCharge() const
{
	int sign = (_pdg > 0) - (_pdg < 0);
	switch (std::abs(_pdg))
	{
		case 1:
		case 3:
		case 5:
		case 7:		//down quarks
			return -1 * sign;
		case 2:
		case 4:
		case 6:
		case 8:		//up quarks
			return  2 * sign;
		case 12:
		case 14:
		case 16:
		case 18:
		case 21:
		case 22:
		case 23:
		case 111:
		case 221:
		case 331:
		case 223:
		case 333:
		case 130:
		case 310:
		case 311:
		case 421:
		case 2112:
		case 2114:
		case 3122:
		case 3212:	//chargeless
			return 0;
		case 11:
		case 13:
		case 15:
		case 17:	//negative particles
			return -3 * sign;
		default:	//positive particles
			return  3 * sign;
			break;
	}
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
	_charge = RealCharge();
}

void Particle::ChargeMisID()
{
	_charge = -RealCharge();
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
	if (E() < mm) {
		TLorentzVector::SetE(mm);
		TLorentzVector::SetRho(0.);
	}
	else {
		TLorentzVector::SetRho(std::sqrt(std::pow(E(), 2) - std::pow(mm, 2)));
	}
}

// set energy keeping same mass
void Particle::SetE(double ee)
{
	if (ee < M()) {
		TLorentzVector::SetE(M());
		TLorentzVector::SetRho(0.);
	}
	else {
		TLorentzVector::SetE(ee);
		TLorentzVector::SetRho(std::sqrt(std::pow(E(), 2) - std::pow(M(), 2)));
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
