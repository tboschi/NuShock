#include "Particle.h"

//manual ctor

Particle::Particle(int pdgCode, double Px, double Py, double Pz, double E, double X, double Y, double Z) :
	pdg(pdgCode)
{
	Init(Px, Py, Pz, E, X, Y, Z);
}

Particle::Particle(int pdgCode, TLorentzVector *vec, TVector3 *pos) :
	pdg(pdgCode)
{
	if (!pos)
		Init(vec->Px(), vec->Py(), vec->Pz(), vec->E(), 0, 0, 0);
	else
		Init(vec->Px(), vec->Py(), vec->Pz(), vec->E(), 
		     pos->X(), pos->Y(), pos->Z());
}

Particle::Particle(int pdgCode, TLorentzVector &Vector) :
	pdg(pdgCode)
{
	Init(Vector.Px(), Vector.Py(), Vector.Pz(), Vector.E(), 0, 0, 0);
}

Particle::Particle(int pdgCode, TLorentzVector &vec, TVector3 &pos) :
	pdg(pdgCode)
{
	Init(vec.Px(), vec.Py(), vec.Pz(), vec.E(), 
	     pos.X(), pos.Y(), pos.Z());
}

//copy ctor
Particle::Particle(const Particle &P) :
	pdg(P.pdg),
	particleVec(P.particleVec),
	particlePos(P.particlePos),
	trackIn(P.trackIn),
	trackBack(P.trackBack),
	trackOut(P.trackOut),
	kShower(P.kShower)
{
}

Particle& Particle::operator=(const Particle &P)
{
	pdg = P.pdg;

	particleVec = P.particleVec;
	particlePos = P.particlePos;

	trackIn   = P.trackIn;
	trackBack = P.trackBack;
	trackOut  = P.trackOut;

	kShower = P.kShower;
}

//destructor
Particle::~Particle()
{
}

void Particle::Init(double Px, double Py, double Pz, double E, double X, double Y, double Z)
{
	particleVec = TLorentzVector(Px, Py, Pz, E);
	particlePos = TVector3(X, Y, Z);

	if (Mass() < -1e8)
	{
		std::cout << "Caution! Irregular four-vector given to particle " << Pdg() << std::endl;
	}
	else

	trackIn   = -1.;
	trackBack = -1.;
	trackOut  = -1.;
	kShower   = false;
}

int Particle::Pdg() const
{
	return pdg;
}

int Particle::Charge() const
{
	switch ( abs(Pdg()) )
	{
		case 1:
		case 3:
		case 5:
		case 7:		//down quarks
			return -1 * ((Pdg() > 0) - (Pdg() < 0));
		case 2:
		case 4:
		case 6:
		case 8:		//up quarks
			return  2 * ((Pdg() > 0) - (Pdg() < 0));
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
			return -3 * ((Pdg() > 0) - (Pdg() < 0));
		default:	//positive particles
			return  3 * ((Pdg() > 0) - (Pdg() < 0));
			break;
	}
}

double Particle::LifeTime() const
{
	switch (Pdg())
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

bool Particle::IsShower() const
{
	return kShower;
}

///four momentum stuff
//
//
//
/*
TLorentzVector* Particle::FourVectorPtr()
{
	return particleVec;
}
*/

TLorentzVector Particle::FourVector() const
{
	return particleVec;
}

double Particle::Mass() const
{
	return particleVec.M();
}

double Particle::Energy() const
{
	return particleVec.E();
}

double Particle::EnergyKin() const
{
	return particleVec.E() - particleVec.M();
}

double Particle::Momentum() const
{
	return particleVec.P();
}

double Particle::Transverse() const
{
	return particleVec.Pt();
}

double Particle::MomentumX() const
{
	return particleVec.Px();
}

double Particle::MomentumY() const
{
	return particleVec.Py();
}

double Particle::MomentumZ() const
{
	return particleVec.Pz();
}

double Particle::Theta() const
{
	return particleVec.Theta();
}

double Particle::Phi() const
{
	return particleVec.Phi();
}

double Particle::Beta() const
{
	return particleVec.Beta();
}

double Particle::Gamma() const
{
	return particleVec.Gamma();
}

TVector3 Particle::Direction() const
{
	return particleVec.Vect();
}

TVector3 Particle::Position() const
{
	return particlePos;
}

double Particle::X() const
{
	return particlePos.X();
}

double Particle::Y() const
{
	return particlePos.Y();
}

double Particle::Z() const
{
	return particlePos.Z();
}

double Particle::Dist() const
{
	return particlePos.Mag2();
}

double Particle::TrackIn() const
{
	return trackIn;
}

double Particle::TrackBack() const
{
	return trackBack;
}

double Particle::TrackOut() const
{
	return trackOut;
}

double Particle::TrackTot() const
{
	return TrackIn() + TrackBack() + TrackOut();
}

//////// non const functions
//
void Particle::SetPdg(int X)
{
	pdg = X;
}

void Particle::SetFourVector(TLorentzVector &V)
{
	particleVec.SetE(V.E());
	particleVec.SetPx(V.Px());
	particleVec.SetPy(V.Py());
	particleVec.SetPz(V.Pz());
}

void Particle::SetFourVector(double Px, double Py, double Pz, double E)
{
	particleVec.SetE(E);
	particleVec.SetPx(Px);
	particleVec.SetPy(Py);
	particleVec.SetPz(Pz);
}

void Particle::SetMass(double dM)
{
	double dE = Energy();
	if (dE < dM)
		dE = dM;

	particleVec.SetE(dE);
	SetRho(sqrt(dE*dE - dM*dM));
}

void Particle::SetEnergy(double dE)
{
	double dM = Mass();
	if (dE < dM)
		dE = dM;


	particleVec.SetE(dE);
	SetRho(sqrt(dE*dE - dM*dM));
}

void Particle::SetEnergyKin(double dE)
{
	SetEnergy(dE + Mass());
}

void Particle::SetMomentum(double dP)
{
	double dM = Mass();
	particleVec.SetE(sqrt(dP*dP + dM*dM));
	SetRho(dP);
}

void Particle::SetRho(double dR)
{
	if (Momentum() == 0.0)
	{
		particleVec.SetPx(0);
		particleVec.SetPy(0);
		particleVec.SetPz(1);
	}

	particleVec.SetRho(dR);
}

void Particle::SetTheta(double Ang)
{
	particleVec.SetTheta(Ang);
}

void Particle::SetPhi(double Ang)
{
	particleVec.SetPhi(Ang);
}

void Particle::SetPosition(TVector3 &V)
{
	particlePos.SetX(V.X());
	particlePos.SetY(V.Y());
	particlePos.SetZ(V.Z());
}

void Particle::SetPosition(double X, double Y, double Z)
{
	particlePos.SetX(X);
	particlePos.SetY(Y);
	particlePos.SetZ(Z);
}

void Particle::SetX(double X)
{
	particlePos.SetX(X);
}

void Particle::SetY(double X)
{
	particlePos.SetY(X);
}

void Particle::SetZ(double X)
{
	particlePos.SetZ(X);
}

void Particle::SetTrackIn(double X)
{
	trackIn = X;
}

void Particle::SetTrackBack(double X)
{
	trackBack = X;
}

void Particle::SetTrackOut(double X)
{
	trackOut = X;
}

void Particle::SetShower(bool X)
{
	kShower = X;
}
