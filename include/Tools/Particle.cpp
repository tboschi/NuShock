#include "Tools/Particle.h"

//manual ctor

Particle::Particle() :
	iPdg(0),
	bShower(false)
{
	Init(0, 0, 0, 0, 0, 0, 0);
}

Particle::Particle(int PdgCode, double Px, double Py, double Pz, double E, double X, double Y, double Z) :
	iPdg(PdgCode)
{
	Init(Px, Py, Pz, E, X, Y, Z);
}

Particle::Particle(int PdgCode, TLorentzVector *Vector, TVector3 *Position) :
	iPdg(PdgCode),
	dTrackIn(-1),
	dTrackBack(-1),
	dTrackOut(-1),
	bShower(false)
{
	Init(Vector->Px(), Vector->Py(), Vector->Pz(), Vector->E(), 
	     Position->X(), Position->Y(), Position->Z());
}

Particle::Particle(int PdgCode, TLorentzVector *Vector) :
	iPdg(PdgCode),
	dTrackIn(-1),
	dTrackBack(-1),
	dTrackOut(-1),
	bShower(false)
{
	Init(Vector->Px(), Vector->Py(), Vector->Pz(), Vector->E(), 0, 0, 0);
}

Particle::Particle(int PdgCode, TLorentzVector &Vector, TVector3 &Position) :
	iPdg(PdgCode),
	dTrackIn(-1),
	dTrackBack(-1),
	dTrackOut(-1),
	bShower(false)
{
	Init(Vector.Px(), Vector.Py(), Vector.Pz(), Vector.E(), 
	     Position.X(), Position.Y(), Position.Z());
}

Particle::Particle(int PdgCode, TLorentzVector &Vector) :
	iPdg(PdgCode),
	dTrackIn(-1),
	dTrackBack(-1),
	dTrackOut(-1),
	bShower(false)
{
	Init(Vector.Px(), Vector.Py(), Vector.Pz(), Vector.E(), 0, 0, 0);
}

//copy ctor
Particle::Particle(const Particle &P) :
	iPdg(P.iPdg),
	//ParticleVec(new TLorentzVector(P.FourVector())),
	//ParticlePos(new TVector3(P.Position())),
	ParticleVec(P.ParticleVec),
	ParticlePos(P.ParticlePos),
	dTrackIn(P.dTrackIn),
	dTrackBack(P.dTrackBack),
	dTrackOut(P.dTrackOut),
	bShower(P.bShower)
{
}

//destructor
Particle::~Particle()
{
	//std::cout << "P delete " << ParticleVec << std::endl;
	//delete ParticleVec; 
	//delete ParticlePos; 
	//ParticleVec = 0;
	//ParticlePos = 0;
}

void Particle::Init(double Px, double Py, double Pz, double E, double X, double Y, double Z)
{
	//ParticleVec = new TLorentzVector(Px, Py, Pz, E);
	//ParticlePos = new TVector3(X, Y, Z);
	ParticleVec = TLorentzVector(Px, Py, Pz, E);
	ParticlePos = TVector3(X, Y, Z);
}

int Particle::Pdg() //const
{
	return iPdg;
}

int Particle::Charge()
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

double Particle::LifeTime()
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

double Particle::LabLifeTime() 
{
	return LifeTime() * Gamma();
}

double Particle::LabSpace()
{
	return LabLifeTime() * Beta();
}

bool Particle::IsShower() //const
{
	return bShower;
}

///four momentum stuff
//
//
//
/*
TLorentzVector* Particle::FourVectorPtr()
{
	return ParticleVec;
}
*/

TLorentzVector Particle::FourVector() //const
{
	return ParticleVec;
}

double Particle::Mass() 
{
	return ParticleVec.M();
}

double Particle::Energy()
{
	return ParticleVec.E();
}

double Particle::EnergyKin()
{
	return ParticleVec.E() - ParticleVec.M();
}

double Particle::Momentum()
{
	return ParticleVec.P();
}

double Particle::Transverse()
{
	return ParticleVec.Pt();
}

double Particle::MomentumX()
{
	return ParticleVec.Px();
}

double Particle::MomentumY()
{
	return ParticleVec.Py();
}

double Particle::MomentumZ()
{
	return ParticleVec.Pz();
}

double Particle::Theta()
{
	return ParticleVec.Theta();
}

double Particle::Phi()
{
	return ParticleVec.Phi();
}

double Particle::Beta()
{
	return ParticleVec.Beta();
}

double Particle::Gamma()
{
	return ParticleVec.Gamma();
}

TVector3 Particle::Direction()
{
	return ParticleVec.Vect();
}

TVector3 Particle::Position() //const
{
	return ParticlePos;
}

double Particle::X()
{
	return ParticlePos.X();
}

double Particle::Y()
{
	return ParticlePos.Y();
}

double Particle::Z()
{
	return ParticlePos.Z();
}

double Particle::Dist()
{
	return ParticlePos.Mag2();
}

double Particle::TrackIn() //const
{
	return dTrackIn;
}

double Particle::TrackBack() //const
{
	return dTrackBack;
}

double Particle::TrackOut() //const
{
	return dTrackOut;
}

double Particle::TrackTot() //const
{
	return TrackIn() + TrackBack() + TrackOut();
}

////////
void Particle::SetPdg(int X)
{
	iPdg = X;
}

void Particle::SetFourVector(TLorentzVector &V)
{
	ParticleVec.SetE(V.E());
	ParticleVec.SetPx(V.Px());
	ParticleVec.SetPy(V.Py());
	ParticleVec.SetPz(V.Pz());
}

void Particle::SetFourVector(double Px, double Py, double Pz, double E)
{
	ParticleVec.SetE(E);
	ParticleVec.SetPx(Px);
	ParticleVec.SetPy(Py);
	ParticleVec.SetPz(Pz);
}

void Particle::SetMass(double dM)
{
	double dE = Energy();
	if (dE < dM)
		dE = dM;

	ParticleVec.SetE(dE);
	SetRho(sqrt(dE*dE - dM*dM));
}

void Particle::SetEnergy(double dE)
{
	double dM = Mass();
	if (dE < dM)
		dE = dM;


	ParticleVec.SetE(dE);
	SetRho(sqrt(dE*dE - dM*dM));
}

void Particle::SetEnergyKin(double dE)
{
	SetEnergy(dE + Mass());
}

void Particle::SetMomentum(double dP)
{
	double dM = Mass();
	ParticleVec.SetE(sqrt(dP*dP + dM*dM));
	SetRho(dP);
}

void Particle::SetRho(double dR)
{
	if (Momentum() == 0.0)
	{
		ParticleVec.SetPx(0);
		ParticleVec.SetPy(0);
		ParticleVec.SetPz(1);
	}

	ParticleVec.SetRho(dR);
}

void Particle::SetTheta(double Ang)
{
	ParticleVec.SetTheta(Ang);
}

void Particle::SetPhi(double Ang)
{
	ParticleVec.SetPhi(Ang);
}

void Particle::SetPosition(TVector3 &V)
{
	ParticlePos.SetX(V.X());
	ParticlePos.SetY(V.Y());
	ParticlePos.SetZ(V.Z());
}

void Particle::SetPosition(double X, double Y, double Z)
{
	ParticlePos.SetX(X);
	ParticlePos.SetY(Y);
	ParticlePos.SetZ(Z);
}

void Particle::SetX(double X)
{
	ParticlePos.SetX(X);
}

void Particle::SetY(double X)
{
	ParticlePos.SetY(X);
}

void Particle::SetZ(double X)
{
	ParticlePos.SetZ(X);
}

void Particle::SetTrackIn(double X)
{
	dTrackIn = X;
}

void Particle::SetTrackBack(double X)
{
	dTrackBack = X;
}

void Particle::SetTrackOut(double X)
{
	dTrackOut = X;
}

void Particle::SetShower(bool X)
{
	bShower = X;
}
