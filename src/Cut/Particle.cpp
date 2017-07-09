#include "Particle.h"

//ctor from gHep
Particle::Particle(genie::GHepParticle *Candidate, double PosX, double PosY, double PosZ)
{
	SetP4(*(Candidate->P4()));
	SetPdg(abs(Candidate->Pdg()));
	SetX(PosX);
	SetY(PosY);
	SetZ(PosZ);
}

//manual ctor
Particle::Particle(int PdgCode, TLorentzVector *Vector, double PosX, double PosY, double PosZ)
{
	SetP4(*Vector);
	SetPdg(abs(PdgCode));
	SetX(PosX);
	SetY(PosY);
	SetZ(PosZ);
}

//copy ctor
Particle::Particle(const Particle &P)
{
	SetP4(P.GetP4());
	SetPdg(abs(P.Pdg()));
	SetPosition(P.Position());
}

int Particle::Pdg()
{
	return iPdg;
}

TLorentzVector Particle::GetP4()
{
	return P4;
}

double Particle::M()
{
	return P4.M();
}

double Particle::E()
{
	return P4.E();
}

double Particle::Ekin()
{
	return P4.E() - P4.M();
}

double Particle::P()
{
	return P4.P();
}

double Particle::Px()
{
	return P4.Px();
}

double Particle::Py()
{
	return P4.Py();
}

double Particle::Pz()
{
	return P4.Pz();
}

double Particle::Theta()
{
	return P4.Theta();
}

double Particle::Phi()
{
	return P4.Phi();
}

TVector3 Particle::Position()
{
	return Pos;
}

double Particle::X()
{
	return Pos.X();
}

double Particle::Y()
{
	return Pos.Y();
}

double Particle::Z()
{
	return Pos.Z();
}

void Particle::SetPdg(int X)
{
	iPdg = X;
}

void Particle::SetP4(TLorentzVector &V)
{
	P4.SetE(V.E());
	P4.SetPx(V.Px());
	P4.SetPy(V.Py());
	P4.SetPz(V.Pz());
}

void Particle::SetE(double E)
{
	P4.SetRho(E*E - P4.M2());
	P4.SetE(E);
}

void Particle::SetMass(double M)
{
	P4.SetRho(P4.E()*P4.E() - M*M);
}

void Particle::SetTheta(double Ang)
{
	P4.SetTheta(Ang);
}

void Particle::SetPhi(double Ang)
{
	P4.SetPhi(Ang);
}

void Particle::SetPosition(TVector3 &V)
{
	Pos.SetX(V.X());
	Pos.SetY(V.Y());
	Pos.SetZ(V.Z());
}

void Particle::SetPosition(double X, double Y, double Z)
{
	Pos.SetX(X);
	Pos.SetY(Y);
	Pos.SetZ(Z);
}

void Particle::SetX(double X)
{
	Pos.SetX(X);
}

void Particle::SetY(double X)
{
	Pos.SetY(X);
}

void Particle::SetZ(double X)
{
	Pos.SetZ(X);
}

