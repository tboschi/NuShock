#include "Particle.h"

//ctor from gHep
Particle::Particle(genie::GHepParticle *Candidate, double PosX, double PosY, double PosZ)
{
	SetP4(*(Candidate->P4()));
	SetPdg(abs(Candidate->Pdg()));
	SetCharge(Candidate->Charge()/3);	//genie use e/3 charge
	SetX(PosX);
	SetY(PosY);
	SetZ(PosZ);
	SetTrack(-1);
}

//manual ctor
Particle::Particle(int PdgCode, double Charge, TLorentzVector &Vector, double PosX, double PosY, double PosZ)
{
	SetP4(Vector);
	SetPdg(abs(PdgCode));
	SetCharge(Charge);
	SetX(PosX);
	SetY(PosY);
	SetZ(PosZ);
	SetTrack(-1);
}

Particle::Particle(int PdgCode, double Charge, TLorentzVector &Vector, TVector3 &Position)
{
	SetP4(Vector);
	SetPdg(abs(PdgCode));
	SetCharge(Charge);
	SetPosition(Position);
	SetTrack(-1);
}

//copy ctor
Particle::Particle(const Particle &P)
{
	SetP4(P.Px(), P.Py(), P.Pz(), P.E());
	SetPdg(abs(P.Pdg()));
	SetCharge(P.Charge());
	SetPosition(P.X(), P.Y(), P.Z());
	SetTrack(-1);
}

int Particle::Pdg() const
{
	return iPdg;
}

int Particle::Charge() const
{
	return iCharge;
}

TLorentzVector Particle::GetP4() const
{
	return P4;
}

double Particle::M() const
{
	return P4.M();
}

double Particle::E() const
{
	return P4.E();
}

double Particle::Ekin() const
{
	return P4.E() - P4.M();
}

double Particle::P() const
{
	return P4.P();
}

double Particle::Pt() const
{
	return P4.Pt();
}

double Particle::Px() const
{
	return P4.Px();
}

double Particle::Py() const
{
	return P4.Py();
}

double Particle::Pz() const
{
	return P4.Pz();
}

double Particle::Theta() const
{
	return P4.Theta();
}

double Particle::Phi() const
{
	return P4.Phi();
}

TVector3 Particle::Position() const
{
	return Pos;
}

double Particle::X() const
{
	return Pos.X();
}

double Particle::Y() const
{
	return Pos.Y();
}

double Particle::Z() const
{
	return Pos.Z();
}

double Particle::Track() const
{
	return dTrack;
}

////////
void Particle::SetPdg(int X)
{
	iPdg = X;
}

void Particle::SetCharge(int X)
{
	iCharge = X;
}

void Particle::SetP4(TLorentzVector &V)
{
	P4.SetE(V.E());
	P4.SetPx(V.Px());
	P4.SetPy(V.Py());
	P4.SetPz(V.Pz());
}

void Particle::SetP4(double Px, double Py, double Pz, double E)
{
	P4.SetE(E);
	P4.SetPx(Px);
	P4.SetPy(Py);
	P4.SetPz(Pz);
}

void Particle::SetE(double E)
{
	if (P() > 0.0)
	{
		P4.SetRho(E*E - P4.M2());
		P4.SetE(E);
	}
	else P4.SetE(E);
}

void Particle::SetMass(double M)
{
	P4.SetRho(sqrt(P4.E()*P4.E() - M*M));
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

void Particle::SetTrack(double X)
{
	dTrack = X;
}
