#include "Particle.h"

//manual ctor
Particle::Particle(int PdgCode, TLorentzVector &Vector, double PosX, double PosY, double PosZ)
{
	SetPdg(abs(PdgCode));
	SetTau();	//in seconds

	SetP4(Vector);

	SetX(PosX);
	SetY(PosY);
	SetZ(PosZ);
	SetTrackIn(-1);
	SetTrackOut(-1);
	SetShower(false);
}

Particle::Particle(int PdgCode, TLorentzVector &Vector, TVector3 &Position)
{
	SetPdg(abs(PdgCode));
	SetTau();	//in seconds

	SetP4(Vector);

	SetPosition(Position);
	SetTrackIn(-1);
	SetTrackOut(-1);
	SetShower(false);
}

//copy ctor
Particle::Particle(const Particle &P)
{
	SetP4(P.Px(), P.Py(), P.Pz(), P.E());
	SetPdg(abs(P.Pdg()));
	SetTau();
	SetPosition(P.X(), P.Y(), P.Z());
	SetTrackIn(-1);
	SetTrackOut(-1);
	SetShower(false);
}

int Particle::Pdg() const
{
	return iPdg;
}

double Particle::Charge() const
{
	return iCharge/3.0;
}

bool Particle::IsShower() const
{
	return bShower;
}

double Particle::Tau() const
{
	return dTau;
}

double Particle::LabTau() const
{
	return Tau() * P4.Gamma();
}

double Particle::LabSpace() const
{
	return LabTau() * P4.Beta();
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

TVector3 Particle::Direction() const
{
	return GetP4().Vect();
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

double Particle::TrackIn() const
{
	return dTrackIn;
}

double Particle::TrackOut() const
{
	return dTrackOut;
}

double Particle::TrackTot() const
{
	return TrackIn() + TrackOut();
}

////////
void Particle::SetPdg(int X)
{
	iPdg = X;
}

void Particle::SetCharge()
{
	switch ((unsigned int) abs(Pdg()))
	{
		case 1:
		case 3:
		case 5:
		case 7:
			iCharge = -1*((Pdg() > 0) - (Pdg() < 0));
			break;
		case 2:
		case 4:
		case 6:
		case 8:
			iCharge = +2*((Pdg() > 0) - (Pdg() < 0));
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
		case 3212:
			iCharge = 0;
			break;
		case 11:
		case 13:
		case 15:
		case 17:
			iCharge = -3*((Pdg() > 0) - (Pdg() < 0));
		default:
			iCharge = 3*((Pdg() > 0) - (Pdg() < 0));
			break;
	}
}


void Particle::SetTau()
{
	if (Pdg() == 13)
		dTau = 2.1969811e-6;
	else if (Pdg() == 211)
		dTau = 2.6033e-8;
	else if (Pdg() == 111)
		dTau = 8.52e-17;
	else dTau = 1e26;
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

void Particle::SetEnergy(double dE)
{
	double dM = M();
	P4.SetE(dE);
	SetRho(sqrt(dE*dE - dM*dM));
}

void Particle::SetEnergyKin(double dE)
{
	SetEnergy(dE + M());
}

void Particle::SetMomentum(double dP)
{
	double dM = M();
	P4.SetE(sqrt(dP*dP + dM*dM));
	SetRho(dP);
}

void Particle::SetMass(double dM)
{
	double dE = E();
	SetRho(sqrt(dE*dE - dM*dM));
}

void Particle::SetRho(double dR)
{
	if (P() != 0.0)
		P4.SetRho(dR);
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

void Particle::SetTrackIn(double X)
{
	dTrackIn = X;
}

void Particle::SetTrackOut(double X)
{
	dTrackOut = X;
}

void Particle::SetShower(bool X)
{
	bShower = X;
}
