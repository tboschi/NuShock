#include "Detector/Tracker.h"

Tracker::Tracker(std::string ConfigName) : Detector(ConfigName)
{
}

void Tracker::TrackProcess(Particle *&P)
{
	if (P->TrackIn() < 0)
		TrackLength(P);
	TrackSmearing(P);

	if (!IsDetectable(P))
	{
		delete P;
		P = 0;
	}
}

void Tracker::TrackSmearing(Particle *P)
{
	double iM = P->Mass();
	double iEkin = P->EnergyKin();
	double iP = P->Momentum();
	double iTheta = P->Theta();
	double iPhi = P->Phi();
	double SigmaE, SigmaP, SigmaA; 

	if (P->Pdg() == 11 || P->Pdg() == 22)
	{
		SigmaA = Get("Angle_Gamma") / Const::fDeg;
		double StatE = Get("Energ_Gamma") / sqrt(iEkin);
		double SystE = Get("Ebias_Gamma");
		SigmaE = sqrt(pow(StatE, 2)+pow(SystE, 2)) * iEkin;

		P->SetEnergyKin(GenMT->Gaus(iEkin, SigmaE));
		P->SetTheta(GenMT->Gaus(iTheta, SigmaA));
		P->SetPhi(GenMT->Gaus(iPhi, SigmaA));
	}
	else if (P->Pdg() == 13)
	{
		SigmaA = Get("Angle_Muon") / Const::fDeg;
		double Ratio = P->TrackIn()/(P->TrackIn()+P->TrackOut());
		//std::cout << "13  " << Ratio << "\t" << P->TrackIn() << "\t" << P->TrackOut() << std::endl;
		if (Ratio > Get("Containment"))	//90% of track inside
		{
			double StatP = Get("Range_Muon") / iP;
			SigmaP = StatP * iP;
		}
		else
		{
			double StatP = Get("Exiti_Muon");
			SigmaP = StatP * iP;
		}


		P->SetMomentum(GenMT->Gaus(iP, SigmaP));
		P->SetTheta(GenMT->Gaus(iTheta, SigmaA));
		P->SetPhi(GenMT->Gaus(iPhi, SigmaA));
	}
	else if (P->Pdg() == 211)
	{
		SigmaA = Get("Angle_Pion") / Const::fDeg;
		double Ratio = P->TrackIn()/(P->TrackIn()+P->TrackOut());
		//std::cout << "211 " << Ratio << "\t" << P->TrackIn() << "\t" << P->TrackOut() << std::endl;
		if (!P->IsShower() && Ratio > Get("Containment"))	//pion track, not shower
		{
			double StatP = Get("Range_Pion") / iP;
			SigmaP = StatP * iP;
		}
		else
		{
			double StatP = Get("Exiti_Pion");
			SigmaP = StatP * iP;
		}

		P->SetMomentum(GenMT->Gaus(iP, SigmaP));
		P->SetTheta(GenMT->Gaus(iTheta, SigmaA));
		P->SetPhi(GenMT->Gaus(iPhi, SigmaA));
	}
	else if (P->Charge() != 0)	//other and protons?
	{
		SigmaA = Get("Angle_Hadron") / Const::fDeg;
		double StatE = Get("Energ_Hadron") / sqrt(iEkin);
		double SystE = Get("Ebias_Hadron");
		SigmaE = sqrt(pow(StatE, 2)+pow(SystE, 2)) * iEkin;

		P->SetEnergyKin(GenMT->Gaus(iEkin, SigmaE));
		P->SetTheta(GenMT->Gaus(iTheta, SigmaA));
		P->SetPhi(GenMT->Gaus(iPhi, SigmaA));
	}

	P->SetTrackIn(GenMT->Gaus(P->TrackIn(), Get("Vertex")));
	if (P->TrackIn() < 0)
		P->SetTrackIn(0);
	P->SetTrackOut(GenMT->Gaus(P->TrackOut(), Get("Vertex")));
	if (P->TrackOut() < 0)
		P->SetTrackOut(0);
}

void Tracker::TrackLength(Particle *P)	//This should not change *P
{
	switch (abs(P->Pdg()))
	{
		case 11: {
			double tmax = log(P->Energy()/CriticalEnergy()) - 1.0;
			double Depth = RadiationLength() * GenMT->Gaus(tmax, 0.3*tmax);	//justify this
			P->SetTrackIn(Depth);
			break; }
		case 22: {
			P->SetTrackIn(GammaDecay());
			break; }
		case 13: {
			double iE = P->Energy();
			double dE;
			double dStep = 0.01;	//cm step
			double dx, LengthIn = 0, LengthOut = 0;
			
			TVector3 Step(P->Direction().Unit());
			TVector3 Pos(P->Position());
			TVector3 Start(P->Position());


			bool InOut;
			while (IsDetectable(P) && !IsDecayed(P, dStep))		//this should quit when particle decays too!
			{
				dx = GenMT->Gaus(dStep, Get("Vertex"));

				dE = dx * EnergyLoss(P, InOut);	//inside material

				if (InOut)
					LengthIn += dx;			//even b2b det?
				else
					LengthOut += dx;

				Pos += Step*dx;		//move by dx
				P->SetEnergy(P->Energy() - dE);
				P->SetPosition(Pos);
			}

			P->SetTrackIn(LengthIn);
			P->SetTrackOut(LengthOut);
			P->SetEnergy(iE);		//reset to original energy
			P->SetPosition(Start);		//reset to original position
			break; }
		case 211: {
			double iE = P->Energy();
			double iM = P->Mass();
			double dE;
			double dStep = 0.01;	//cm step
			double dx, LengthIn = 0, LengthOut = 0;
			
			TVector3 Step(P->Direction().Unit());
			TVector3 Pos(P->Position());
			TVector3 Start(P->Position());

			double TotTrack = 0;
			int Layer = 0;
			bool InOut;
			while (IsDetectable(P) && !IsDecayed(P, dStep))		//this should quit when particle decays too!
			{
				dx = GenMT->Gaus(dStep, Get("Vertex"));
				double Length = GenMT->Exp(RadiationLength(1));
				double Cover = 0;

				while (IsDetectable(P) && !IsDecayed(P, dStep) && Cover < Length)	//this should quit when particle decays too!
				{
					dx = GenMT->Gaus(dStep, Get("Vertex"));

					dE = dx * EnergyLoss(P, InOut);	//inside material
					Cover += dx;

					if (InOut)
						LengthIn += dx;			//even b2b det?
					else
						LengthOut += dx;

					Pos += Step*dx;		//move by dx
					P->SetEnergy(P->Energy() - dE);
					P->SetPosition(Pos);
				}

				int Mult = ceil(pow(P->EnergyKin() * 1e6, 0.18) * 0.15);
				if (Mult > 1)
				{
					++Layer;
					P->SetEnergyKin(P->EnergyKin()/Mult);
				}
			}

			if (GenMT->Rndm() > 1.0/Layer)		//need something better
				P->SetShower(true);
			else P->SetShower(false);

			P->SetTrackIn(LengthIn);
			P->SetTrackOut(LengthOut);
			P->SetEnergy(iE);			//reset to original energy
			P->SetPosition(Start);		//reset to original position
			break; }
		default:
			break;
	}
}

double Tracker::GammaDecay()
{
	double Path = 9.0/7.0 * RadiationLength(0);
	return GenMT->Exp(Path);
}

double Tracker::CriticalEnergy()	//assuming same for positron and electron
{
	switch (GetMaterial("InTarget"))
	{
		case LAr:
			return 0.03284;	//GeV
		case GasAr:
			return 0.03803;	//GeV
		case Fe:
			return 0.02168;	//GeV
		default:
			return 0;
	}
}

double Tracker::RadiationLength(bool Nuclear)
{
	if (!Nuclear)
		switch (GetMaterial("InTarget"))
		{
			case 1:
				return 19.55/1.3945 / 100;
			case 2:
				return 19.55/0.1020 / 100;
			case 3:
				return 13.84/7.874 / 100;
			default:
				return 0;
		}
	else
		switch (GetMaterial("InTarget"))
		{
			case 1:
				return 119.7/1.3945 / 100;
			case 2:
				return 119.7/0.1020 / 100;
			case 3:
				return 132.1/7.874 / 100;
			default:
				return 0;
		}
}

double Tracker::EnergyLoss(Particle *P, bool &Contained)
{
	if (IsContained(P) && IsInside(P))
	{
		Contained = true;
		return BetheLoss(P, GetMaterial("InTarget"));
	}
	else if (IsContained(P) && GetMaterial("BackTarget") != 0)
	{
		Contained = true;
		return BetheLoss(P, GetMaterial("BackTarget"));
	}
	else
	{
		Contained = false;
		return BetheLoss(P, GetMaterial("OutTarget"));
	}
}

double Tracker::BetheLoss(Particle *P, Material Target)
{
	switch (Target)
	{
		case 1:
			return Bethe(P, 1.3945, 188.0, 18, 40);	//LAr, code 1
		case 2:
			return Bethe(P, 0.1020, 188.0, 18, 40);	//GasAr, code 2
		case 3:
			return Bethe(P, 7.874, 286.0, 26, 56);	//Fe, code 3
		default:
			0;
	}
}

double Tracker::Bethe(Particle *P, double Density, double I, int Z, int A)
{
	double K = 0.307075;		//From PDG MeV mol-1 cm2
	double e2M = Const::fMElectron / P->Mass();
	double Beta2 = pow(P->Beta(), 2);
	double Gamma2 = pow(P->Gamma(), 2);

	double Wmax = (2000 * Const::fMElectron * Beta2 * Gamma2) / 
		(1 + 2 * P->Gamma() * e2M + e2M * e2M);
	double LogArg = 2000 * Const::fMElectron * Beta2 * Gamma2 * Wmax /
		(1e-12 * I * I); 	//Everything in MeV

	return 0.1 * Density * (K * Z) / (A * Beta2 * (0.5 * log (LogArg) - Beta2));	//stopping power in GeV/m (0.1*)
}

bool Tracker::IsDecayed(Particle *P, double dx)	//Threshold check
{
	return GenMT->Rndm() > exp(-dx/(Const::fC * P->LabSpace()));	//this should quit when particle decays too!
}

bool Tracker::IsDetectable(Particle *P)	//Threshold check
{
	double Threshold = 0.0;

	switch (abs(P->Pdg()))
	{
		case 11:
		case 22:
			Threshold = Get("Thres_Gamma");
			break;
		case 13:
			Threshold = Get("Thres_Muon");
			break;
		case 211:
			Threshold = Get("Thres_Pion");
			break;
		default:
			if (P->Charge() != 0)
				Threshold = Get("Thres_Hadron");
			else
				Threshold = 2 * P->EnergyKin();		//not detectable
			break;
	}

	return (P->EnergyKin() > Threshold);
}

void Tracker::Pi0Decay(Particle *Pi0, Particle *&PA, Particle *&PB)
{
	//in rest frame
	double M_Pion0 = Const::fMPion0;
	TLorentzVector GammaA(0, 0, M_Pion0/2.0, M_Pion0/2.0); 
	TLorentzVector GammaB(0, 0, -M_Pion0/2.0, M_Pion0/2.0); 

	TVector3 vBoost(Pi0->FourVector().BoostVector());
	TVector3 Start(Pi0->Position());		//starting point is vertex
	double Theta = GenMT->Uniform(-Const::fPi, Const::fPi);
	double Phi = GenMT->Uniform(-Const::fPi, Const::fPi);

	GammaA.SetTheta(Theta);
	GammaB.SetTheta(Theta + Const::fPi);
	GammaA.SetPhi(Phi);
	GammaB.SetPhi(Phi + Const::fPi);

	GammaA.Boost(vBoost);
	GammaB.Boost(vBoost);

	TVector3 MoveA(GammaA.Vect().Unit());
	TVector3 MoveB(GammaB.Vect().Unit());
	MoveA *= GammaDecay();
	MoveB *= GammaDecay();
	MoveA += Start;
	MoveB += Start;

	delete PA, PB;
	PA = 0, PB = 0;

	PA = new Particle(22, GammaA, MoveA);	//here are the photons
	PB = new Particle(22, GammaB, MoveB);	//position should be different
}

void Tracker::Focus(Particle *P)
{
	double Radius = sqrt(pow(Get("Width"), 2) + pow(Get("Height"), 2));
	double SigmaT = atan2(Radius, Get("Baseline"));					//3sigma will be inside the detector (better distribution needed)

	P->SetTheta( abs(GenMT->Gaus(0, SigmaT)) );	
	P->SetPhi( GenMT->Uniform(-Const::fPi, Const::fPi) );
}

bool Tracker::IsContained(Particle *P)
{
	if (P->X() > Xstart() && P->X() < Xend() &&
	    P->Y() > Ystart() && P->Y() < Yend() &&
	    P->Z() > Zstart() )
		return true;
	else return false;
}

bool Tracker::IsInside(Particle *P)
{
	if (P->X() > Xstart() && P->X() < Xend() &&
	    P->Y() > Ystart() && P->Y() < Yend() &&
	    P->Z() > Zstart() && P->Z() < Zend() )
		return true;
	else return false;
}
