#include "Tracker.h"

Tracker::Tracker(std::string configName, std::string mod) :
	Detector(configName, mod)
{
}

bool Tracker::Reconstruct(Particle &P)
{
	if (P.Pdg() && IsDetectable(P))	//valid particle
	{
		//if (P.Dist() == 0.0)
		//	Vertex(P);

		if (P.TrackOut() < 0)
			Length(P);

		//with approximation
		double track2 = pow(P.TrackFGT() * cos(P.Theta()), 2) +
				pow(P.TrackFGT() * P.Theta() * P.Phi(), 2);
		double sag = 0.1 * std::abs(P.RealCharge()) * Get("MagneticField") *
			     track2 / P.Momentum();

		//without approximation
		//double rad = 0.1 * std::abs(P.RealCharge()) *
		//	     Get("MagneticField") / P.Momentum();
		//double track = P.TrackFGT() * sqrt(pow(cos(P.Theta()), 2) +
		//				   pow(P.Theta() * P.Phi(), 2));
		//double sag = rad * ( 1 - cos(track/rad) );

		P.SetSagitta(sag);
		//if SagittaRes <= 0 then this is false
		if (sag > Get("SagittaRes"))
			P.ChargeID();
		else
			P.ChargeID(0);

		Smearing(P);

		if (IsDetectable(P))
			return true;
		else
			return false;
	}
	else
		return false;
}

void Tracker::Vertex(Particle &P)
{
	double Z = GenMT->Uniform(Zstart(), Zend());
	double X, Y;
	if (Z > ZstartLAr() && Z < ZendLAr())
	{
		X = GenMT->Uniform(XstartLAr(), XendLAr());
		Y = GenMT->Uniform(YstartLAr(), YendLAr());
	}
	else if (Z > ZstartFGT() && Z < ZendFGT())
	{
		X = GenMT->Uniform(XstartFGT(), XendFGT());
		Y = GenMT->Uniform(YstartFGT(), YendFGT());
	}

	P.SetPosition(X, Y, Z);
}

void Tracker::Smearing(Particle &P)
{
	double iM = P.Mass();
	double iEkin = P.EnergyKin();
	double iP = P.Momentum();
	double iTheta = P.Theta();
	double iPhi = P.Phi();
	double SigmaE, SigmaP, SigmaA; 
	double StatE, SystE;
	double Ratio;
	bool hiRes;

	switch (std::abs(P.Pdg()))
	{
		case 11:
		case 22:
			SigmaA = Get("Angle_Gamma") / Const::Deg;
			StatE = Get("Energ_Gamma") / sqrt(iEkin);
			SystE = Get("Ebias_Gamma");
			SigmaE = sqrt(pow(StatE, 2)+pow(SystE, 2)) * iEkin;

			P.SetEnergyKin(GenMT->Gaus(iEkin, SigmaE));
			P.SetTheta(GenMT->Gaus(iTheta, SigmaA));
			P.SetPhi(GenMT->Gaus(iPhi, SigmaA));

			break;
		case 13:
			SigmaA = Get("Angle_Muon") / Const::Deg;
			Ratio = P.TrackTot() ? P.TrackIn()/P.TrackTot() : 0.0;

			//if (Ratio > Get("Containment"))	//80% of track inside
			hiRes = false;
			if ( IsInsideLAr(P) && P.TrackLAr() < ZsizeLAr()
					    && Ratio > Get("Containment") )
				hiRes = true;	//the particle is contained in LAr
			else if (P.TrackIn() > ZsizeFGT())
				hiRes = true;

			if (hiRes)
			{
				//std::cout << "muon in range, ";
				SigmaP = Get("Range_Muon") * iP;
				//double StatP = Get("Range_Muon") / iP;
				//SigmaP = StatP * iP;
			}
			else
			{
				//std::cout << "muon escaping, ";
				SigmaP = Get("Exiti_Muon") * iP;
				//double StatP = Get("Exiti_Muon");
				//SigmaP = StatP * iP;
			}
			//std::cout << P.TrackIn() << " - " << P.TrackOut() << " - " << P.TrackTot() << "; " << Ratio << std::endl;

			P.SetMomentum(GenMT->Gaus(iP, SigmaP));
			P.SetTheta(GenMT->Gaus(iTheta, SigmaA));
			P.SetPhi(GenMT->Gaus(iPhi, SigmaA));

			break;
		case 211:
			SigmaA = Get("Angle_Pion") / Const::Deg;
			Ratio = P.TrackIn()/P.TrackTot();

			hiRes = false;
			if ( IsInsideLAr(P) && P.TrackIn() < ZsizeLAr() && Ratio > Get("Containment") )
				hiRes = true;	//the particle is contained in LAr
			else if (P.TrackIn() > ZsizeFGT())
				hiRes = true;

			//std::cout << "HiRes " << std::boolalpha << hiRes << std::endl;
			if (!P.IsShower() && hiRes)	//pion track, not shower
			{
				//std::cout << "pion in range, ";
				SigmaP = Get("Range_Pion") * iP;
				//double StatP = Get("Range_Pion") / iP;
				//SigmaP = StatP * iP;
			}
			else
			{
				//std::cout << "pion escaping, ";
				SigmaP = Get("Exiti_Pion") * iP;
				//double StatP = Get("Exiti_Pion");
				//SigmaP = StatP * iP;
			}
			//std::cout << P.TrackIn() << " - " << P.TrackOut() << " - " << P.TrackTot() << "; " << Ratio << std::endl;

			P.SetMomentum(GenMT->Gaus(iP, SigmaP));
			P.SetTheta(GenMT->Gaus(iTheta, SigmaA));
			P.SetPhi(GenMT->Gaus(iPhi, SigmaA));

			break;
		default:	//I don't care about other particles
			//if (std::abs(P.Pdg()) == 2112 || P.Charge() != 0)	//other and protons?
			//{
			//	SigmaA = Get("Angle_Hadron") / Const::Deg;
			//	StatE = Get("Energ_Hadron") / sqrt(iEkin);
			//	SystE = Get("Ebias_Hadron");
			//	SigmaE = sqrt(pow(StatE, 2)+pow(SystE, 2)) * iEkin;
			//
			//	P.SetEnergyKin(GenMT->Gaus(iEkin, SigmaE));
			//	P.SetTheta(GenMT->Gaus(iTheta, SigmaA));
			//	P.SetPhi(GenMT->Gaus(iPhi, SigmaA));
			//}
			break;
	}

	//P.SetTrackIn(GenMT->Gaus(P.TrackIn(), Get("Vertex")));

	//if (P.TrackIn() < 0)
	//	P.SetTrackIn(0);

	//P.SetTrackOut(GenMT->Gaus(P.TrackOut(), Get("Vertex")));

	//if (P.TrackOut() < 0)
	//	P.SetTrackOut(0);
}

void Tracker::Length(Particle &P)	//This should not change P
{
	double tmax, Depth;
	double iE, iM, dE;
	//double step = 0.05 + , dx;	//step of 5cm
	double step = GenMT->Gaus(0.01, Get("Vertex"));	//step of 1cm
	double LengthLAr = 0.0, LengthFGT = 0.0, LengthOut = 0.0;
	double length = -1, cover = -2;
	int layer = 0;

	//TVector3 Step(P.Direction().Unit());
	//TVector3 Pos(P.Position());
	//TVector3 Start(P.Position());

	Particle clone = P;
	//TVector3 pos   = P.Position();
	//TVector3 start = P.Position();
	TVector3 dir  = P.Direction().Unit();

	switch (std::abs(P.Pdg()))
	{
		case 11:
			tmax = log(P.Energy()/CriticalEnergy(P)) - 1.0;
			Depth = RadiationLength(P) * GenMT->Gaus(tmax, 0.3*tmax);

			if (IsInsideLAr(P))
			{
				P.SetTrackIn(Depth, 1);
				P.SetTrackIn(0, 0);
			}
			else
			{
				P.SetTrackIn(Depth, 0);
				P.SetTrackIn(0, 1);
			}

			P.SetTrackOut(0);
			break;
		case 22:
			if (IsInsideLAr(P))
			{
				P.SetTrackIn(GammaDecay(P), 1);
				P.SetTrackIn(0, 0);
			}
			else
			{
				P.SetTrackIn(GammaDecay(P), 0);
				P.SetTrackIn(0, 1);
			}

			P.SetTrackOut(0);
			break;
		case 13:
			//std::cout << "position " << clone.X() << " - " << clone.Y() << " - " << clone.Z() << std::endl;
			while (IsDetectable(clone) && !IsDecayed(clone, LengthLAr+LengthFGT+LengthOut))
				//quit when particle decays too!
			{
				int InOut;
				//double dx = GenMT->Gaus(step, Get("Vertex"));
				double dx = step;
				double loss = EnergyLoss(clone, InOut);
				double dE = dx * loss;

				//std::cout << "Energy " << dE << "(" << InOut << ")" << std::endl;
				if (InOut == 1)	//inside detector
					LengthLAr += dx;
				else if (InOut == 2)	//inside detector
					LengthFGT += dx;
				else		//outside detector
					LengthOut += dx;


				TVector3 pos = clone.Position() + dx * dir;	//move by dx
				clone.SetEnergy(clone.Energy() - dE);
				clone.SetPosition(pos);
			}
			//std::cout << "position " << clone.X() << " - " << clone.Y() << " - " << clone.Z() << std::endl;

			P.SetTrackIn(LengthLAr, 1);
			P.SetTrackIn(LengthFGT, 0);
			P.SetTrackOut(LengthOut);
			//std::cout << "Setto " << LengthLAr << ", " << LengthFGT << ", " << LengthOut << std::endl;
			break;
		case 211:
			while (IsDetectable(clone) && !IsDecayed(clone, LengthLAr+LengthFGT+LengthOut))
			{
				if (cover >= length)	//distance covered more than interaction length
				{
					int mult = ceil(pow(clone.EnergyKin() * 1e6, 0.18) * 0.15);
					if (mult > 1)	//multiplicity of hadronic interaction
					{
						++layer;	//number of interactions
						//clone.SetEnergyKin(clone.EnergyKin()/mult);
						//divide energy by number (mult) of duaghter particles
					}
					length = -1;	//repeat
				}

				if (length < 0)
				{
					length = GenMT->Exp(RadiationLength(clone, 1));
					cover = 0;
				}

				int InOut;
				//double dx = GenMT->Gaus(step, Get("Vertex"));
				double dx = step;
				double dE = dx * EnergyLoss(clone, InOut);	//inside material

				cover += dx;	//total distance covered in one layer

				if (InOut == 1)	//inside detector
					LengthLAr += dx;
				else if (InOut == 2)	//inside detector
					LengthFGT += dx;
				else		//outside detector
					LengthOut += dx;

				TVector3 pos = clone.Position() + dx * dir;	//move by dx
				clone.SetEnergy(clone.Energy() - dE);
				clone.SetPosition(pos);
			}

			if (GenMT->Rndm() > 1.0/layer)		//need something better
				P.SetShower(true);
			else
				P.SetShower(false);

			P.SetTrackIn(LengthLAr, 1);
			P.SetTrackIn(LengthFGT, 0);
			P.SetTrackOut(LengthOut);
			break;
		default:
			break;
	}
}

double Tracker::GammaDecay(const Particle &P)
{
	double Path = 9.0/7.0 * RadiationLength(P, 0);
	return GenMT->Exp(Path);
}

double Tracker::CriticalEnergy(const Particle &P)
{
	Detector::Material Element;

	if (IsInsideLAr(P))
		Element = GetMaterial("TargetLAr");
	else if (IsInsideFGT(P))
		Element = GetMaterial("TargetFGT");
	else
		Element = GetMaterial("TargetOut");

	return CriticalEnergy(Element);
}


double Tracker::RadiationLength(const Particle &P, bool Nuclear)
{
	Detector::Material Element;

	if (IsInsideLAr(P))
		Element = GetMaterial("TargetLAr"), Nuclear;
	else if (IsInsideFGT(P))
		Element = GetMaterial("TargetFGT"), Nuclear;
	else
		Element = GetMaterial("TargetOut"), Nuclear;

	return RadiationLength(Element, Nuclear);
}

double Tracker::EnergyLoss(const Particle &P, int &inout)
{
	if (IsInsideLAr(P))
	{
		inout = 1;
		return BetheLoss(P, GetMaterial("TargetLAr"));
	}
	else if (IsInsideFGT(P))
	{
		inout = 2;
		return BetheLoss(P, GetMaterial("TargetFGT"));
	}
	else
	{
		inout = 0;
		return BetheLoss(P, GetMaterial("TargetOut"));
	}
}

double Tracker::BetheLoss(const Particle &P, Material Target)
{
	switch (Target)
	{
		case LAr:
			return Bethe(P, 1.3945, 188.0, 18, 40);
		case GasAr:
			//return Bethe(P, 0.1020, 188.0, 18, 40);
			return Bethe(P, 0.024, 188.0, 18, 40);	//10 atm, 293 K
		case Fe:
			return Bethe(P, 7.875, 286.0, 26, 56);
		case Pb:
			return Bethe(P, 11.34, 823.0, 82, 207);
		default:
			0;
	}
}

// Bethe formula is missing charge of particle
							//I is in eV
double Tracker::Bethe(const Particle &P, double Density, double I, int Z, int A)	//GeV/m
{
	double K = 0.307075;		//From PDG MeV mol-1 cm2
	double e2M = Const::MElectron / P.Mass();
	double Beta2  = pow(P.Beta(), 2);
	double Gamma2 = pow(P.Gamma(), 2);
			//electron mass in MeV -> *1000
	double Wmax = (2000 * Const::MElectron * Beta2 * Gamma2) / 
		(1 + 2 * P.Gamma() * e2M + e2M * e2M);

	double LogArg = 2000 * Const::MElectron * Beta2 * Gamma2 * Wmax /
		(1e-12 * I * I); 	//Everything in MeV

	return 0.1 * Density * K * Z / (A * Beta2) * (0.5 * log (LogArg) - Beta2);	//stopping power in GeV/m (0.1*)
}

bool Tracker::IsDecayed(const Particle &P, double dx)	//Threshold check
{
	return GenMT->Rndm() > exp(-dx/(Const::C * P.LabSpace()));	//this should quit when particle decays too!
}

bool Tracker::IsDetectable(const Particle &P, bool print)	//Threshold check
{
	double Threshold = 0.0;

	switch (std::abs(P.Pdg()))
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
		case 311:	//neutral kaons
		case 2112:	//neutrons
			if (GenMT->Rndm() < 0.1)	//10% chance escape
				Threshold = 1 + P.EnergyKin();
			else
				Threshold = Get("Thres_Hadron");
			break;
		case 2212:	//protons
			Threshold = Get("Thres_Hadron");
			break;
		default:
			if (P.RealCharge() != 0)
				Threshold = Get("Thres_Hadron");
			else
				Threshold = 1 + P.EnergyKin();		//not detectable
			break;
	}

	if (print)
		std::cout << "detectable " << P.Pdg() << ": "
			  << P.EnergyKin() << " > " << Threshold << std::endl;
	return (P.EnergyKin() > Threshold);
}

void Tracker::Pi0Decay(Particle &pi0, Particle &pA, Particle &pB)
{
	//Vertex(Pi0);
	TVector3 bst(pi0.FourVector().BoostVector());
	TVector3 start(pi0.Position());		//starting point is vertex

	TLorentzVector gammaA(0, 0,  Const::MPion0/2.0, Const::MPion0/2.0); 
	TLorentzVector gammaB(0, 0, -Const::MPion0/2.0, Const::MPion0/2.0); 
	double the = GenMT->Uniform(-Const::pi, Const::pi);
	double phi = GenMT->Uniform(-Const::pi, Const::pi);
	gammaA.SetTheta(the);
	gammaB.SetTheta(the + Const::pi);
	gammaA.SetPhi(phi);
	gammaB.SetPhi(phi + Const::pi);
	gammaA.Boost(bst);
	gammaB.Boost(bst);

	pA = Particle(22, gammaA, start);	//here are the photons
	pB = Particle(22, gammaB, start);	//position should be different

	TVector3 moveA(gammaA.Vect().Unit());
	TVector3 moveB(gammaB.Vect().Unit());
	moveA *= GammaDecay(pA);
	moveB *= GammaDecay(pB);
	moveA += start;
	moveB += start;

	pA.SetPosition(moveA);
	pB.SetPosition(moveB);
}

//special use for neutrino beam
void Tracker::Focus(Particle &P)
{
	double radius = sqrt(pow(P.X(), 2) + pow(P.Y(), 2));
	double zlength = Zstart() + P.Z();

	P.SetTheta( atan2(radius, zlength) );
	P.SetPhi( atan2(P.Y(), P.X()) );
}
