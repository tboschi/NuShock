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
			if ( IsInsideLAr(P) && P.TrackIn() < ZsizeLAr() && Ratio > Get("Containment") )
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
	double step = 0.05, dx;	//step of 5cm
	double LengthIn = 0.0, LengthBack = 0.0, LengthOut = 0.0;
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

			P.SetTrackIn(Depth);
			P.SetTrackOut(0);
			break;
		case 22:
			P.SetTrackIn(GammaDecay(P));
			P.SetTrackOut(0);
			break;
		case 13:
			while (IsDetectable(clone) && !IsDecayed(clone, LengthIn+LengthOut))
				//quit when particle decays too!
			{
				bool InOut;
				double dx = GenMT->Gaus(step, Get("Vertex"));
				double loss = EnergyLoss(clone, InOut);
				double dE = dx * loss;

				if (InOut)	//inside detector
					LengthIn += dx;
				else		//outside detector
					LengthOut += dx;

				TVector3 pos = clone.Position() + dx * dir;	//move by dx
				clone.SetEnergy(clone.Energy() - dE);
				clone.SetPosition(pos);
			}

			P.SetTrackIn(LengthIn);
			P.SetTrackOut(LengthOut);
			break;
		case 211:
			while (IsDetectable(clone) && !IsDecayed(clone, LengthIn+LengthOut))
			{
				if (cover >= length)	//distance covered more than interaction length
				{
					int mult = ceil(pow(clone.EnergyKin() * 1e6, 0.18) * 0.15);
					if (mult > 1)	//multiplicity of hadronic interaction
					{
						++layer;	//number of interactions
						clone.SetEnergyKin(clone.EnergyKin()/mult);
						//divide energy by number (mult) of duaghter particles
					}
					length = -1;	//repeat
				}

				if (length < 0)
				{
					length = GenMT->Exp(RadiationLength(clone));
					cover = 0;
				}

				bool InOut;
				double dx = GenMT->Gaus(step, Get("Vertex"));
				double dE = dx * EnergyLoss(P, InOut);	//inside material

				cover += dx;	//total distance covered in one layer

				if (InOut)	//inside detector
					LengthIn += dx;
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

			P.SetTrackIn(LengthIn);
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

double Tracker::CriticalEnergy(Detector::Material Element)	//assuming same for positron and electron
{
	switch (Element)
	{
		case LAr:
			return 0.03284;	//GeV
		case GasAr:
			return 0.03803;	//GeV
		case Fe:
			return 0.02168;	//GeV
		case Pb:
			return 0.00743;	//GeV
		default:
			return 0;
	}
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

double Tracker::RadiationLength(Detector::Material Element, bool Nuclear)
{
	if (!Nuclear)
		switch (Element)
		{
			case LAr:
				return 19.55/1.3945 / 100;
			case GasAr:
				return 19.55/0.1020 / 100;
			case Fe:
				return 13.84/7.874 / 100;
			case Pb:
				return 6.37 /11.34 / 100;
			default:
				return 0;
		}
	else
		switch (Element)
		{
			case LAr:
				return 119.7/1.3945 / 100;
			case GasAr:
				return 119.7/0.1020 / 100;
			case Fe:
				return 132.1/7.874 / 100;
			case Pb:
				return 199.6/11.34 / 100;
			default:
				return 0;
		}
}

double Tracker::EnergyLoss(const Particle &P, bool &inout)
{
	if (IsInsideLAr(P))
	{
		inout = true;
		return BetheLoss(P, GetMaterial("TargetLAr"));
	}
	else if (IsInsideFGT(P))
	{
		inout = true;
		return BetheLoss(P, GetMaterial("TargetFGT"));
	}
	else
	{
		inout = false;
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
			if (P.Charge() != 0)
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
