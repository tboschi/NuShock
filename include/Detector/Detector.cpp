#include "Detector.h"

Detector::Detector(std::string ConfigName) :
	GenMT(new TRandom3(0))
{
	std::ifstream ConfigFile(ConfigName.c_str());

	std::string Line, Key, Name;
	std::stringstream ssL;
	double Element;
	Eff = false;
	EnergyEfficiency EnEff;

	while (std::getline(ConfigFile, Line))
	{
		if (Line[0] == '#') continue;

		ssL.str("");
		ssL.clear();
		ssL << Line;

		if (Line[0] == '&')
		{
			ssL >> Key >> Name;
			Key.erase(Key.begin());
			mapEfficiencyE[Key] = Name;
		}
		else if (Line[0] == '$')
		{
			ssL >> Key >> Name;
			Key.erase(Key.begin());
			mapEfficiencyM[Key] = Name;
		}
		else
		{
			ssL >> Key >> Element;
			mapDetector[Key] = Element;
		}
	}

	FuncFile = new TFile();
}

std::vector<std::string> Detector::ListKey()
{
	std::vector<std::string> List;
	std::map<std::string, double>::iterator it;
	for (it = mapDetector.begin(); it != mapDetector.end(); ++it)
		List.push_back(it->first);
	return List;
}

std::vector<std::string> Detector::ListChannel()
{
	std::vector<std::string> List;
	std::map<std::string, std::string>::iterator it;
	for (it = mapEfficiencyE.begin(); it != mapEfficiencyE.end(); ++it)
		List.push_back(it->first);
	for (it = mapEfficiencyM.begin(); it != mapEfficiencyM.end(); ++it)
		List.push_back(it->first);
	return List;
}

double Detector::Get(std::string Key, bool S)
{
	if (S && mapDetector.count("Scale") > 0)	//there is a scale element
		return mapDetector["Scale"] * mapDetector[Key];
	else
		return mapDetector[Key];
}

double Detector::Efficiency(double Energy, double Mass)
{
	//std::cout << "eff0" << std::endl;
	int Ebin = hhFunc->GetXaxis()->FindBin(Energy);
	//std::cout << "eff1" << std::endl;
	int Mbin = hhFunc->GetYaxis()->FindBin(Mass);
	//std::cout << "eff2" << std::endl;

	return hhFunc->GetBinContent(Ebin, Mbin);
}

void Detector::SetEfficiency(std::string Channel, char Coupling)
{
	if (FuncFile != 0 && FuncFile->IsOpen())
		FuncFile->Close();

	if (Coupling == 'E')
		FuncFile = new TFile(mapEfficiencyE[Channel].c_str(), "OPEN");
	else if (Coupling == 'M')
		FuncFile = new TFile(mapEfficiencyM[Channel].c_str(), "OPEN");

 	hhFunc = dynamic_cast<TH2D*> (FuncFile->Get("hhfunc"));
}

/*
double Detector::Efficiency(std::string Channel, double Energy)
{
	double diff = 10.0, sdiff = 10.0;
	double EA = 0, EB = 0;
	double fA, fB;

	std::map<std::string, std::vector<EnergyEfficiency> >::iterator imap;
	std::vector<EnergyEfficiency>::iterator it, iA, iB;

	if ((imap = mapEfficiency.find(Channel)) != mapEfficiency.end())
	{
		for (it = mapEfficiency[Channel].begin(); it != mapEfficiency[Channel].end(); ++it)
		{
			if (Energy == it->E)
				return it->f;
			else if (Energy > it->E)
			{
				if (diff > (Energy - it->E))
				{
					diff = Energy - it->E;
					EA = it->E;
					fA = it->f;
				}
			}
			else 
			{
				if (sdiff > (it->E - Energy))
				{
					sdiff = it->E - Energy;
					EB = it->E;
					fB = it->f;
				}
			}
		}
		if (EA == 0)
			return fB;
		else if (EB == 0)
			return fA;
		else
		{
			double mFactor = (fB - fA)/(EB - EA);
			double qFactor = fA - mFactor * EA;
		
			return Energy*mFactor+qFactor;
		}
	}
	else return 1.0;
}
*/

//detector effect and smearing

void Detector::SignalSmearing(Particle *P)
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

void Detector::TrackLength(Particle *P)	//This should not change *P
{
	if (P->Pdg() == 11)
	{
		double tmax = log(P->Energy()/CriticalEnergy()) - 1.0;
		double Depth = RadiationLength() * GenMT->Gaus(tmax, 0.3*tmax);	//justify this
		P->SetTrackIn(Depth);
	}
	else if (P->Pdg() == 22)
	{
		P->SetTrackIn(GammaDecay());
	}
	else if (P->Pdg() == 13)
	{
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
	}
	else if (P->Pdg() == 211)
	{
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
	}
}

double Detector::GammaDecay()
{
	double Path = 9.0/7.0 * RadiationLength(0);
	return GenMT->Exp(Path);
}

double Detector::CriticalEnergy()	//assuming same for positron and electron
{
	int Target = Get("InTarget");
	if (Target == 1)
		return 0.03284;	//GeV
	else if (Target == 2)
		return 0.03803;	//GeV
	else if (Target == 2)
		return 0.02168;	//GeV
	else return 0;
}

double Detector::RadiationLength(bool Nuclear)
{
	int Target = Get("InTarget");
	double Ret;
	if (!Nuclear) switch (Target)
	{
		case 1:
			Ret = 19.55/1.3945;
			break;
		case 2:
			Ret = 19.55/0.1020;
			break;
		case 3:
			Ret = 13.84/7.874;
			break;
	}
	else switch (Target)
	{
		case 1:
			Ret = 119.7/1.3945;
			break;
		case 2:
			Ret = 119.7/0.1020;
			break;
		case 3:
			Ret = 132.1/7.874;
			break;
	}

	return Ret / 100.0;
}

double Detector::EnergyLoss(Particle *P, bool &Contained)
{
	if (IsContained(P) && IsInside(P))
	{
		Contained = true;
		return BetheLoss(P, Get("InTarget"));
	}
	else if (IsContained(P) && Get("BackTarget") != 0)
	{
		Contained = true;
		return BetheLoss(P, Get("BackTarget"));
	}
	else
	{
		Contained = false;
		return BetheLoss(P, Get("OutTarget"));
	}
}

double Detector::BetheLoss(Particle *P, double Target)
{
	if (Target == 1)
		return Bethe(P->Beta(), P->Mass(), 1.3945, 188.0, 18, 40);	//LAr, code 1
	else if (Target== 2)
		return Bethe(P->Beta(), P->Mass(), 0.1020, 188.0, 18, 40);	//GasAr, code 2
	else if (Target == 3)
		return Bethe(P->Beta(), P->Mass(), 7.874, 286.0, 26, 56);	//Fe, code 3
	else return 0.0;
}

double Detector::Bethe(double Beta, double Mass, double Density, double I, int Z, int A)
{
	double K = 0.307075;	//From PDG MeV mol-1 cm2
	double Beta2 = Beta*Beta;
	double Gamma = 1.0/sqrt(1-Beta2);
	double Gamma2 = Gamma*Gamma;
	double e2M = Const::fMElectron / Mass;
	double Wmax = (2000 * Const::fMElectron * Beta2 * Gamma2) / (1 + 2*Gamma*e2M + e2M*e2M);
	double LogArg = 2000 * Const::fMElectron * Beta2 * Gamma2 * Wmax / (1e-12*I*I); 	//Everything in MeV
	return 0.1 * Density * (K * Z) / (A * Beta2) * (0.5 * log (LogArg) - Beta2);	//stopping power in GeV/m (0.1*)
}

bool Detector::IsDecayed(Particle *P, double dx)	//Threshold check
{
	return GenMT->Rndm() > exp(-dx/(Const::fC * P->LabSpace()));	//this should quit when particle decays too!
}

bool Detector::IsDetectable(Particle *P)	//Threshold check
{
	bool Ret = false;

	if (P->Pdg() == 11 || P->Pdg() == 22)	//electron or photon
	{
		if (P->EnergyKin() > Get("Thres_Gamma"))
			Ret = true;
	}
	else if (P->Pdg() == 13)	//muon
	{
		if (P->EnergyKin() > Get("Thres_Muon"))
			Ret = true;
	}
	else if (P->Pdg() == 211)	//pion
	{
		if (P->EnergyKin() > Get("Thres_Pion"))
			Ret = true;
	}
	else if (P->Charge() != 0)	//hadron
	{
		if (P->EnergyKin() > Get("Thres_Hadron"))
			Ret = true;
	}

	return Ret;
}

bool Detector::IsContained(Particle *P)
{
	if (P->X() > Xstart() && P->X() < Xend() &&
	    P->Y() > Ystart() && P->Y() < Yend() &&
	    P->Z() > Zstart() )
		return true;
	else return false;
}

bool Detector::IsInside(Particle *P)
{
	if (P->X() > Xstart() && P->X() < Xend() &&
	    P->Y() > Ystart() && P->Y() < Yend() &&
	    P->Z() > Zstart() && P->Z() < Zend() )
		return true;
	else return false;
}

double Detector::Xsize()
{
	double F = Get("Fiducial");
	return F > 0.0 ? F *Get("Width", 1) : Get("Width", 1);
}

double Detector::Xstart()
{
	double F = Get("Fiducial");
	return F > 0.0 ? 0.5*Get("Width", 1)*(1-F) : 0.0;
}

double Detector::Xend()
{
	return Xsize() + Xstart();
}

double Detector::Ysize()
{
	double F = Get("Fiducial");
	return F > 0.0 ? F * Get("Height", 1) : Get("Height", 1);
}

double Detector::Ystart()
{
	double F = Get("Fiducial");
	return F > 0.0 ? 0.5*Get("Height", 1)*(1-F) : 0.0;
}

double Detector::Yend()
{
	return Ysize() + Ystart();
}

double Detector::Zsize()
{
	double F = Get("Fiducial");
	return F > 0.0 ? F * Get("Length", 1) : Get("Length", 1);
}

double Detector::Zstart()
{
	double F = Get("Fiducial");
	return F > 0.0 ? 0.5*Get("Length", 1)*(1-F) : 0.0;
}

double Detector::Zend()
{
	return Zsize() + Zstart();
}


double Detector::DecayProb(Neutrino *N)	//reaching the detector and decaying
{
	if (TheN->EnergyKin() < 0.0)
		return 0.0;
	else
	{
		//double Total = GetUserData() < 0.0 ? TheGamma->Total() : GetUserData();
		double Total = TheN->DecayTotal();
		double Ratio = TheN->DecayBranch(); 
		//std::cout << "Ratio " << Ratio << std::endl;
		double Length = Const::fM2GeV * Get("Baseline");
		double Lambda = Const::fM2GeV * Zsize();
		double Lorentz = TheN->Beta() * TheN->Gamma();
			//sqrt(pow(TheN->Energy()/TheN->Mass(), 2) - 1.0);	//betagamma, to invert

		return exp(- Total * Length / Lorentz) * (1 - exp(- Total * Lambda / Lorentz)) * Ratio;
	}
}

void Detector::Focus(Neutrino *Nu)
{
	double Radius = sqrt(pow(Get("Width"), 2) + pow(Get("Height"), 2));
	double SigmaT = atan2(Radius, Get("Baseline"));					//3sigma will be inside the detector (better distribution needed)

	Nu->SetTheta( abs(GenMT->Gaus(0, SigmaT)) );	
	Nu->SetPhi( GenMT->Uniform(-Const::fPi, Const::fPi) );
}
