#include "Detector.h"

Detector::Detector(std::string ConfigName, TRandom3 *Random)
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
	RanGen = Random;
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

double Detector::GetElement(std::string Key, bool S)
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
	double iM = P->M();
	double iEkin = P->Ekin();
	double iP = P->P();
	double iTheta = P->Theta();
	double iPhi = P->Phi();
	double SigmaE, SigmaP, SigmaA; 

	if (P->Pdg() == 11 || P->Pdg() == 22)
	{
		SigmaA = GetElement("Angle_Gamma") / Const::fDeg;
		double StatE = GetElement("Energ_Gamma") / sqrt(iEkin);
		double SystE = GetElement("Ebias_Gamma");
		SigmaE = sqrt(pow(StatE, 2)+pow(SystE, 2)) * iEkin;

		P->SetEnergyKin(RanGen->Gaus(iEkin, SigmaE));
		P->SetTheta(RanGen->Gaus(iTheta, SigmaA));
		P->SetPhi(RanGen->Gaus(iPhi, SigmaA));
	}
	else if (P->Pdg() == 13)
	{
		SigmaA = GetElement("Angle_Muon") / Const::fDeg;
		double Ratio = P->TrackIn()/(P->TrackIn()+P->TrackOut());
		//std::cout << "13  " << Ratio << "\t" << P->TrackIn() << "\t" << P->TrackOut() << std::endl;
		if (Ratio > GetElement("Containment"))	//90% of track inside
		{
			double StatP = GetElement("Range_Muon") / iP;
			SigmaP = StatP * iP;
		}
		else
		{
			double StatP = GetElement("Exiti_Muon");
			SigmaP = StatP * iP;
		}


		P->SetMomentum(RanGen->Gaus(iP, SigmaP));
		P->SetTheta(RanGen->Gaus(iTheta, SigmaA));
		P->SetPhi(RanGen->Gaus(iPhi, SigmaA));
	}
	else if (P->Pdg() == 211)
	{
		SigmaA = GetElement("Angle_Pion") / Const::fDeg;
		double Ratio = P->TrackIn()/(P->TrackIn()+P->TrackOut());
		//std::cout << "211 " << Ratio << "\t" << P->TrackIn() << "\t" << P->TrackOut() << std::endl;
		if (!P->IsShower() && Ratio > GetElement("Containment"))	//pion track, not shower
		{
			double StatP = GetElement("Range_Pion") / iP;
			SigmaP = StatP * iP;
		}
		else
		{
			double StatP = GetElement("Exiti_Pion");
			SigmaP = StatP * iP;
		}

		P->SetMomentum(RanGen->Gaus(iP, SigmaP));
		P->SetTheta(RanGen->Gaus(iTheta, SigmaA));
		P->SetPhi(RanGen->Gaus(iPhi, SigmaA));
	}
	else if (P->Charge() != 0)	//other and protons?
	{
		SigmaA = GetElement("Angle_Hadron") / Const::fDeg;
		double StatE = GetElement("Energ_Hadron") / sqrt(iEkin);
		double SystE = GetElement("Ebias_Hadron");
		SigmaE = sqrt(pow(StatE, 2)+pow(SystE, 2)) * iEkin;

		P->SetEnergyKin(RanGen->Gaus(iEkin, SigmaE));
		P->SetTheta(RanGen->Gaus(iTheta, SigmaA));
		P->SetPhi(RanGen->Gaus(iPhi, SigmaA));
	}

	P->SetTrackIn(RanGen->Gaus(P->TrackIn(), GetElement("Vertex")));
	if (P->TrackIn() < 0)
		P->SetTrackIn(0);
	P->SetTrackOut(RanGen->Gaus(P->TrackOut(), GetElement("Vertex")));
	if (P->TrackOut() < 0)
		P->SetTrackOut(0);
}

void Detector::TrackLength(Particle *P)	//This should not change *P
{
	if (P->Pdg() == 11)
	{
		double tmax = log(P->E()/CriticalEnergy()) - 1.0;
		double Depth = RadiationLength() * RanGen->Gaus(tmax, 0.3*tmax);	//justify this
		P->SetTrackIn(Depth);
	}
	else if (P->Pdg() == 22)
	{
		P->SetTrackIn(GammaDecay());
	}
	else if (P->Pdg() == 13)
	{
		double iE = P->E();
		double dE;
		double dStep = 0.01;	//cm step
		double dx, LengthIn = 0, LengthOut = 0;
		
		TVector3 Step(P->Direction().Unit());
		TVector3 Pos(P->Position());
		TVector3 Start(P->Position());


		bool InOut;
		while (IsDetectable(P) && !IsDecayed(P, dStep))		//this should quit when particle decays too!
		{
			dx = RanGen->Gaus(dStep, GetElement("Vertex"));

			dE = dx * EnergyLoss(P, InOut);	//inside material

			if (InOut)
				LengthIn += dx;			//even b2b det?
			else
				LengthOut += dx;

			Pos += Step*dx;		//move by dx
			P->SetEnergy(P->E() - dE);
			P->SetPosition(Pos);
		}

		P->SetTrackIn(LengthIn);
		P->SetTrackOut(LengthOut);
		P->SetEnergy(iE);		//reset to original energy
		P->SetPosition(Start);		//reset to original position
	}
	else if (P->Pdg() == 211)
	{
		double iE = P->E();
		double iM = P->M();
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
			dx = RanGen->Gaus(dStep, GetElement("Vertex"));
			double Length = RanGen->Exp(RadiationLength(1));
			double Cover = 0;

			while (IsDetectable(P) && !IsDecayed(P, dStep) && Cover < Length)	//this should quit when particle decays too!
			{
				dx = RanGen->Gaus(dStep, GetElement("Vertex"));

				dE = dx * EnergyLoss(P, InOut);	//inside material
				Cover += dx;

				if (InOut)
					LengthIn += dx;			//even b2b det?
				else
					LengthOut += dx;

				Pos += Step*dx;		//move by dx
				P->SetEnergy(P->E() - dE);
				P->SetPosition(Pos);
			}

			int Mult = ceil(pow(P->Ekin() * 1e6, 0.18) * 0.15);
			if (Mult > 1)
			{
				++Layer;
				P->SetEnergyKin(P->Ekin()/Mult);
			}
		}

		if (RanGen->Rndm() > 1.0/Layer)		//need something better
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
	return RanGen->Exp(Path);
}

double Detector::CriticalEnergy()	//assuming same for positron and electron
{
	int Target = GetElement("InTarget");
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
	int Target = GetElement("InTarget");
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
		return BetheLoss(P, GetElement("InTarget"));
	}
	else if (IsContained(P) && GetElement("BackTarget") != 0)
	{
		Contained = true;
		return BetheLoss(P, GetElement("BackTarget"));
	}
	else
	{
		Contained = false;
		return BetheLoss(P, GetElement("OutTarget"));
	}
}

double Detector::BetheLoss(Particle *P, double Target)
{
	if (Target == 1)
		return Kine::Bethe(P->GetP4().Beta(), P->M(), 1.3945, 188.0, 18, 40);	//LAr, code 1
	else if (Target== 2)
		return Kine::Bethe(P->GetP4().Beta(), P->M(), 0.1020, 188.0, 18, 40);	//GasAr, code 2
	else if (Target == 3)
		return Kine::Bethe(P->GetP4().Beta(), P->M(), 7.874, 286.0, 26, 56);	//Fe, code 3
	else return 0.0;
}

bool Detector::IsDecayed(Particle *P, double dx)	//Threshold check
{
	return RanGen->Rndm() > exp(-dx/(Const::fC * P->LabSpace()));	//this should quit when particle decays too!
}

bool Detector::IsDetectable(Particle *P)	//Threshold check
{
	bool Ret = false;

	if (P->Pdg() == 11 || P->Pdg() == 22)	//electron or photon
	{
		if (P->Ekin() > GetElement("Thres_Gamma"))
			Ret = true;
	}
	else if (P->Pdg() == 13)	//muon
	{
		if (P->Ekin() > GetElement("Thres_Muon"))
			Ret = true;
	}
	else if (P->Pdg() == 211)	//pion
	{
		if (P->Ekin() > GetElement("Thres_Pion"))
			Ret = true;
	}
	else if (P->Charge() != 0)	//hadron
	{
		if (P->Ekin() > GetElement("Thres_Hadron"))
			Ret = true;
	}

	return Ret;
}

bool Detector::IsContained(Particle *P)
{
	if (P->X() > GetXstart() && P->X() < GetXend() &&
	    P->Y() > GetYstart() && P->Y() < GetYend() &&
	    P->Z() > GetZstart() )
		return true;
	else return false;
}

bool Detector::IsInside(Particle *P)
{
	if (P->X() > GetXstart() && P->X() < GetXend() &&
	    P->Y() > GetYstart() && P->Y() < GetYend() &&
	    P->Z() > GetZstart() && P->Z() < GetZend() )
		return true;
	else return false;
}

double Detector::GetXsize()
{
	double F = GetElement("Fiducial");
	return F > 0.0 ? F *GetElement("Width", 1) : GetElement("Width", 1);
}

double Detector::GetXstart()
{
	double F = GetElement("Fiducial");
	return F > 0.0 ? 0.5*GetElement("Width", 1)*(1-F) : 0.0;
}

double Detector::GetXend()
{
	return GetXsize() + GetXstart();
}

double Detector::GetYsize()
{
	double F = GetElement("Fiducial");
	return F > 0.0 ? F * GetElement("Height", 1) : GetElement("Height", 1);
}

double Detector::GetYstart()
{
	double F = GetElement("Fiducial");
	return F > 0.0 ? 0.5*GetElement("Height", 1)*(1-F) : 0.0;
}

double Detector::GetYend()
{
	return GetYsize() + GetYstart();
}

double Detector::GetZsize()
{
	double F = GetElement("Fiducial");
	return F > 0.0 ? F * GetElement("Length", 1) : GetElement("Length", 1);
}

double Detector::GetZstart()
{
	double F = GetElement("Fiducial");
	return F > 0.0 ? 0.5*GetElement("Length", 1)*(1-F) : 0.0;
}

double Detector::GetZend()
{
	return GetZsize() + GetZstart();
}
