#include "Detector.h"

Detector::Detector(std::string ConfigName)
{
	std::ifstream ConfigFile(ConfigName.c_str());

	std::string Line, Key;
	std::stringstream ssL;
	double Element;
	EnergyEfficiency EnEff;

	while (std::getline(ConfigFile, Line))
	{
		if (Line[0] == '#') continue;

		ssL.str("");
		ssL.clear();
		ssL << Line;

		if (Line[0] == '_')
		{
			ssL >> Key >> EnEff.E >> EnEff.f;
			Key.erase(Key.begin());
			mapEfficiency[Key].push_back(EnEff);
		}
		else
		{
			ssL >> Key >> Element;
			mapDetector[Key] = Element;
		}
	 }
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
	std::map<std::string, std::vector<EnergyEfficiency> >::iterator it;
	for (it = mapEfficiency.begin(); it != mapEfficiency.end(); ++it)
		List.push_back(it->first);
	return List;
}

double Detector::GetElement(std::string Key)
{
	return mapDetector[Key];
}

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

//detector effect and smearing

double Detector::GausSmearing(TRandom3 *RanGen, double Mean, double Sigma)
{
	return RanGen->Gaus(Mean, Sigma);
}

bool Detector::IsDetectable(Particle *P)	//Threshold check
{
	bool Ret = false;
	if (P->Pdg() == 14)	//muon
	{
		if (P->Ekin() > GetElement("ThrMuon"))
			Ret = true;
	}
	else if (P->Pdg() == 211)	//pion
	{
		if (P->Ekin() > GetElement("ThrPion"))
			Ret = true;
	}
	else if (P->Pdg() == 2212)	//proton
	{
		if (P->Ekin() > GetElement("ThrNuclear"))
			Ret = true;
	}
	else if (P->Pdg() == 22)	//photon
	{
		if (P->Ekin() > GetElement("ThrGamma"))
			Ret = true;
	}

	return Ret;
}

void Detector::SignalSmearing(TRandom3 *RanGen, Particle *P)
{
	double iE = P->E();
	double iTheta = P->Theta();
	double iPhi = P->Phi();

	double SigmaEz;
	double SigmaA = GetElement("ResAngle");
	double SigmaZt = GetElement("Vertex");			//2 is for fiducial volume of 35%

	double Length = TrackLength(P);
	switch(P->Pdg())
	{
		case 14:
			SigmaEz = sqrt(iE) * GetElement("ResMuon");
			break;
		case 211:
			SigmaEz = sqrt(iE) * GetElement("ResPion");
			break;
		case 2212:
			SigmaEz = sqrt(iE) * GetElement("ResNuclear");
			break;
		default:
			break;
	}
	//double SigmaP = sqrt(1.5) * TheBox->GetElement("Vertex") * 8.0*iP*iP * 2.0 / (TheBox->GetElement("Bfield") * pow(TheBox->GetElement("Length"),2));

	double SigmaE = sqrt(pow(SigmaEz,2)*pow(GetElement("Length"),2) + pow(SigmaZt,2)*iE);

	P->SetE(RanGen->Gaus(iE, SigmaE));
	P->SetTheta(RanGen->Gaus(iTheta, SigmaA));
	P->SetPhi(RanGen->Gaus(iPhi, SigmaA));
	P->SetX(RanGen->Gaus(P->X(), SigmaZt));
	P->SetY(RanGen->Gaus(P->Y(), SigmaZt));
	P->SetZ(RanGen->Gaus(P->Z(), SigmaZt));
}

double Detector::TrackLength(Particle *P)	//This should not change *P
{
	double dx = 0.01;	//1 cm
	TVector3 Step(P->Position());
	TVector3 Proj(P->Position());
	Step.SetMag(dx);	//move in centimeter
	
	double Length = 0.0;
	double dE = P->E();
	double iM = P->M();

	double Beta;
	while (IsInside(Proj) && IsDetectable(P))
	{
		Beta = sqrt(1-iM*iM/dE/dE);
		dE -= dx * EnergyLoss(Beta, iM);

		Proj += Step;
		Length += dx;
	}

	//*Where = Proj;
	//P->SetPosition(Proj);

	return Length;
}

double Detector::EnergyLoss(double Beta, double Mass)
{
	if (GetElement("Target") == 1)
		return Kine::Bethe(Beta, Mass, 1.3945, 188.0, 18, 40);	//LAr, code 1
	else return 0.0;
}

/*
std::string Detector::ParticleID(TLorentzVector *Track)
{
	std::string Ret = "undefined";

	std::map<std::string, double>::iterator it;
	for (it = mapMass.begin(); it != mapMass.end(); ++it)
	{
		if (fabs(Track->M() - it->second) < 0.005)
			Ret.assign(it->first);
	}

	return Ret;
}
*/

bool Detector::IsInside(Particle *P)
{
	if (P->X() < GetXsize() &&
	    P->Y() < GetYsize() &&
	    P->Z() < GetZsize())
		return true;
	else return false;
}

bool Detector::IsInside(TVector3 &P)
{
	if (P.X() < GetXsize() &&
	    P.Y() < GetYsize() &&
	    P.Z() < GetZsize())
		return true;
	else return false;
}

double Detector::GetXsize()
{
	return GetElement("Fiducial") > 0.0 ? GetElement("Width")*GetElement("Fiducial") : GetElement("Width");
}

double Detector::GetYsize()
{
	return GetElement("Fiducial") > 0.0 ? GetElement("Heigth")*GetElement("Fiducial") : GetElement("Heigth");
}

double Detector::GetZsize()
{
	return GetElement("Fiducial") > 0.0 ? GetElement("Length")*GetElement("Fiducial") : GetElement("Length");
}
