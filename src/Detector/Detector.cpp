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
/*
double Detector::GausSmearing(TRandom3 *RanGen, double Mean, double Sigma)
{
	return RanGen->Gaus(Mean, Sigma);
}
*/


void Detector::SignalSmearing(TRandom3 *RanGen, Particle *P)
{
	double iM = P->M();
	double iEkin = P->Ekin();
	double iTheta = P->Theta();
	double iPhi = P->Phi();

	double SigmaA = GetElement("ResAngle");
	double SigmaZt = GetElement("Vertex");

	double Length = TrackLength(P);
	double Resolution;
	if (P->Pdg() == 11 || P->Pdg() == 22)
		Resolution = GetElement("ResGamma");
	else if (P->Pdg() == 13)
		Resolution = GetElement("ResMuon");
	else if (P->Pdg() == 211)
		Resolution = GetElement("ResPion");
	else if (P->Charge() != 0)
		Resolution = GetElement("ResHadron");
	else Resolution = 0;
	double SigmaEz = sqrt(iEkin) * Resolution;
	//double SigmaP = sqrt(1.5) * TheBox->GetElement("Vertex") * 8.0*iP*iP * 2.0 / (TheBox->GetElement("Bfield") * pow(TheBox->GetElement("Length"),2));
	double SigmaE = sqrt(pow(SigmaEz,2)*pow(Length,2) + pow(SigmaZt,2)*pow(iEkin, 2)	);

	if (IsDetectable(P))
	{
		P->SetE(iM + RanGen->Gaus(iEkin, SigmaE));

		P->SetTheta(RanGen->Gaus(iTheta, SigmaA));
		P->SetPhi(RanGen->Gaus(iPhi, SigmaA));

		P->SetX(RanGen->Gaus(P->X(), SigmaZt));
		P->SetY(RanGen->Gaus(P->Y(), SigmaZt));
		P->SetZ(RanGen->Gaus(P->Z(), SigmaZt));
	}
}

double Detector::TrackLength(Particle *P)	//This should not change *P
{
	if (P->Track() < 0.0)
	{
		double dx = 0.01;	//1 cm
		TVector3 Step(P->Position());
		TVector3 Proj(P->Position());
		Step.SetMag(dx);	//move in centimeter
		
		double Length = 0.0;
		double dE = P->E();
		double iM = P->M();

		double Beta;
		while (IsInside(Proj) && IsDetectable(P->Pdg(), P->Charge(), dE-iM))
		{
			Beta = sqrt(1-iM*iM/dE/dE);
			dE -= dx * EnergyLoss(Beta, iM);
	
			Proj += Step;
			Length += dx;
		}
	
		//*Where = Proj;
		//P->SetPosition(Proj);
	
		P->SetTrack(Length);
		return Length;
	}
	else
		return P->Track();
}

double Detector::InteractionLength()
{
	if (GetElement("Target") == 1)	//Argon
		return Kine::RadiationLength(18, 40);
	else return 0.0;
}

double Detector::EnergyLoss(double Beta, double Mass)
{
	if (GetElement("Target") == 1)
	{
		if (Mass > 10*Const::fMElectron) 
		{
			double Ret = Kine::Bethe(Beta, Mass, 1.3945, 188.0, 18, 40);	//LAr, code 1
			return Ret;
		}
		else return Mass/sqrt(1-Beta*Beta)/InteractionLength();
	}
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

bool Detector::IsDetectable(Particle *P)	//Threshold check
{
	bool Ret = false;

	if (P->Pdg() == 11)
	{
		if (P->Ekin() > GetElement("ThrGamma"))
			Ret = true;
	}
	else if (P->Pdg() == 22)	//photon
	{
		if (P->Ekin() > 2*Const::fMElectron)
			Ret = true;
	}
	else if (P->Pdg() == 13)	//muon
	{
		if (P->Ekin() > GetElement("ThrMuon"))
			Ret = true;
	}
	else if (P->Pdg() == 211)	//pion
	{
		if (P->Ekin() > GetElement("ThrPion"))
			Ret = true;
	}
	else if (P->Charge() != 0)	//proton
	{
		if (P->Ekin() > GetElement("ThrHadron"))
			Ret = true;
	}

	return Ret;
}

bool Detector::IsDetectable(int Pdg, int Charge, double Ekin)	//Threshold check
{
	bool Ret = false;

	if (Pdg == 11)	//photon
	{
		if (Ekin > GetElement("ThrGamma"))
			Ret = true;
	}
	else if (Pdg == 22)	//photon
	{
		if (Ekin > 2*Const::fMElectron)
			Ret = true;
	}
	else if (Pdg == 13)	//muon
	{
		if (Ekin > GetElement("ThrMuon"))
			Ret = true;
	}
	else if (Pdg == 211)	//pion
	{
		if (Ekin > GetElement("ThrPion"))
			Ret = true;
	}
	else if (Charge != 0)	//proton
	{
		if (Ekin > GetElement("ThrHadron"))
			Ret = true;
	}

	return Ret;
}

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
	bool Ret;
	if (P.X() < GetXsize() &&
	    P.Y() < GetYsize() &&
	    P.Z() < GetZsize())
		return true;
	else return false;
}

double Detector::GetXsize()
{
	double F = GetElement("Fiducial");
	return F > 0.0 ? F *GetElement("Width") : GetElement("Width");
}

double Detector::GetXstart()
{
	double F = GetElement("Fiducial");
	return F > 0.0 ? 0.5*GetElement("Width")*(1-F) : 0.0;
}

double Detector::GetXend()
{
	return GetXsize() + GetXstart();
}

double Detector::GetYsize()
{
	double F = GetElement("Fiducial");
	return F > 0.0 ? F * GetElement("Height") : GetElement("Height");
}

double Detector::GetYstart()
{
	double F = GetElement("Fiducial");
	return F > 0.0 ? 0.5*GetElement("Height")*(1-F) : 0.0;
}

double Detector::GetYend()
{
	return GetYsize() + GetYstart();
}

double Detector::GetZsize()
{
	double F = GetElement("Fiducial");
	return F > 0.0 ? F * GetElement("Length") : GetElement("Length");
}

double Detector::GetZstart()
{
	double F = GetElement("Fiducial");
	return F > 0.0 ? 0.5*GetElement("Length")*(1-F) : 0.0;
}

double Detector::GetZend()
{
	return GetZsize() + GetZstart();
}
