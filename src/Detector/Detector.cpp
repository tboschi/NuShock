#include "Detector.h"

Detector::Detector(std::string ConfigName)
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

		if (Line == "EffFile")
		{
			ssL >> Key >> Name;
			EffFile = new TFile(Name.c_str(), "OPEN");
			Eff = true;

			hTemp = (TH1D*) EffFile->Get("efficiency");
			hEfficiency = (TH1D*) hTemp->Clone();
			hEfficiency->SetDirectory(0);
			EffFile->Close();
		}
		else if (Line[0] == '_')
		{
			ssL >> Key >> Key >> Name;
			Key.erase(Key.begin());
			mapEfficiency[Key] = Name;
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
	std::map<std::string, std::string>::iterator it;
	for (it = mapEfficiency.begin(); it != mapEfficiency.end(); ++it)
		List.push_back(it->first);
	return List;
}

double Detector::GetElement(std::string Key)
{
	return mapDetector[Key];
}

double Detector::Efficiency(std::string Channel, double Energy, double Mass)
{
	TFile FuncFile(mapEfficiency[Channel].c_str(), "OPEN");
 	TH2D *hhEff = dynamic_cast<TH2D*> (FuncFile.Get("hhfunc"));

	int Ebin = hhEff->GetXaxis()->FindBin(Energy);
	int Mbin = hhEff->GetYaxis()->FindBin(Mass);

	return hhEff->GetBinContent(Ebin, Mbin);
}

double Detector::Background(std::string Channel, double Energy)
{
	TFile FuncFile(mapEfficiency[Channel].c_str(), "OPEN");
 	TH1D *hBack = dynamic_cast<TH1D*> (FuncFile.Get("hcut"));

	int Ebin = hBack->GetXaxis()->FindBin(Energy);
	double POT = GetElement("POT/s") * pow(10, GetElement("Years"));
	double Norm = GetElement("Events") * GetElement("Weigth") * POT / 1.0e6;

	return hBack->GetBinContent(Ebin) * Norm;
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
/*
double Detector::GausSmearing(TRandom3 *RanGen, double Mean, double Sigma)
{
	return RanGen->Gaus(Mean, Sigma);
}
*/

double Detector::EnergySigma(Particle *P)
{
	double SigmaE;

	if (P->Pdg() == 11)
	{
		double SysErr = P->E()*GetElement("SysElectron");
		double StatErr = sqrt(P->E())*GetElement("ResElectron");
		SigmaE = sqrt(pow(SysErr, 2)+pow(StatErr, 2));
	}
	else if (P->Pdg() == 22)
	{
		double SysErr = P->E()*GetElement("SysGamma");
		double StatErr = sqrt(P->E())*GetElement("ResGamma");
		SigmaE = sqrt(pow(SysErr, 2)+pow(StatErr, 2));
	}
	else if (P->Pdg() == 13)
	{
		double SysErr = P->E()*GetElement("SysMuon");
		double StatErr = sqrt(P->E())*GetElement("ResMuon");
		double Ratio = P->TrackIn()/(P->TrackIn()+P->TrackOut());
		SigmaE = sqrt(pow(SysErr, 2)+pow(StatErr, 2));
		if (Ratio < 0.90)
			SigmaE *= 4;
	}
	else if (P->Pdg() == 211)
	{
		double SysErr = P->E()*GetElement("SysPion");
		double StatErr = sqrt(P->E())*GetElement("ResPion");
		double Ratio = P->TrackIn()/(P->TrackIn()+P->TrackOut());
		SigmaE = sqrt(pow(SysErr, 2)+pow(StatErr, 2));
		if (Ratio < 0.90)
			SigmaE *= 4;
	}
	else if (P->Charge() != 0)
	{
		double SysErr = P->E()*GetElement("SysHadron");
		double StatErr = sqrt(P->E())*GetElement("ResHadron");
		SigmaE = sqrt(pow(SysErr, 2)+pow(StatErr, 2));
	}
	else SigmaE = 0;

	return SigmaE;
}

void Detector::SignalSmearing(TRandom3 *RanGen, Particle *P)
{
	double iM = P->M();
	double iE = P->E();
	double iTheta = P->Theta();
	double iPhi = P->Phi();

	double SigmaA = GetElement("ResAngle");
	double SigmaZt = GetElement("Vertex");

	double Length = TrackLength(RanGen, P);

	double SigmaE = EnergySigma(P);
	//double SigmaEz = sqrt(iEkin) * Resolution;
	//double SigmaP = sqrt(1.5) * TheBox->GetElement("Vertex") * 8.0*iP*iP * 2.0 / (TheBox->GetElement("Bfield") * pow(TheBox->GetElement("Length"),2));
	//double SigmaE = sqrt(pow(SigmaEz,2)*pow(Length,2) + pow(SigmaZt,2)*pow(iEkin, 2));

	//double SigmaP = GetElement("SysMomentum")*P->P()/sqrt(Length)
	//	        GetElement("ResMomentum")*pow(P->P(), 2)/pow(Length, 2.5);

	//P->SetEnergy(iM + RanGen->Gaus(iEkin, SigmaE));
	P->SetEnergy(RanGen->Gaus(iE, SigmaE));

	P->SetTheta(RanGen->Gaus(iTheta, SigmaA));
	P->SetPhi(RanGen->Gaus(iPhi, SigmaA));

	P->SetX(RanGen->Gaus(P->X(), SigmaZt));
	P->SetY(RanGen->Gaus(P->Y(), SigmaZt));
	P->SetZ(RanGen->Gaus(P->Z(), SigmaZt));
}

double Detector::TrackLength(TRandom3 *RanGen, Particle *P)	//This should not change *P
{
	if (P->TrackIn() < 0.0 && P->Charge() != 0)
	{
		//double dx = GetElement("Vertex");	//1 mm
		double dx = 0.001;	//1 mm
		TVector3 Step(P->Direction().Unit());
		Step *= dx;		//move by dx
		TVector3 Pos(P->Position());
		
		double LengthIn = 0.0;
		double LengthOut = 0.0;
		double dE = P->E();
		double iM = P->M();

		double Beta = P->GetP4().Beta();
		double Gamma = P->GetP4().Gamma();
		double Ran = 1, Exp = 0;

		while ( IsDetectable(P->Pdg(), P->Charge(), dE-iM) &&
			RanGen->Rndm() < exp(-dx/(Const::fC * Gamma * Beta * P->Tau())))	//this should quit when particle decays too!
		{
			Beta = sqrt(1-iM*iM/dE/dE);
			Gamma = dE/iM;

			if (IsInside(Pos))
			{
				dE -= dx * EnergyLoss(1, Beta, iM);	//LAr	argon inside
				LengthIn += dx;
			}
			else
			{
				dE -= dx * EnergyLoss(3, Beta, iM);	//Fe	iron outside
				LengthOut += dx;
			}
			Pos += Step;
		}
	
		P->SetTrackIn(LengthIn);
		P->SetTrackOut(LengthOut);
		return LengthIn;
	}
	else
		return P->TrackIn();
}

double Detector::InteractionLength(int Target)
{
	if (Target == 1)	//Liquid Argon
		return Kine::RadiationLength(1.3945, 18, 40);
	else if (Target == 3)	//Steel
		return Kine::RadiationLength(7.874, 26, 56);
	else return 0.0;
}

double Detector::EnergyLoss(int Target, double Beta, double Mass)
{
	if (Mass > 10*Const::fMElectron) 
	{
		if (Target == 1)
			return Kine::Bethe(Beta, Mass, 1.3945, 188.0, 18, 40);	//LAr, code 1
		else if (Target == 3)
			return Kine::Bethe(Beta, Mass, 7.874, 286.0, 26, 56);	//Fe, code 3
		else return 0.0;
	}
	else
		return Mass/sqrt(1-Beta*Beta)/InteractionLength(Target);	//first bit is energy
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

	if (P->Pdg() == 11 || P->Pdg() == 22)	//electron or photon
	{
		if (P->Ekin() > GetElement("ThrGamma"))
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
	else if (P->Charge() != 0)	//hadron
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
	if (P->X() > GetXstart() && P->X() < GetXend() &&
	    P->Y() > GetYstart() && P->Y() < GetYend() &&
	    P->Z() > GetZstart() && P->Z() < GetZend() )
		return true;
	else return false;
}

bool Detector::IsInside(TVector3 &P)
{
	if (P.X() > GetXstart() && P.X() < GetXend() &&
	    P.Y() > GetYstart() && P.Y() < GetYend() &&
	    P.Z() > GetZstart() && P.Z() < GetZend() )
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
