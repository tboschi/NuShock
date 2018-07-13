#include "Detector/Detector.h"

Detector::Detector(std::string ConfigName) :
	GenMT(new TRandom3(0))
{
	std::ifstream ConfigFile(ConfigName.c_str());

	std::string Line, Key, Name, Channel;
	std::stringstream ssL;
	double Element;

	while (std::getline(ConfigFile, Line))
	{
		if (Line[0] == '#') continue;

		ssL.str("");
		ssL.clear();
		ssL << Line;

		if (ssL >> Key)
		{
			if (Key == "E")			//electon channels
			{
				ssL >> Channel >> Name;
				mEfficiencyE[Channel] = Name;
			}
			else if (Key == "M")
			{
				ssL >> Channel >> Name;
				mEfficiencyM[Channel] = Name;
			}
			else if (Key == "T")
			{
				ssL >> Channel >> Name;
				mEfficiencyT[Channel] = Name;
			}
			else if (Key.find("Target") != std::string::npos && ssL >> Name)
				mMaterial[Key] = FindMaterial(Name);
			else if (ssL >> Element)
				mDetector[Key] = Element;
		}
	}

	FuncFile = new TFile();
}

std::vector<std::string> Detector::ListKey()
{
	std::vector<std::string> List;
	std::map<std::string, double>::iterator it;
	for (it = mDetector.begin(); it != mDetector.end(); ++it)
		List.push_back(it->first);
	return List;
}

std::vector<std::string> Detector::ListChannel()
{
	std::vector<std::string> List;
	std::map<std::string, std::string>::iterator it;

	for (it = mEfficiencyE.begin(); it != mEfficiencyE.end(); ++it)
		List.push_back(it->first);

	for (it = mEfficiencyM.begin(); it != mEfficiencyM.end(); ++it)
		List.push_back(it->first);

	for (it = mEfficiencyT.begin(); it != mEfficiencyT.end(); ++it)
		List.push_back(it->first);

	return List;
}

double Detector::Get(std::string Key)
{
	return mDetector[Key];
}

Detector::Material Detector::GetMaterial(std::string Key)
{
	return mMaterial[Key];
}

Detector::Material Detector::FindMaterial(std::string Name)
{
	if (Name == "LAr" || Name == "LiquidArgon")
		return LAr;
	else if (Name == "GasAr" || Name == "GasseousArgon")
		return GasAr;
	else if (Name == "Fe" || Name == "Iron")
		return Fe;
}

double Detector::Efficiency(Neutrino *Nu)
{
	return Efficiency(Nu->Energy(), Nu->Mass());
}

double Detector::Efficiency(double Energy, double Mass)
{
	if (hhFunc)
	{
		int Ebin = hhFunc->GetXaxis()->FindBin(Energy);
		int Mbin = hhFunc->GetYaxis()->FindBin(Mass);

		return hhFunc->GetBinContent(Ebin, Mbin);
	}
	else
		return -1.0;
}

void Detector::SetEfficiency(std::string Channel, Coupling U)
{
	if (FuncFile != 0 && FuncFile->IsOpen())
		FuncFile->Close();

	switch (U)
	{
		case E:
			FuncFile = new TFile(mEfficiencyE[Channel].c_str(), "OPEN");
			hhFunc = dynamic_cast<TH2D*> (FuncFile->Get("hhfunc"));
			break;
		case M:
			FuncFile = new TFile(mEfficiencyM[Channel].c_str(), "OPEN");
			hhFunc = dynamic_cast<TH2D*> (FuncFile->Get("hhfunc"));
			break;
		case T:
			FuncFile = new TFile(mEfficiencyT[Channel].c_str(), "OPEN");
			hhFunc = dynamic_cast<TH2D*> (FuncFile->Get("hhfunc"));
			break;
		default:
			hhFunc = 0;
			break;
	}
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

double Detector::XsizeLAr()
{
	double F = pow(Get("Fiducial"), 1.0/3.0);
	return Get("WidthLAr") * (F > 0.0 ?  F : 1.0);
}

double Detector::XstartLAr()
{
	double F = pow(Get("Fiducial"), 1.0/3.0);
	return - 0.5 * Get("WidthLAr") * (F > 0.0 ? F : 1.0);
}

double Detector::XendLAr()
{
	return XsizeLAr() + XstartLAr();
}

double Detector::YsizeLAr()
{
	double F = pow(Get("Fiducial"), 1.0/3.0);
	return Get("HeightLAr") * (F > 0.0 ? F : 1.0);
}

double Detector::YstartLAr()
{
	double F = pow(Get("Fiducial"), 1.0/3.0);
	return - 0.5 * Get("HeightLAr") * (F > 0.0 ? F : 1.0);
}

double Detector::YendLAr()
{
	return YsizeLAr() + YstartLAr();
}

double Detector::ZsizeLAr()
{
	double F = pow(Get("Fiducial"), 1.0/3.0);
	return Get("LengthLAr") * (F > 0.0 ? F : 1.0);
}

double Detector::ZstartLAr()
{
	double F = (1 - pow(Get("Fiducial"), 1.0/3.0)) / 2.0;
	return F > 0.0 ? F * Get("LengthFGT") : 0.0;
}

double Detector::ZendLAr()
{
	return ZsizeLAr() + ZstartLAr();
}

double Detector::XsizeFGT()
{
	double F = pow(Get("Fiducial"), 1.0/3.0);
	return Get("WidthFGT") * (F > 0.0 ?  F : 1.0);
}

double Detector::XstartFGT()
{
	double F = pow(Get("Fiducial"), 1.0/3.0);
	return - 0.5 * Get("WidthFGT") * (F > 0.0 ? F : 1.0);
}

double Detector::XendFGT()
{
	return XsizeFGT() + XstartFGT();
}

double Detector::YsizeFGT()
{
	double F = pow(Get("Fiducial"), 1.0/3.0);
	return Get("HeightFGT") * (F > 0.0 ? F : 1.0);
}

double Detector::YstartFGT()
{
	double F = pow(Get("Fiducial"), 1.0/3.0);
	return - 0.5 * Get("HeightFGT") * (F > 0.0 ? F : 1.0);
}

double Detector::YendFGT()
{
	return YsizeFGT() + YstartFGT();
}

double Detector::ZsizeFGT()
{
	double F = pow(Get("Fiducial"), 1.0/3.0);
	return Get("LengthFGT") * (F > 0.0 ? F : 1.0);
}

double Detector::ZstartFGT()
{
	double F = (1 - pow(Get("Fiducial"), 1.0/3.0)) / 2.0;
	return ZstartLAr() + (F > 0.0 ? F * Get("LengthFGT") : 0.0);
}

double Detector::ZendFGT()
{
	return ZsizeFGT() + ZstartFGT();
}

double Detector::Xsize()
{
	return std::max(XsizeLAr(), XsizeFGT());
}

double Detector::Ysize()
{
	return std::max(YsizeLAr(), YsizeFGT());
}

double Detector::Zstart()
{
	return ZstartLAr();
}

double Detector::Zend()
{
	return ZendFGT();
}

double Detector::Zsize()
{
	return ZsizeLAr() + ZsizeFGT();
}

double Detector::AreaLAr()
{
	return XsizeLAr() * YsizeLAr();
}

double Detector::AreaFGT()
{
	return XsizeFGT() * YsizeFGT();
}

double Detector::Area()
{
	return std::max(AreaLAr(), AreaFGT());
}

bool Detector::IsInside(Particle *P)
{
	if (P->Z() > ZstartLAr() && P->Z() < ZendLAr())
		return IsInsideLAr(P);
	else if (P->Z() > ZstartFGT() && P->Z() < ZendFGT())
		return IsInsideFGT(P);
}

bool Detector::IsInsideLAr(Particle *P)
{
	return (P->Z() > ZstartLAr() && P->Z() < ZendLAr() &&
		P->Y() > YstartLAr() && P->Y() < YendLAr() &&
		P->X() > XstartLAr() && P->X() < XendLAr());
}

bool Detector::IsInsideFGT(Particle *P)
{
	return (P->Z() > ZstartFGT() && P->Z() < ZendFGT() &&
		P->Y() > YstartFGT() && P->Y() < YendFGT() &&
		P->X() > XstartFGT() && P->X() < XendFGT());
}

//special, neutrino class required
double Detector::DecayProb(Neutrino *Nu)
{
	return DecayProb(Nu, Nu->DecayTotal(), Nu->DecayBranch());
}

double Detector::DecayProb(Particle *P, double Total, double Branch)	//reaching the detector and decaying
{
	if (P->EnergyKin() < 0.0)
		return 0.0;
	else
	{
		double Length = Const::fM2GeV * Get("Baseline");
		double Lambda = Const::fM2GeV * Zsize();
		double Lorentz = P->Beta() * P->Gamma();

		return exp(- Total * Length / Lorentz) * (1 - exp(- Total * Lambda / Lorentz)) * Branch;
	}
}

