#include "Detector.h"

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
			if (Key == "D")			//dirac
			{
				ssL >> Channel >> Name;
				mEfficiencyD[Channel] = Name;
			}
			else if (Key == "M")		//majorana
			{
				ssL >> Channel >> Name;
				mEfficiencyM[Channel] = Name;
			}
			else if (Key.find("Target") != std::string::npos && ssL >> Name)
				mMaterial[Key] = FindMaterial(Name);
			else if (ssL >> Element)
				mDetector[Key] = Element;
		}
	}

	FuncFile = new TFile();
	EffSet = false;
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

	for (it = mEfficiencyD.begin(); it != mEfficiencyD.end(); ++it)
		List.push_back(it->first);

	for (it = mEfficiencyM.begin(); it != mEfficiencyM.end(); ++it)
		List.push_back(it->first);

	return List;
}

double Detector::Get(std::string key)
{
	if (mDetector.count(key))
		return mDetector[key];
	else
		return 0.0;
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

double Detector::Efficiency(const Neutrino &Nu)
{
	if (EffSet)
		return Efficiency(Nu.Energy(), Nu.Mass());
	else
		return 1.0;
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

void Detector::SetEfficiency(std::string Channel, bool U)
{
	if (FuncFile != 0 && FuncFile->IsOpen())
		FuncFile->Close();

	if (U)
	{
		FuncFile = new TFile(mEfficiencyD[Channel].c_str(), "OPEN");
		hhFunc = dynamic_cast<TH2D*> (FuncFile->Get("hhfunc"));
	}
	else
	{
		FuncFile = new TFile(mEfficiencyM[Channel].c_str(), "OPEN");
		hhFunc = dynamic_cast<TH2D*> (FuncFile->Get("hhfunc"));
	}

	EffSet = true;
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
	double F = pow(Get("FiducialLAr"), 1.0/3.0);
	return Get("WidthLAr") * (F > 0.0 ?  F : 1.0);
}

double Detector::XstartLAr()
{
	double F = pow(Get("FiducialLAr"), 1.0/3.0);
	return - 0.5 * Get("WidthLAr") * (F > 0.0 ? F : 1.0);
}

double Detector::XendLAr()
{
	return XsizeLAr() + XstartLAr();
}

double Detector::YsizeLAr()
{
	double F = pow(Get("FiducialLAr"), 1.0/3.0);
	return Get("HeightLAr") * (F > 0.0 ? F : 1.0);
}

double Detector::YstartLAr()
{
	double F = pow(Get("FiducialLAr"), 1.0/3.0);
	return - 0.5 * Get("HeightLAr") * (F > 0.0 ? F : 1.0);
}

double Detector::YendLAr()
{
	return YsizeLAr() + YstartLAr();
}

double Detector::ZsizeLAr()
{
	double F = pow(Get("FiducialLAr"), 1.0/3.0);
	return Get("LengthLAr") * (F > 0.0 ? F : 1.0);
}

double Detector::ZstartLAr()
{
	double F = pow(Get("FiducialLAr"), 1.0/3.0);
	return - Get("LengthLAr") * (F > 0.0 ? (1 + F) / 2.0 : 1.0);
}

double Detector::ZendLAr()
{
	return ZsizeLAr() + ZstartLAr();
}

double Detector::XsizeFGT()
{
	double F = pow(Get("FiducialFGT"), 1.0/3.0);
	return Get("WidthFGT") * (F > 0.0 ?  F : 1.0);
}

double Detector::XstartFGT()
{
	double F = pow(Get("FiducialFGT"), 1.0/3.0);
	return - 0.5 * Get("WidthFGT") * (F > 0.0 ? F : 1.0);
}

double Detector::XendFGT()
{
	return XsizeFGT() + XstartFGT();
}

double Detector::YsizeFGT()
{
	double F = pow(Get("FiducialFGT"), 1.0/3.0);
	return Get("HeightFGT") * (F > 0.0 ? F : 1.0);
}

double Detector::YstartFGT()
{
	double F = pow(Get("FiducialFGT"), 1.0/3.0);
	return - 0.5 * Get("HeightFGT") * (F > 0.0 ? F : 1.0);
}

double Detector::YendFGT()
{
	return YsizeFGT() + YstartFGT();
}

double Detector::ZsizeFGT()
{
	double F = pow(Get("FiducialFGT"), 1.0/3.0);
	return Get("LengthFGT") * (F > 0.0 ? F : 1.0);
}

double Detector::ZstartFGT()
{
	double F = pow(Get("FiducialFGT"), 1.0/3.0);
	return Get("LengthFGT") * (F > 0.0 ? (1 - F) / 2.0 : 1.0);
}

double Detector::ZendFGT()
{
	return ZsizeFGT() + ZstartFGT();
}

double Detector::Xsize()
{
	return std::min(XsizeLAr(), XsizeFGT());
}

double Detector::Ysize()
{
	return std::min(YsizeLAr(), YsizeFGT());
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

double Detector::Radius()
{
	return sqrt(Area() / Const::pi);
}

double Detector::AngularAcceptance()
{
	return atan2(Radius(), Get("Baseline"));
}

bool Detector::IsInside(const Particle &P)
{
	if (P.Z() > ZstartLAr() && P.Z() < ZendLAr())
		return IsInsideLAr(P);
	else if (P.Z() > ZstartFGT() && P.Z() < ZendFGT())
		return IsInsideFGT(P);
}

bool Detector::IsInsideLAr(const Particle &P)
{
	return (P.Z() > ZstartLAr() && P.Z() < ZendLAr() &&
		P.Y() > YstartLAr() && P.Y() < YendLAr() &&
		P.X() > XstartLAr() && P.X() < XendLAr());
}

bool Detector::IsInsideFGT(const Particle &P)
{
	return (P.Z() > ZstartFGT() && P.Z() < ZendFGT() &&
		P.Y() > YstartFGT() && P.Y() < YendFGT() &&
		P.X() > XstartFGT() && P.X() < XendFGT());
}

double Detector::DecayProb(Neutrino &Nu)
{
	return DecayProb(Nu, Nu.DecayTotal(), Nu.DecayBranch());
}

double Detector::DecayProb(const Particle &P, double Total, double Branch)	//reaching the detector and decaying
{
	if (P.EnergyKin() < 0.0)
		return 0.0;
	else if (fabs(P.Beta() - 1.0) < 1e-9)
		return 1.0;
	else
	{
		double Length = Const::M2GeV * Get("Baseline");
		double Lambda = Const::M2GeV * Zsize();
		double Lorentz = P.Beta() * P.Gamma();

		//std::cout << "Length " << Length << "\tLambda " << Lambda << "\tLorentz " << Lorentz << "\tTotal " << Total << "\tBranch " << Branch << std::endl;

		return exp(- Total * Length / Lorentz) * (1 - exp(- Total * Lambda / Lorentz)) * Branch;
	}
}
