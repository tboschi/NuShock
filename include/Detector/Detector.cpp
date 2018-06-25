#include "Detector/Detector.h"

Detector::Detector(std::string ConfigName) :
	GenMT(new TRandom3(0))
{
	std::ifstream ConfigFile(ConfigName.c_str());

	std::string Line, Key, Name, Id;
	std::stringstream ssL;
	double Element;

	while (std::getline(ConfigFile, Line))
	{
		if (Line[0] == '#') continue;

		ssL.str("");
		ssL.clear();
		ssL << Line;

		if (ssL >> Key >> Element)
		{
			std::cout << "here 0" << std::endl;
			mDetector[Key] = Element;
		}
		else if (ssL >> Key >> Name)
		{
			std::cout << "here 1" << std::endl;
			if (Key.find("Target") != std::string::npos)
				mMaterial[Key] = FindMaterial(Name);
		}
		else if (ssL >> Id >> Key >> Name)
		{
			std::cout << "here 2" << std::endl;
			if (Id == "E")			//electon channels
				mEfficiencyE[Key] = Name;
			else if (Id == "M")
				mEfficiencyM[Key] = Name;
			else if (Id == "T")
				mEfficiencyT[Key] = Name;
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
	std::cout << "found " << Name << std::endl;
	if (Name == "LAr" || Name == "LiquidArgon")
	{
		std::cout << "return LAr" << std::endl;
		return LAr;
	}
	else if (Name == "GasAr" || Name == "GasseousArgon")
	{
		std::cout << "return GasAr" << std::endl;
		return GasAr;
	}
	else if (Name == "Fe" || Name == "Iron")
	{
		std::cout << "return Fe" << std::endl;
		return Fe;
	}
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

double Detector::Xsize()
{
	double F = pow(Get("Fiducial"), 1.0/3.0);
	return Get("Width") * (F > 0.0 ?  F : 1.0);
}

double Detector::Xstart()
{
	double F = pow(Get("Fiducial"), 1.0/3.0);
	return - 0.5 * Get("Width") * (F > 0.0 ? F : 1.0);
}

double Detector::Xend()
{
	return Xsize() + Xstart();
}

double Detector::Ysize()
{
	double F = pow(Get("Fiducial"), 1.0/3.0);
	return Get("Height") * (F > 0.0 ? F : 1.0);
}

double Detector::Ystart()
{
	double F = pow(Get("Fiducial"), 1.0/3.0);
	return - 0.5 * Get("Height") * (F > 0.0 ? F : 1.0);
}

double Detector::Yend()
{
	return Ysize() + Ystart();
}

double Detector::Zsize()
{
	double F = pow(Get("Fiducial"), 1.0/3.0);
	return Get("Length") * (F > 0.0 ? F : 1.0);
}

double Detector::Zstart()
{
	double F = (1 - pow(Get("Fiducial"), 1.0/3.0)) / 2.0;
	return F > 0.0 ? F * Get("Length") : 0.0;
}

double Detector::Zend()
{
	return Zsize() + Zstart();
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
			//sqrt(pow(TheN->Energy()/TheN->Mass(), 2) - 1.0);	//betagamma, to invert

		return exp(- Total * Length / Lorentz) * (1 - exp(- Total * Lambda / Lorentz)) * Branch;
	}
}

