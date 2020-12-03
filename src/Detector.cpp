#include "Detector.h"

Detector::Detector(std::string configName, std::string mod) :
	GenMT(new TRandom3(0)),
	module(mod)
{
	dc.ReadCard(configName);

	/*
	std::ifstream ConfigFile(ConfigName.c_str());

	std::string Line, Key, Name, Channel;
	std::stringstream ssL;
	double Element;

	dc.ReadCard(configName);
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
	*/

	effSet = false;
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
	double val;
	if (dc.Get(key, val))
		return val;
	else
		return 0.0;
}

Detector::Material Detector::GetMaterial(std::string key)
{
	std::string mat;
	if (dc.Get(key, mat))
		return FindMaterial(mat);
	else
		return undefined;
}

Detector::Material Detector::FindMaterial(std::string name)
{
	if (name == "LAr" || name == "LiquidArgon")
		return LAr;
	else if (name == "GasAr" || name == "GasseousArgon")
		return GasAr;
	else if (name == "Fe" || name == "Iron")
		return Fe;
	else if (name == "Pb" || name == "Lead")
		return Pb;
	else
		return undefined;
}

double Detector::Efficiency(const Neutrino &Nu)
{
	if (effSet)
		return Efficiency(Nu.Energy(), Nu.Mass(), Nu.DecayChannelName());
	else
		return 1.0;
}

double Detector::Efficiency(double energy, double mass, std::string channel)
{
	if (mhFunc.count(channel))
	{
		int eBin = mhFunc[channel]->GetXaxis()->FindBin(energy);
		int mBin = mhFunc[channel]->GetYaxis()->FindBin(mass);

		return mhFunc[channel]->GetBinContent(eBin, mBin);
	}
	else
	{
		std::cout << "efficiency for " << channel << " not set yet" << std::endl;
		return -1.0;
	}
}

//key will be a combination such as CHANNEL_MODULE_FERMION
//e.g. MPI_LAr_dirac
void Detector::SetEfficiency(std::string key, std::string channel)
{
	double rat = 1;
	if (module.empty())
	{
		if (key.find("LAr") != std::string::npos)
			rat = RatioLAr();
		else if (key.find("FGT") != std::string::npos)
			rat = RatioFGT();
	}

	std::string file;
	std::cout << "Setting key " << key << std::endl;
	if (dc.Get(key, file))
	{
		std::cout << "Found" << std::endl;
		TFile funcFile(file.c_str(), "READ");

		TH2D *hist = dynamic_cast<TH2D*> (funcFile.Get("hhfunc"));
		if (!mhFunc.count(channel))
		{
			mhFunc[channel] = dynamic_cast<TH2D*> (hist->Clone());
			mhFunc[channel]->SetDirectory(0);
			mhFunc[channel]->Scale(rat);
		}
		else	//adding module with rate
			mhFunc[channel]->Add(hist, rat);

		effSet = true;
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
	//double F = (Get("FiducialLAr"), 1.0/3.0);
	//return Get("WidthLAr") * (F > 0.0 ?  F : 1.0);
	double F = Get("FiducialLAr");
	return Get("WidthLAr") - 2*F;
}

double Detector::XstartLAr()
{
	//double F = pow(Get("FiducialLAr"), 1.0/3.0);
	//return - 0.5 * Get("WidthLAr") * (F > 0.0 ? F : 1.0);
	double F = Get("FiducialLAr");
	return - 0.5 * (Get("WidthLAr") - 2*F);
}

double Detector::XendLAr()
{
	return XsizeLAr() + XstartLAr();
}

double Detector::YsizeLAr()
{
	//double F = pow(Get("FiducialLAr"), 1.0/3.0);
	//return Get("HeightLAr") * (F > 0.0 ? F : 1.0);
	double F = Get("FiducialLAr");
	return Get("HeightLAr") - 2*F;
}

double Detector::YstartLAr()
{
	//double F = pow(Get("FiducialLAr"), 1.0/3.0);
	//return - 0.5 * Get("HeightLAr") * (F > 0.0 ? F : 1.0);
	double F = Get("FiducialLAr");
	return - 0.5 * (Get("HeightLAr") - 2*F);
}

double Detector::YendLAr()
{
	return YsizeLAr() + YstartLAr();
}

double Detector::ZsizeLAr()
{
	//double F = pow(Get("FiducialLAr"), 1.0/3.0);
	//return Get("LengthLAr") * (F > 0.0 ? F : 1.0);
	double F = Get("FiducialLAr");
	return Get("LengthLAr") - 2*F;
}

double Detector::ZstartLAr()
{
	//double F = pow(Get("FiducialLAr"), 1.0/3.0);
	//return - Get("LengthLAr") * (F > 0.0 ? (1 + F) / 2.0 : 1.0);
	double F = Get("FiducialLAr");
	return F + Get("Baseline");;
}

double Detector::ZendLAr()
{
	return ZsizeLAr() + ZstartLAr();
}

double Detector::XsizeFGT()
{
	//double F = pow(Get("FiducialFGT"), 1.0/3.0);
	//return Get("WidthFGT") * (F > 0.0 ?  F : 1.0);
	double F = Get("FiducialFGT");
	return Get("WidthFGT") - 2*F;
}

double Detector::XstartFGT()
{
	//double F = pow(Get("FiducialFGT"), 1.0/3.0);
	//return - 0.5 * Get("WidthFGT") * (F > 0.0 ? F : 1.0);
	double F = Get("FiducialFGT");
	return - 0.5 * (Get("WidthFGT") - 2*F);
}

double Detector::XendFGT()
{
	return XsizeFGT() + XstartFGT();
}

double Detector::YsizeFGT()
{
	//double F = pow(Get("FiducialFGT"), 1.0/3.0);
	//return Get("HeightFGT") * (F > 0.0 ? F : 1.0);
	double F = Get("FiducialFGT");
	//return Get("HeightFGT") - 2*F;
	return 2 * (Get("RadiusFGT") - F);
}

double Detector::YstartFGT()
{
	//double F = pow(Get("FiducialFGT"), 1.0/3.0);
	//return - 0.5 * Get("HeightFGT") * (F > 0.0 ? F : 1.0);
	double F = Get("FiducialFGT");
	return - 0.5 * (Get("HeightFGT") - 2*F);
}

double Detector::YendFGT()
{
	return YsizeFGT() + YstartFGT();
}

double Detector::YcentreFGT()
{
	return 0;
}

double Detector::ZsizeFGT()
{
	//double F = pow(Get("FiducialFGT"), 1.0/3.0);
	//return Get("LengthFGT") * (F > 0.0 ? F : 1.0);
	double F = Get("FiducialFGT");
	//return Get("LengthFGT") - 2*F;
	return 2 * (Get("RadiusFGT") - F);
}

double Detector::ZstartFGT()
{
	//double F = pow(Get("FiducialFGT"), 1.0/3.0);
	//return Get("LengthFGT") * (F > 0.0 ? (1 - F) / 2.0 : 1.0);
	double F = Get("FiducialFGT");
	return ZstartLAr() + ZsizeLAr() + F;
}

double Detector::ZendFGT()
{
	return ZsizeFGT() + ZstartFGT();
}

double Detector::ZcentreFGT()
{
	return ZstartFGT() + 0.5 * ZsizeFGT();
}

double Detector::Xstart()
{
	if (module == "LAr")
		return XstartLAr();
	else if (module == "FGT")
		return XstartFGT();
	else 
		return XstartLAr();
}

double Detector::Xend()
{
	if (module == "LAr")
		return XendLAr();
	else if (module == "FGT")
		return XendFGT();
	else 
		return XendFGT();
}

double Detector::Xsize()
{
	if (module == "LAr")
		return XsizeLAr();
	else if (module == "FGT")
		return XsizeFGT();
	else 
		return std::min(XsizeLAr(), XsizeFGT());
}

double Detector::Ystart()
{
	if (module == "LAr")
		return YstartLAr();
	else if (module == "FGT")
		return YstartFGT();
	else 
		return YstartLAr();
}

double Detector::Yend()
{
	if (module == "LAr")
		return YendLAr();
	else if (module == "FGT")
		return YendFGT();
	else 
		return YendFGT();
}

double Detector::Ysize()
{
	if (module == "LAr")
		return YsizeLAr();
	else if (module == "FGT")
		return YsizeFGT();
	else 
		return std::min(YsizeLAr(), YsizeFGT());
}

double Detector::Zstart()
{
	if (module == "LAr")
		return ZstartLAr();
	else if (module == "FGT")
		return ZstartFGT();
	else 
		return ZstartLAr();
}

double Detector::Zend()
{
	if (module == "LAr")
		return ZendLAr();
	else if (module == "FGT")
		return ZendFGT();
	else 
		return ZendFGT();
}

double Detector::Zsize()
{
	if (module == "LAr")
		return ZsizeLAr();
	else if (module == "FGT")
		return ZsizeFGT();
	else 
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
	if (module == "LAr")
		return AreaLAr();
	else if (module == "FGT")
		return AreaFGT();
	else 
		return std::max(AreaLAr(), AreaFGT());
}

double Detector::Radius()
{
	return sqrt(Area() / Const::pi);
}

double Detector::VolumeLAr()
{
	return AreaLAr() * ZsizeLAr();
}

double Detector::VolumeFGT()
{
	return AreaFGT() * ZsizeFGT() * Const::pi / 4;
}

double Detector::Volume()
{
	if (module == "LAr")
		return VolumeLAr();
	else if (module == "FGT")
		return VolumeFGT();
	else 
		return VolumeLAr() + VolumeFGT();
}

double Detector::RatioLAr()
{
	return VolumeLAr() / Volume();
}

double Detector::RatioFGT()
{
	return VolumeFGT() / Volume();
}

double Detector::Ratio()
{
	if (module == "LAr")
		return RatioLAr();
	else if (module == "FGT")
		return RatioFGT();
}

double Detector::WeightLAr()
{
	return Get("WeightLAr");
}

double Detector::WeightFGT()
{
	return Get("WeightFGT");
}

double Detector::Weight()
{
	if (module == "LAr")
		WeightLAr();
	else if (module == "FGT")
		WeightFGT();
	else
		return WeightLAr() + WeightFGT();
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
	return ( sqrt(pow(P.Z() - ZcentreFGT(), 2) +
		      pow(P.Y() - YcentreFGT(), 2)) < 0.5 * ZsizeFGT() &&
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
	else if (std::abs(P.Beta() - 1.0) < 1e-9)
		return 0.0;
	else
	{
		double Length = Const::M2GeV * Zstart();
		double Lambda = Const::M2GeV * Zsize();
		double Lorentz = P.Beta() * P.Gamma();

		//std::cout << "Length " << Length << "\tLambda " << Lambda << "\tLorentz " << Lorentz << "\tTotal " << Total << "\tBranch " << Branch << std::endl;

		return exp(- Total * Length / Lorentz) * (1 - exp(- Total * Lambda / Lorentz)) * Branch;
	}
}
