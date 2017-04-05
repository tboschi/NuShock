#include "EventGenerator.h"

EventGenerator::EventGenerator(std::string SMConfig, std::string FluxConfig, std::string DetectorConfig)
{
	std::string Line, Key, Name;
	std::stringstream ssL;
	double Element;

	std::ifstream ConfigFile(SMConfig.c_str());
	while (std::getline(ConfigFile, Line))
	{
		if (Line[0] == '#') continue;

		ssL.str("");
		ssL.clear();
		ssL << Line;

		ssL >> Key >> Element;
		if (Key == "M_Sterile") SetMSterile(Element);
		if (Key == "U_e") SetUe(Element);
		if (Key == "U_m") SetUm(Element);
		if (Key == "U_t") SetUt(Element);
	}
	ConfigFile.close();


	TheBox = new Detector(DetectorConfig);
	TheGamma = new Decay(GetMSterile(), GetUe(), GetUm(), GetUt());
	hTotalFlux = new TH1D("htotalflux", "Sterile total flux", 100, 0, 20);

	StandardEnergy = new Flux(0);
	SterileEnergy = new Flux(0);

	GenMT = new TRandom3(0);	//19937 Mersenne Twister generator;
}

EventGenerator::~EventGenerator()
{
	delete TheFlux;
	delete TheBox;
	delete TheGamma;
	delete StandardEnergy;
	delete SterileEnergy;
	delete GenMT;
}

double EventGenerator::DrawEnergy()	//returns energy from distribution and also store components
{
	return hTotalFlux->GetRandom();
}

double EventGenerator::SetEnergy(double X)
{
	E_Sterile = X;
}

void EventGenerator::MakeFlux()
{
	if (NuMuFlux != 0)
	{
		NuMuFlux->MakeSterileFlux(GetMSterile(), M_Muon, GetUm(), GetUm());
		hTotalFlux->Add(NuMuFlux->GetHist());
	}
	if (NuMu_Flux != 0)
	{
		NuMu_Flux->MakeSterileFlux(GetMSterile(), M_Muon, GetUm(), GetUm());
		hTotalFlux->Add(NuMuFlux->GetHist());
	}
	if (NuEFlux != 0)
	{
		NuEFlux->MakeSterileFlux(GetMSterile(), M_Electron, GetUe(), GetUm());
		hTotalFlux->Add(NuMuFlux->GetHist());
	}
	if (NuE_Flux != 0)
	{
		NuE_Flux->MakeSterileFlux(GetMSterile(), M_Electron, GetUe(), GetUm());
		hTotalFlux->Add(NuMuFlux->GetHist());
	}
}

double EventGenerator::Probability(std::string Channel, double ESterile)
{
	double Length = Const::fM2GeV * TheBox->GetElement("Baseline");
	double Lambda = Const::fM2GeV * TheBox->GetElement("Length");
	double Ratio = TheGamma->Branch(Channel); 
	double Total = TheGamma->Total();
	double Lorentz = M_Sterile/sqrt(ESterile*ESterile - M_Sterile*M_Sterile);
	return exp(-Total * Length / Lorentz) * (1-exp(- Total * Lambda / Lorentz)) * Ratio;
}

int EventGenerator::Detectable(std::string Channel, double ESterile)
{
	double Num = GenMT->Rndm();
	if (Num <= Probability(Channel, ESterile))
		return 1;
	else return 0;
}

double EventGenerator::RandomDetectionEvent(std::string Channel)
{
	double Energy = DrawEnergy();
	return Energy*Detectable(Channel, Energy)*TheBox->Efficiency(Channel, Energy);
}

double EventGenerator::NumberOfDetected(std::string Channel)
{
	double TotNumber = 0;
	for (int i = 1; i <= hTotalflux->GetNbinsX(); ++i)
	{
		double E = hTotalFlux->GetBinContent(i);
		TotNumber += E*Probability(Channel, E)*TheBox->Efficiency(Channel, E);
	}
	return TotNumber;
}

double EventGenerator::RandomEvent()
{
	double Energy = DrawEnergy();
}

std::string EventGenerator::RandomChannel(double Energy)
{
	double Num = GenMT->Rndm();
	double Sum = 0;
	std::string Channel;
	std::vector<std::string> vChan = TheGamma->ListChannels();

	for (int i = 0; i < vChan.size(); ++i)
	{
		if (vChan.at(i) == "ALL") continue;

		Sum += TheGamma->Branch(vChan.at(i));
		if (Num <= Sum)
		{
			Channel = vChan.at(i);
			break;
		}
	}
	return Channel;
}

//Get functions
double EventGenerator::GetMSterile()
{
	return M_Sterile;
}

double EventGenerator::GetUe()
{
	return U_e;
}

double EventGenerator::GetUm()
{
	return U_m;
}

double EventGenerator::GetUt()
{
	return U_t;
}

//Set functions
void EventGenerator::SetMSterile(double X)
{
	M_Sterile = X;
}

void EventGenerator::SetUe(double X)
{
	U_e = X;
}

void EventGenerator::SetUm(double X)
{
	U_m = X;
}

void EventGenerator::SetUt(double X)
{
	U_t = X;
}
