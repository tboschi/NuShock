#include "EventGenerator.h"

EventGenerator::EventGenerator(std::string SMConfig, std::string FluxConfig, std::string DetectorConfig)
{
	std::ifstream ConfigFile(SMConfig.c_str());
	std::string Line, Key;
	double Element;
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

	TheFlux = new FluxDriver(FluxConfig)
	TheBox = new Detector(DetectorConfig);
	TheGamma = new DecayRates(GetMSterile(), GetUe(), GetUm(), GetUt());

	StandardEnergy = Flux(0);
	SterileEnergy = Flux(0);

	GenMT = new TRandom3(0);	//19937 Mersenne Twister generator;
}

double EventGenerator::DrawEnergy()
{
	TheFlux->SampleEnergy(StandardEnergy, SterileEnergy);
}

double EventGenerator::SetEnergy(double X)	//this will draw energy from distribution
{
	E_Sterile = X;
}

double EventGenerator::Probability(std::string Channel, double ESterile)
{
	double Length = Tools::Const::fM2GeV * TheBox->GetElement("Baseline");
	double Lambda = Tools::Const::fM2GeV * TheBox->GetElement("Length");
	double Ratio = TheGamma->Branch(Channel); 
	double Total = TheGamma->Total();
	double Lorentz = M_Sterile/sqrt(ESterile*ESterile - M_Sterile*M_Sterile);
	return exp(- Total * Length / Lorentz) * (1-exp(- Total * Lambda / Lorentz)) * Ratio;
}

double EventGenerator::RandomDetectionEvent(std::string Channel)
{
	double Energy = DrawEnergy();
	return DrawEnergy()*Probability(Channel, Energy)*TheBox->Efficiency(Channel, Energy);
}

double EventGenerator::NumberOfDetected(std::string Channel)
{
	double TotNumber = 0;
	for (double E = 0; E < TheFlux->GetMaxEnergy(); E += 0.1)
		TotNumber += E*Probability(E)*TheBox->Efficiency(Channel, E);
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
		Sum += TheGamma->Branch(vChan.at(i));
		if (Num <= Sum)
		{
			Decay = vChan.at(i);
			break;
		}
	}
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
