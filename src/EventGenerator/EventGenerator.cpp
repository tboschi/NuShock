#include "EventGenerator.h"

EventGenerator::EventGenerator(std::string SMConfig, std::string DetectorConfig, std::string FluxConfig) : 
	M_Electron(Const::fMElectron),
	M_Muon(Const::fMMuon),
	M_Pion(Const::fMPion),        	
	M_Kaon(Const::fMKaon)
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
		if (Key == "M_Sterile") SetMass(Element);
		if (Key == "U_e") SetUe(Element);
		if (Key == "U_m") SetUm(Element);
		if (Key == "U_t") SetUt(Element);
	}
	ConfigFile.close();

	TheBox = new Detector(DetectorConfig);
	TheGamma = new Decay(GetMass(), GetUe(), GetUm(), GetUt());
	TheFlux = new FluxDriver(FluxConfig);

	GenMT = new TRandom3(0);	//19937 Mersenne Twister generator
}

EventGenerator::~EventGenerator()
{
	delete TheBox;
	delete TheGamma;
	delete TheFlux;

	delete GenMT;
}

Detector* EventGenerator::GetDetectorPtr()
{
	return TheBox;
}

Decay* EventGenerator::GetDecayPtr()
{
	return TheGamma;
}

FluxDriver* EventGenerator::GetFluxDriverPtr()
{
	return TheFlux;
}

double EventGenerator::Probability(std::string Channel)		//of reaching the detectore and decaying inside
{								//using given energy
	double Length = Const::fM2GeV * TheBox->GetElement("Baseline");
	double Lambda = Const::fM2GeV * TheBox->GetElement("Length");
	double Ratio = TheGamma->Branch(Channel); 
	double Total = TheGamma->Total();
	double Lorentz = GetMass()/sqrt(GetEnergy()*GetEnergy() - GetMass()*GetMass());
	return exp(-Total * Length * Lorentz) * (1-exp(- Total * Lambda * Lorentz)) * Ratio;
}

bool EventGenerator::Detectable(std::string Channel)	//defaul is random
{
	if (Channel == "R")
		return (GenMT->Rndm() <= Probability(RandomChannel()));		//Return bool according distribution of probability
	else return (GenMT->Rndm() <= Probability(Channel));		//fixed channel
}

bool EventGenerator::RandomDetectionEvent(std::string Channel)	//Defaul is random channel
{
	if (Detectable(Channel))
	{
		if (Channel == "R")
			return (GenMT->Rndm() <= TheBox->Efficiency(RandomChannel(), GetEnergy()));
		else return (GenMT->Rndm() <= TheBox->Efficiency(Channel, GetEnergy()));
	}
	else return false;
}

std::string EventGenerator::RandomChannel()
{
	double Num = GenMT->Rndm();
	double Sum = 0;
	std::vector<std::string> vChan = TheGamma->ListChannels();

	int i = 0;
	for ( ; i < vChan.size(); ++i)
	{
		if (vChan.at(i) == "ALL") continue;

		Sum += TheGamma->Branch(vChan.at(i));
		if (Num <= Sum)
			break;
	}
	return vChan.at(i);
}

void EventGenerator::MakeSterileFlux()
{
	TheFlux->SetBaseline(TheBox->GetElement("Baseline"));
	TheFlux->MakeSterileFlux(GetMass(), GetUe(), GetUm(), GetUt());
}

void EventGenerator::MakeStandardFlux()
{
	TheFlux->SetBaseline(TheBox->GetElement("Baseline"));
	TheFlux->MakeStandardFlux();
}

//Real MC generation
double EventGenerator::SampleEnergy()
{
	double Energy = TheFlux->SampleEnergy();
	SetEnergy(Energy);
	return Energy;
}

void EventGenerator::Decay2Body(TLorentzVector &D1, TLoretnzVector &D2)
{
	TLorentzVector N_vec(0, 0, GetMomentum(), GetEnergy());
	std::vector<TLorentzVector> vDecay;

	TGenPhaseSpace Event;
	double Mass[2] = {D1.M(), D2.M()};
	Event.SetDecay(N_vec, 2, Mass);

	double Weight = Event.Generate();
	D1 = Event.GetDecay(0);	
	D2 = Event.GetDecay(1);	
}

void EventGenerator::Decay3Body(TLorentzVector &D1, TLoretnzVector &D2, TLoretnzVector &D3)
{
	TLorentzVector N_vec(0, 0, GetMomentum(), GetEnergy());
	std::vector<TLorentzVector> vDecay;

	TGenPhaseSpace Event;
	double Mass[3] = {D1.M(), D2.M(), D2.M()};
	Event.SetDecay(N_vec, 3, Mass);

	double Weight = Event.Generate();
	D1 = Event.GetDecay(0);	
	D2 = Event.GetDecay(1);	
	D3 = Event.GetDecay(2);	
}

//Get functions
double EventGenerator::GetMass()
{
	return M_Sterile;
}

double EventGenerator::GetEnergy()
{
	return E_Sterile;
}

double EventGenerator::GetMomentum()
{
	return sqrt(GetEnergy()*GetEnergy()-GetMass()*GetMass());
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
void EventGenerator::SetMass(double X)
{
	M_Sterile = X;
}

void EventGenerator::SetEnergy(double X)
{
	E_Sterile = X;
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
