#include "EventGenerator.h"

EventGenerator::EventGenerator(std::string SMConfig, std::string DetectorConfig, std::string FluxConfig)
//	M_Electron(Const::fMElectron),
//	M_Muon(Const::fMMuon),
//	M_Pion(Const::fMPion),        	
//	M_Kaon(Const::fMKaon)
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

/*
EventGenerator::~EventGenerator()
{
	delete TheBox;
	delete TheGamma;
	delete TheFlux;

	delete GenMT;
}
*/

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


//MC procedures to detection of event after sampling energy
std::string EventGenerator::RandomChannel()	//First step: define decay mode
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

double EventGenerator::Probability(std::string Channel)	//reaching the detector and decaying
{							//using sampled energy
	double Length = Const::fM2GeV * TheBox->GetElement("Baseline");
	double Lambda = Const::fM2GeV * TheBox->GetElement("Length");
	double Ratio = TheGamma->Branch(Channel); 
	double Total = TheGamma->Total();
	double Lorentz = GetMass()/sqrt(GetEnergy()*GetEnergy() - GetMass()*GetMass());
	return exp(-Total * Length * Lorentz) * (1-exp(- Total * Lambda * Lorentz)) * Ratio;
}

bool EventGenerator::Detectable(std::string Channel)	//Second step: can I detect the decay?
{
	if (Channel == "R")
		return (GenMT->Rndm() <= Probability(RandomChannel()));
	else return (GenMT->Rndm() <= Probability(Channel));
}

bool EventGenerator::RandomDetectionEvent(std::string Channel)	//Third step: do I have enough efficiency?
{
	if (Detectable(Channel))
	{
		if (Channel == "R")
			return (GenMT->Rndm() <= TheBox->Efficiency(RandomChannel(), GetEnergy()));
		else return (GenMT->Rndm() <= TheBox->Efficiency(Channel, GetEnergy()));
	}
	else return false;
}

int EventGenerator::SimulateDecay(std::string Channel)	//Fourth step: simulate the phase space of the decay
{
	TLorentzVector N_vec(0, 0, GetMomentum(), GetEnergy());		//Heavy neutrino is along z-axis

	double Weight;
	int Return;
	TheGamma->SetNvec(N_vec);
	if (Channel == "R")
		Return = TheGamma->GetPhaseSpace(RandomChannel(), Weight);
	else Return = TheGamma->GetPhaseSpace(Channel, Weight);
	return Return;
}

TLorentzVector *EventGenerator::GetDecayProduct(int i)
{
	return TheGamma->GetDecayProduct(i);
}

//Flux as PDF for MC
void EventGenerator::MakeSterileFlux()	//Generate the flux for heavy neutrinos
{
	TheFlux->SetBaseline(TheBox->GetElement("Baseline"));
	TheFlux->MakeSterileFlux(GetMass(), GetUe(), GetUm(), GetUt());
}

void EventGenerator::MakeStandardFlux()	//Generate the flux of SM neutrinos, just for comparison
{
	TheFlux->SetBaseline(TheBox->GetElement("Baseline"));
	TheFlux->MakeStandardFlux();
}

double EventGenerator::SampleEnergy()	//Sample Energy according to PDF distribution
{
	double Energy = TheFlux->SampleEnergy();
	SetEnergy(Energy);
	return Energy;
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
