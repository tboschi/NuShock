#include "EventGenerator.h"

EventGenerator::EventGenerator(std::string SMConfig, std::string DetectorConfig, std::string FluxConfig)
//	M_Electron(Const::fMElectron),
//	M_Muon(Const::fMMuon),
//	M_Pion(Const::fMPion),        	
//	M_Kaon(Const::fMKaon)
{
	TheBox = new Detector(DetectorConfig);
	TheGamma = new Decay();
	TheFlux = new FluxDriver(FluxConfig);

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

	GenMT = new TRandom3(0);	//19937 Mersenne Twister generator

	SetChannel();			//Channel is initialised randomly
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
		{
			++i;
			break;
		}
	}
	return vChan.at(--i);
}

double EventGenerator::EventProbability()	//reaching the detector and decaying
{						//using sampled energy
	if (GetEnergy() < GetMass())
		return 0.0;
	else
	{
		double Length = Const::fM2GeV * TheBox->GetElement("Baseline");
		double Lambda = Const::fM2GeV * TheBox->GetElement("Length");
		double Ratio = TheGamma->Branch(GetChannel()); 
		double Total = TheGamma->Total();
		double Lorentz = GetMass()/sqrt(GetEnergy()*GetEnergy() - GetMass()*GetMass());
		return exp(- Total * Length * Lorentz) * (1 - exp(- Total * Lambda * Lorentz)) * Ratio;
	}
}

bool EventGenerator::EventInDetector()		//Second step: is the decay inside the detector?
{
	return (GenMT->Rndm() <= EventProbability());
}

bool EventGenerator::EventDetectable()	//Third step: is the detector able to detect it?
{
	return (GenMT->Rndm() <= TheBox->Efficiency(GetChannel(), GetEnergy()));
}

int EventGenerator::EventKinematics()	//Fourth step: simulate the phase space of the decay
{
	TLorentzVector N_vec(0, 0, GetMomentum(), GetEnergy());		//Heavy neutrino is along z-axis

	double Weight;
	int Return;
	TheGamma->SetNvec(N_vec);

	return TheGamma->PhaseSpace(GetChannel(), Weight);
}

TLorentzVector *EventGenerator::GetDecayProduct(int i)
{
	return TheGamma->GetDecayProduct(i);
}

//Flux as PDF for MC
void EventGenerator::MakeSterileFlux()	//Generate the flux for heavy neutrinos
{
	if (IsChanged())
	{
		TheFlux->MakeSterileFlux(GetMass(), GetUe(), GetUm(), GetUt());
		TheFlux->SetBaseline(TheBox->GetElement("Baseline"));
		TheFlux->SetPOT(TheBox->GetElement("POT"));
	}
}

void EventGenerator::MakeStandardFlux()	//Generate the flux of SM neutrinos, just for comparison
{
	TheFlux->MakeStandardFlux();
	TheFlux->SetBaseline(TheBox->GetElement("Baseline"));
}

double EventGenerator::SampleEnergy()	//Sample Energy according to PDF distribution
{
	double Energy = 0.0;
	while (Energy < GetMass())
		Energy = TheFlux->SampleEnergy();

	SetEnergy(Energy);
	return Energy;
}

double EventGenerator::FluxIntensity()	//Get the flux intensity at given energy
{
	return TheFlux->GetIntensity(GetEnergy());
}

//Get functions

bool EventGenerator::IsChanged()
{
	return ( M_Sterile != M_Sterile_prev || 
		 U_e != U_e_prev ||
		 U_m != U_m_prev ||
		 U_t != U_t_prev );
}

std::string EventGenerator::GetChannel()
{
	return sChannel;
}

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
void EventGenerator::SetChannel(std::string Ch)
{
	if (Ch == "R")
		sChannel.assign(RandomChannel());
	else sChannel.assign(Ch); 
}

void EventGenerator::SetMass(double X)
{
	M_Sterile_prev = M_Sterile;
	M_Sterile = X;
	TheGamma->SetMass(X);
}

void EventGenerator::SetEnergy(double X)
{
	E_Sterile_prev = E_Sterile;
	E_Sterile = X;
}

void EventGenerator::SetUe(double X)
{
	U_e_prev = U_e;
	U_e = X;
	TheGamma->SetUe(X);
}

void EventGenerator::SetUm(double X)
{
	U_m_prev = U_m;
	U_m = X;
	TheGamma->SetUm(X);
}

void EventGenerator::SetUt(double X)
{
	U_t_prev = U_t;
	U_t = X;
	TheGamma->SetUt(X);
}
