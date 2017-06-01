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

double EventGenerator::EventProbability()	//reaching the detector and decaying
{						//using sampled energy
	if (GetEnergy() < GetMass())
		return 0.0;
	else
	{
		double Total = TheGamma->Total();
		double Ratio = TheGamma->Branch(GetChannel()); 
		double Length = Const::fM2GeV * TheBox->GetElement("Baseline");
		double Lambda = Const::fM2GeV * TheBox->GetElement("Length");
		double Lorentz = sqrt(GetEnergy(2)/GetMass(2) - 1.0);
		return exp(- Total * Length / Lorentz) * (1 - exp(- Total * Lambda / Lorentz)) * Ratio;
	}
}

double EventGenerator::EventEfficiency(double Efficiency)
{
	if (Efficiency < 0.0)
		return TheBox->Efficiency(GetChannel(), GetEnergy());
	else return Efficiency;
}

//modifications:
//hMod, hProb
//even in .h file
double EventGenerator::EventTotalNumber(TH1D* hMod, TH1D* hProb, double Efficiency)
{
	double A = TheFlux->GetStartRange();
	double B = TheFlux->GetEndRange();
	double EnStep = (B-A)/(TheFlux->GetBinNumber()-1);

	double Total = 0;
	for (double Energy = A; Energy < B; Energy += EnStep)
	{
		SetEnergy(Energy);
		if (Efficiency < 0.0)
		{
			//hProb->Fill(Energy+1e-6, EventProbability());
			//hMod->Fill(Energy+1e-6, FluxIntensity() * EventProbability());
			Total += FluxIntensity() * EventProbability();
		}
		else Total += FluxIntensity() * EventProbability() * EventEfficiency(Efficiency);
	}
	return Total * EnStep;
}

//Random generator

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

bool EventGenerator::EventInDetector()		//Second step: is the decay inside the detector?
{
	return (GenMT->Rndm() <= EventProbability());
}

bool EventGenerator::EventDetectable()	//Third step: is the detector able to detect it?
{
	return (GenMT->Rndm() <= EventEfficiency());
}

//Kinematics
int EventGenerator::EventKinematics()	//Fourth step: simulate the phase space of the decay
{
	TLorentzVector N_vec(0, 0, 0, GetMass());	//Rest frame fro the heavy neutrino

	double Weight;
	int Return;
	TheGamma->SetNvec(N_vec);

	double Products = TheGamma->PhaseSpace(GetChannel(), Weight);
	if (GenMT->Rndm() <= Weight)
		return Products;
	else return 0.0;
}

TLorentzVector *EventGenerator::GetDecayProduct(int i, bool Smear)
{
	TLorentzVector * Vi = TheGamma->GetDecayProduct(i);

	double beta = GetMomentum()/GetEnergy();
	Vi->Boost(0,0,beta);	//boost along z axis

	if (Smear)
		SmearVector(Vi);

	return Vi;
}


//Flux as PDF for MC
//modifications
//baseline
//pot
//area
void EventGenerator::MakeSterileFlux(bool TotalPOT)	//Generate the flux for heavy neutrinos
{
	if (TheFlux->MakeSterileFlux(GetMass(), GetUe(), GetUm(), GetUt()))
	{
		TheFlux->SetBaseline(TheBox->GetElement("Baseline"));
		//TheFlux->SetBaseline(1/1e4);

		double Y = TotalPOT ? 1.0e7 * TheBox->GetElement("Years") : 1.0;
		TheFlux->SetPOT(Y * TheBox->GetElement("POT/s"));
		//TheFlux->SetPOT(1.0e21);

		TheFlux->SetArea(TheBox->GetElement("Height")*TheBox->GetElement("Width")*1.0e4);
	}
}

/*
void EventGenerator::MakeStandardFlux(bool TotalPOT)	//Generate the flux of SM neutrinos, just for comparison
{
	if (IsChanged())
	{
		TheFlux->MakeStandardFlux();
		TheFlux->SetBaseline(TheBox->GetElement("Baseline"));

		double Y = TotalPOT ? 1.0e7 * TheBox->GetElement("Years") : 1.0;
		TheFlux->SetPOT(Y * TheBox->GetElement("POT/s"));

		TheFlux->SetArea(TheBox->GetElement("Height")*TheBox->GetElement("Width")*1.0e4);
	}
}
*/

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

//Statistics
void EventGenerator::SmearVector(TLorentzVector* N)
{
	double iE = N->E();
	double iP = N->P();
	double iTheta = N->Theta();
	double iPhi = N->Phi();

	double SigmaE = sqrt(iE) * TheBox->GetElement("Energy");
	double SigmaA = TheBox->GetElement("Angle");

	double oE = GenMT->Gaus(iE, SigmaE);
	double oP = sqrt(oE*oE - iE*iE + iP*iP);
	double oTheta = GenMT->Gaus(iTheta, SigmaA);
	double oPhi = GenMT->Gaus(iPhi, SigmaA);
	
	N->SetE(oE);
	N->SetPx(oP * cos(oPhi) * sin(oTheta));
	N->SetPy(oP * sin(oPhi) * sin(oTheta));
	N->SetPz(oP * cos(oTheta));
}

//Get functions


std::string EventGenerator::GetChannel()
{
	return sChannel;
}

double EventGenerator::GetMass(int Pow)
{
	return pow(M_Sterile, Pow);
}

double EventGenerator::GetEnergy(int Pow)
{
	return pow(E_Sterile, Pow);
}

double EventGenerator::GetMomentum(int Pow)
{
	double P = sqrt(GetEnergy(2) - GetMass(2));
	return pow(P, Pow);
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
	M_Sterile = X;
	TheGamma->SetMass(X);
}

void EventGenerator::SetEnergy(double X)
{
	E_Sterile = X;
}

void EventGenerator::SetUe(double X)
{
	U_e = X;
	TheGamma->SetUe(X);
}

void EventGenerator::SetUm(double X)
{
	U_m = X;
	TheGamma->SetUm(X);
}

void EventGenerator::SetUt(double X)
{
	U_t = X;
	TheGamma->SetUt(X);
}
