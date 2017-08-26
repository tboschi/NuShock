#include "EventGenerator.h"

EventGenerator::EventGenerator(std::string SMConfig, std::string DetectorConfig, std::string FluxConfig)
//	M_Electron(Const::fMElectron),
//	M_Muon(Const::fMMuon),
//	M_Pion(Const::fMPion),        	
//	M_Kaon(Const::fMKaon)
{
	TheBox = new Detector(DetectorConfig);
	TheGamma = new Decay();
	TheFlux = new FluxDriver(FluxConfig);	//I have decided this is kinetic energy
	InDetector = 0;

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

	SetChannel("ALL");			//Channel is not initialised randomly!
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

long double EventGenerator::EventProbability()	//reaching the detector and decaying
{						//using sampled energy
	if (GetEnergy() < GetMass())
		return 0.0;
	else
	{
		double Total = TheGamma->Total();
		double Ratio = TheGamma->Branch(GetChannel()); 
		double Length = Const::fM2GeV * (TheBox->GetElement("Baseline") + TheBox->GetZstart());
		double Lambda = Const::fM2GeV * TheBox->GetZsize();
		long double Lorentz = sqrt(GetEnergy(2)/GetMass(2) - 1.0);	//betagamma, to invert
		return exp(- Total * Length / Lorentz) * (1 - exp(- Total * Lambda / Lorentz)) * Ratio;
	}
}

double EventGenerator::EventEfficiency()
{
	return TheBox->Efficiency(GetChannel(), GetEnergy(), GetMass());
}

//modifications:
//hMod, hProb
//even in .h file
void EventGenerator::EventTotalNumber(bool Efficiency)
{
	double Range = TheFlux->GetRange();
	//double A = TheFlux->GetStartRange();
	//double B = TheFlux->GetEndRange();
	//double EnStep = Range/500.0;
	double EnStep = Range/100.0;

	long double Signal = 0.0L;
	long double Background = 0.0L; 
	long double ChiSquared = 0.0L; 
	long double Sgn, Bkg, Chi2;
	for (double EnKin = 0.0; EnKin < Range; EnKin += EnStep)
	{
		SetEnergyKin(EnKin);
	//	if (EventProbability() <= 0.0) continue;	//speed up?

		Sgn = EnStep * FluxIntensity() * EventProbability();
		if (Efficiency)
		{
			Sgn *= EventEfficiency();
			Bkg = BackgroundIntensity();
			Chi2 = Sgn*Sgn / Bkg;		//X2 = Sum (S^2 / B)
		}

		//std::cout << "M " << GetMass() << "\t" << Energy << "\t" << Sgn << std::endl;
		//std::cout << GetMass() << "\t" << GetEnergy() << "\t" << FluxIntensity() << "\t" << EventProbability() << std::endl;

		Signal += Sgn;
		Background += Background;
		ChiSquared += Chi2;
	}

	SetSignalNumber(Signal);
	SetBackgroundNumber(Background);
	SetReducedChi2(Chi2 / 99.0);	//not so quite correct
	//return Signal * EnStep;
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
	TLorentzVector N_vec(0, 0, 0, GetMass());	//Rest frame for the heavy neutrino

	double Weight;
	TheGamma->SetNvec(N_vec);

	double Products = TheGamma->PhaseSpace(GetChannel(), Weight);
	if (GenMT->Rndm() <= Weight)
		return Products;
	else return 0.0;
}

Particle *EventGenerator::GetDecayProduct(int i, bool Smear)
{
	int pdg;
	TLorentzVector * Vi = TheGamma->GetDecayProduct(i, pdg);
	//Particle * Pi = TheGamma->GetDecayProduct(pdg, i);

	double BetaZ = GetMomentum()/GetEnergy();
	double Radius = sqrt(pow(TheBox->GetElement("Width"),2) + pow(TheBox->GetElement("Height"),2));
	double SigmaT = atan2(Radius, TheBox->GetElement("Baseline"));	//3sigma will be inside the detector (better distribution needed)

	TVector3 Parent(0, 0, GetMomentum()/GetEnergy());
	Parent.SetTheta(abs(GenMT->Gaus(0, SigmaT)));	//abs needed to avoid degeneracy with phi (?)
	Parent.SetPhi(GenMT->Uniform(-Const::fPi, Const::fPi));

	Vi->Boost(Parent);	//boost along z axis

	double PosX = GenMT->Uniform(TheBox->GetXsize());
	double PosY = GenMT->Uniform(TheBox->GetYsize());
	double PosZ = GenMT->Uniform(TheBox->GetZsize());
	TVector3 Pos(PosX, PosY, PosZ);
	Particle *P = new Particle(pdg, *Vi, Pos);	//mainly charged particles

	if (Smear)
		TheBox->SignalSmearing(GenMT, P);

	return P;
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

void EventGenerator::MakeInDetector(bool Efficiency)
{
	delete InDetector;
	InDetector = 0;

	double A = TheFlux->GetStartRange();
	double B = TheFlux->GetEndRange();
	double Range = TheFlux->GetRange();
	//double EnStep = (B-A)/(TheFlux->GetBinNumber()-1);
	double EnStep = Range/500.0;

	InDetector = new TH1D("InDetector", "Neutrinos in detector", TheFlux->GetBinNumber(), A, B);

	long double Add;
	for (double EnKin = 0.0; EnKin < Range; EnKin += EnStep)
	//for (double Energy = A; Energy < B; Energy += EnStep)
	{
		//SetEnergy(Energy);
		SetEnergyKin(EnKin);
		//if (EventProbability() <= 0.0) continue;	//speed up?

		Add = FluxIntensity() * EventProbability();
		if (Efficiency)
			Add *= EventEfficiency();

		InDetector->Fill(EnKin+1e-6, Add);
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

double EventGenerator::SampleEnergy(bool Set)	//Sample Energy according to PDF distribution
{
	double Energy = TheFlux->SampleEnergy();
	if (Set)
		SetEnergyKin(Energy);

	return Energy + GetMass();
}

double EventGenerator::SampleInDetector(bool Set)
{
	double Energy = InDetector->GetRandom();
	if (Set)
		SetEnergyKin(Energy);

	return Energy + GetMass();
}

long double EventGenerator::FluxIntensity()	//Get the flux intensity at given energy
{
	return TheFlux->GetIntensity(GetEnergyKin());
}

long double EventGenerator::BackgroundIntensity()	//Get the background intensity at given energy
{
	return TheBox->Background(GetChannel(), GetEnergy());
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

double EventGenerator::GetEnergyKin(int Pow)
{
	return pow(E_Sterile - M_Sterile, Pow);
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

long double EventGenerator::GetSignal()
{
	return lSignal;
}

long double EventGenerator::GetBackground()
{
	return lBackground;
}

long double EventGenerator::GetReducedChi2()
{
	return lRedChi2;
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

void EventGenerator::SetEnergyKin(double X)
{
	E_Sterile = X + M_Sterile;
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

void EventGenerator::SetSignalNumber(long double X)
{
	lSignal = X;
}

void EventGenerator::SetBackgroundNumber(long double X)
{
	lBackground = X;
}

void EventGenerator::SetReducedChi2(long double X)
{
	lRedChi2 = X;
}
