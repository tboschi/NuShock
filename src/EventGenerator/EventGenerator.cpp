#include "EventGenerator.h"

EventGenerator::EventGenerator(std::string SMConfig, std::string DetectorConfig, std::string FluxConfig)
//	M_Electron(Const::fMElectron),
//	M_Muon(Const::fMMuon),
//	M_Pion(Const::fMPion),        	
//	M_Kaon(Const::fMKaon)
{
	GenMT = new TRandom3(0);	//19937 Mersenne Twister generator
	M_Sterile_prev = -1.0;
	E_Sterile_prev = -1.0;

	fUserData = -1.0;
	fUserData_prev = -1.0;

	TheBox = new Detector(DetectorConfig, GenMT);
	TheGamma = new Decay();
	TheFlux = new FluxDriver(FluxConfig);	//I have decided this is kinetic energy
	TheProton  = new Nucleon(1, 1);
	TheNeutron = new Nucleon(1, 0);

	Position = new TVector3;
	hSampler = 0;

	SyncUu();	//both parameters are set

	std::string Line, Key;
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

void EventGenerator::SetUserData(double X)
{
	fUserData = X;
}

double EventGenerator::GetUserData()
{
	return fUserData;
}

//MC procedures to detection of event after sampling energy

double EventGenerator::DecayProb()	//reaching the detector and decaying
{						//using sampled energy
	if (GetEnergy() < GetMass())
		return 0.0;
	else
	{
		//double Total = GetUserData() < 0.0 ? TheGamma->Total() : GetUserData();
		double Total = TheGamma->Total();
		double Ratio = TheGamma->Branch(GetChannel()); 
		//std::cout << "Ratio " << Ratio << std::endl;
		double Length = Const::fM2GeV * TheBox->GetElement("Baseline");
		double Lambda = Const::fM2GeV * TheBox->GetZsize();
		double Lorentz = sqrt(GetEnergy(2)/GetMass(2) - 1.0);	//betagamma, to invert
		return exp(- Total * Length / Lorentz) * (1 - exp(- Total * Lambda / Lorentz)) * Ratio;
	}
}

double EventGenerator::ScatterProb(double Eh)	//decaying inside the detector
{							//using heavy neutrino energy
	if (Eh < GetMass())
		return 0.0;
	else
	{
		double Total = TheGamma->Total(GetEnhancement());
		//std::cout << "T " << Total << std::endl;
		//double Ratio = TheGamma->Branch(GetChannel());
		double Lambda = Const::fM2GeV * TheBox->GetZsize();
		double Lorentz = sqrt(Eh*Eh/GetMass(2) - 1.0);	//betagamma, to invert
		return (1 - exp(- Total * Lambda / Lorentz));
	}
}

void EventGenerator::SetEnhancement(double X)
{
	fEnhance = X;
}

double EventGenerator::GetEnhancement()
{
	return fEnhance;
}

//double EventGenerator::ScatterProb(double Eh)	//decaying inside the detector
//{							//using heavy neutrino energy
//	if (Eh < GetMass())
//		return 0.0;
//	else
//	{
//		double Ratio = TheGamma->Branch(GetChannel()); 
//		double Lambda = Const::fM2GeV * TheBox->GetZsize();
//		double Lorentz = sqrt(Eh*Eh/GetMass(2) - 1.0);	//betagamma, to invert
//		return (1 - exp(- Total * Lambda / Lorentz)) * Ratio;
//	}
//}


double EventGenerator::EventEfficiency()
{
	////std::cout << "call" << std::endl;
	return TheBox->Efficiency(GetEnergy(), GetMass());
}

//modifications:
//hMod, hProb
//even in .h file
double EventGenerator::DecayNumber(double EnergyKin, bool Efficiency)
{
	SetEnergyKin(EnergyKin);

	//std::cout << " I " << FluxIntensity(1) << "\t P " << DecayProb() << std::endl;
	double Signal = FluxIntensity(1) * DecayProb();
	if (Efficiency)
		Signal *= EventEfficiency();

	return Signal;
}

double EventGenerator::ScatterNumber(double Energy, bool Efficiency)	//watchout for units!!! not set
{
	SetEnergy(Energy);
	TLorentzVector Nu(0, 0, Energy, Energy);

	TheProton->SetProbe(Nu);
	TheNeutron->SetProbe(Nu);

	if (IsChanged())
		fTotalXSec = Kine::BooleIntegration(this);	//booleint will integrate the Variable func

	double Signal = (GetUe()*GetUe() + GetUm()*GetUm() + GetUt()*GetUt()) * fTotalXSec;

	if (Efficiency)
		Signal *= EventEfficiency();

	return Signal;
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

bool EventGenerator::DecayInDetector()		//Second step: is the decay inside the detector?
{
	return (GenMT->Rndm() <= DecayProb());
}

bool EventGenerator::EventDetectable()	//Third step: is the detector able to detect it?
{
	return (GenMT->Rndm() <= EventEfficiency());
}

//Kinematics
bool EventGenerator::IsAllowed()
{
	return TheGamma->IsAllowed(GetChannel());
}

int EventGenerator::EventKinematics()	//Fourth step: simulate the phase space of the decay
{
	TLorentzVector N_vec(0, 0, GetMomentum(), GetEnergy());	//Lab frame for the heavy neutrino

	double Radius = sqrt(pow(TheBox->GetElement("Width", 1),2) + pow(TheBox->GetElement("Height", 1),2));
	double SigmaT = atan2(Radius, TheBox->GetElement("Baseline"));	//3sigma will be inside the detector (better distribution needed)
	N_vec.SetTheta(abs(GenMT->Gaus(0, SigmaT)));	//abs needed to avoid degeneracy with phi (?)
	N_vec.SetPhi(GenMT->Uniform(-Const::fPi, Const::fPi));

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
	TLorentzVector Vi = TheGamma->GetDecayProduct(i, pdg);
	//Particle * Pi = TheGamma->GetDecayProduct(pdg, i);

	TVector3 Coord = *Position;
	Particle *P = new Particle(pdg, Vi, Coord);	//mainly charged particles

	if (Smear)
	{
		if (P->TrackIn() < 0)
			TheBox->TrackLength(P);
		TheBox->SignalSmearing(P);
	}

	if (pdg != 12 && pdg != 111 && !TheBox->IsDetectable(P))
	{
		delete P;
		P = 0;
	}

	return P;
}

void EventGenerator::Pi0Decay(Particle *Pi0, Particle *&PA, Particle *&PB, bool Smear)
{
	//in rest frame
	double M_Pion0 = Const::fMPion0;
	TLorentzVector GammaA(0, 0, M_Pion0/2.0, M_Pion0/2.0); 
	TLorentzVector GammaB(0, 0, -M_Pion0/2.0, M_Pion0/2.0); 

	TVector3 vBoost(Pi0->GetP4().BoostVector());
	TVector3 Start(Pi0->Position());		//starting point is vertex
	double Theta = GenMT->Uniform(-Const::fPi, Const::fPi);
	double Phi = GenMT->Uniform(-Const::fPi, Const::fPi);

	GammaA.SetTheta(Theta);
	GammaB.SetTheta(Theta + Const::fPi);
	GammaA.SetPhi(Phi);
	GammaB.SetPhi(Phi + Const::fPi);

	GammaA.Boost(vBoost);
	GammaB.Boost(vBoost);

	TVector3 MoveA(GammaA.Vect().Unit());
	TVector3 MoveB(GammaB.Vect().Unit());
	MoveA *= TheBox->GammaDecay();
	MoveB *= TheBox->GammaDecay();
	MoveA += Start;
	MoveB += Start;

	delete PA, PB;
	PA = 0, PB = 0;

	PA = new Particle(22, GammaA, MoveA);	//here are the photons
	PB = new Particle(22, GammaB, MoveB);	//position should be different

	if (Smear)
	{
		if (PA->TrackIn() < 0)
			TheBox->TrackLength(PA);
		if (PB->TrackIn() < 0)
			TheBox->TrackLength(PB);
		TheBox->SignalSmearing(PA);
		TheBox->SignalSmearing(PB);
	}

	if (!TheBox->IsDetectable(PA))
	{
		delete PA;
		PA = 0;
	}
	if (!TheBox->IsDetectable(PB))
	{
		delete PB;
		PB = 0;
	}
}
	
void EventGenerator::GeneratePosition()
{
	double PosX = GenMT->Uniform(TheBox->GetXsize());
	double PosY = GenMT->Uniform(TheBox->GetYsize());
	double PosZ = GenMT->Uniform(TheBox->GetZsize());
	Position->SetX(PosX);
	Position->SetY(PosY);
	Position->SetZ(PosZ);
}

//Flux as PDF for MC
//modifications
//baseline
//pot
//area
void EventGenerator::MakeFlux(bool Mass, bool TotalPOT)	//Generate the flux for heavy neutrinos
{
	bool Ret;
	if (Mass)	//heavy neutrino flux
		Ret = TheFlux->MakeFlux(GetMass());
	else		//light neutrino
		Ret = TheFlux->MakeFlux(0);

	if (Ret)
	{
		TheFlux->SetBaseline(TheBox->GetElement("Baseline"));
		double Y = TotalPOT ? 1.0e7 * TheBox->GetElement("Years") : 1.0;
		TheFlux->SetPOT(Y * TheBox->GetElement("POT/s"));
		TheFlux->SetArea(TheBox->GetElement("Height", 1)*TheBox->GetElement("Width", 1)*1.0e4);
	}
}

void EventGenerator::MakeSampler()
{
	delete hSampler;
	hSampler = 0;

	double Start, End;
	double EnStep = GetRange(Start, End)/GetBinNumber();

	hSampler = new TH1D("sampler", "Neutrinos in detector", GetBinNumber(), Start, End);

	double Sgn;
	for (double EnKin = Start; EnKin < End; EnKin += EnStep)
	{
		SetEnergyKin(EnKin);

		Sgn = EnStep * FluxIntensity(1) * DecayProb();
		hSampler->Fill(GetEnergy()+1e-6, Sgn);
	}
}

double EventGenerator::SampleEnergy(bool Set)	//Sample Energy according to PDF distribution
{
	double Energy = hSampler->GetRandom();
	if (Set)
		SetEnergy(Energy);

	return Energy;
}

double EventGenerator::NeutIntensity(bool Mass)	//Get the neutrino flux intensity at given energy
{
	if (Mass)
		return TheFlux->GetIntensity(GetEnergyKin(), 1, 1);
	else 
		return TheFlux->GetIntensity(GetEnergy(), 1, 0);
}

double EventGenerator::AntiIntensity(bool Mass)	//Get the antineutrino flux intensity at given energy
{
	if (Mass)
		return TheFlux->GetIntensity(GetEnergyKin(), 0, 1);
	else
		return TheFlux->GetIntensity(GetEnergy(), 0, 0);
}

double EventGenerator::FluxIntensity(bool Mass)	//Get the flux intensity at given energy
{
	return NeutIntensity(Mass) + AntiIntensity(Mass);
}

double EventGenerator::ScaleXSec(double IntP, double IntN)	//including probability (one single integration)
{
	double MolarMass = 1.0;
	int A = 1.0, Z = 1.0;
	
	if (TheBox->GetElement("InTarget") == 1 || TheBox->GetElement("InTarget") == 2)	//Argon
	{
		Z = 18;
		A = 40;
	}
	else if (TheBox->GetElement("InTarget") == 3)	//Iron
	{
		Z = 26;
		A = 56;
	}

	MolarMass = double(A);	//grams
	double Area = TheBox->GetElement("Width", 1) * TheBox->GetElement("Height", 1);
	double RhoTarget = Const::fNa * TheBox->GetElement("Weight", 1) * 1e6 / MolarMass / Area;	//surface density of targets

	return RhoTarget * Const::fGeV2cm * (Z * IntP + (A-Z) * IntN);	//xs in cm2
}

double EventGenerator::Variable(double dt)
{
	double A_P, A_N;
	double B_P, B_N;

	double dQ2_P =  TheProton->Q2Lim(A_P, B_P);
	double dQ2_N = TheNeutron->Q2Lim(A_N, B_N);

	if (dQ2_P > 0 && dQ2_N > 0)
	{
		//set Q2
		TheProton->SetQ2( dQ2_P*dt + A_P);
		TheNeutron->SetQ2(dQ2_N*dt + A_N);

		//compute integrand
		double Int_P = 0.0, Int_N = 0.0;
		if (NeutIntensity(0) > 0.0)
		{
			Int_P += NeutIntensity(0) *  TheProton->dSigmadQ2(1);
			Int_N += NeutIntensity(0) * TheNeutron->dSigmadQ2(1);
		}
		if (AntiIntensity(0) > 0.0)
		{
			Int_P += AntiIntensity(0) *  TheProton->dSigmadQ2(-1);
			Int_N += AntiIntensity(0) * TheNeutron->dSigmadQ2(-1);
		}

		Int_P *= dQ2_P * ScatterProb(TheProton->GetHeavyE());
		Int_N *= dQ2_N * ScatterProb(TheNeutron->GetHeavyE());

		//return TotalXSec(Int_Pn, Int_Nn, Int_Pa, Int_Na);
		return ScaleXSec(Int_P, Int_N);
	}
	else 
		return 0.0;
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

double EventGenerator::GetUe(int Pow)	//it is the last value
{
	return pow(Ue, Pow);
}

double EventGenerator::GetUm(int Pow)
{
	return pow(Um, Pow);
}

double EventGenerator::GetUt(int Pow)
{
	return pow(Ut, Pow);
}

double EventGenerator::GetRange(double &Start, double &End)
{
	Start = TheFlux->GetRangeStart();
	End = TheFlux->GetRangeEnd();
	return End-Start;
}

int EventGenerator::GetBinNumber()
{
	return TheFlux->GetBinNumber();
}


//Set functions
void EventGenerator::SetChannel(std::string Ch, bool Eff, char Couple)
{
	if (Ch == "R")
		sChannel.assign(RandomChannel());
	else sChannel.assign(Ch); 

	if (Eff)
		TheBox->SetEfficiency(Ch, Couple);
}

void EventGenerator::SetMass(double X)
{
	M_Sterile = X;
	TheGamma->SetMass(X);

	TheProton->SetPS(X);
	TheNeutron->SetPS(X);
}

void EventGenerator::SetEnergy(double X)
{
	E_Sterile = X;
}

void EventGenerator::SetEnergyKin(double X)
{
	E_Sterile = X + M_Sterile;
}

void EventGenerator::SetUe(double X, bool GvF)	//GvF == 1 -> Gamma, GvF == 0 -> Flux
{
	Ue = X;

	if (Sync || GvF)
		TheGamma->SetUe(X);
	if (Sync || !GvF)
		TheFlux->SetUe(X);
}

void EventGenerator::SetUm(double X, bool GvF) 	//GvF == 1 -> Gamma, GvF == 0 -> Flux
{
	Um = X;

	if (Sync || GvF)
		TheGamma->SetUm(X);
	if (Sync || !GvF)
		TheFlux->SetUm(X);
}

void EventGenerator::SetUt(double X, bool GvF)	//GvF == 1 -> Gamma, GvF == 0 -> Flux
{
	Ut = X;

	if (Sync || GvF)
		TheGamma->SetUt(X);
	if (Sync || !GvF)
		TheFlux->SetUt(X);
}

void EventGenerator::SyncUu(bool B)	//B == 1 -> Gamma&Flux same U, B == 0 -> Indipendent (Gamma default)
{
	Sync = B;
}

bool EventGenerator::IsChanged()
{
	bool Ret = (fabs(M_Sterile - M_Sterile_prev) > 1e-9 || 
		    fabs(E_Sterile - E_Sterile_prev) > 1e-9 || 
		    fabs(fUserData - fUserData_prev) > 1e-30 );

	M_Sterile_prev = M_Sterile;
	E_Sterile_prev = E_Sterile;
	fUserData_prev = fUserData;

	return Ret;
}
