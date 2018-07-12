#include "Flux/Engine.h"

Engine::Engine(std::string FluxConfig, unsigned int nFHC, unsigned int nRHC)
{
	for (unsigned int i = 0; i < nFHC; ++i)
		vDriverFHC.push_back(new Driver(FluxConfig, FHC));

	for (unsigned int i = 0; i < nRHC; ++i)
		vDriverRHC.push_back(new Driver(FluxConfig, RHC));

	vNeutrinoFHC.resize(vDriver(FHC));
	vNeutrinoRHC.resize(vDriver(RHC));

	vSampleNuFHC.resize(vDriver(FHC));
	vSampleNuRHC.resize(vDriver(RHC));
}

Engine::~Engine()
{
	for (unsigned int i = 0; i < vDriver(FHC); ++i)
		delete vDriver(FHC, i);

	for (unsigned int i = 0; i < vDriver(RHC); ++i)
		delete vDriver(RHC, i);
}

//load neutrino to driver with specific current. ID is user input and must be remembered
void Engine::BindNeutrino(Neutrino *N, Current Horn, unsigned int ID)
{
	delete vNeutrino(Horn, ID);

	switch (Horn)
	{
		case FHC:
			vNeutrinoFHC.at(ID) = N;
			break;
		case RHC:
			vNeutrinoRHC.at(ID) = N;
			break;
	}
}

std::vector<double> Engine::SampleEnergy(bool Set)
{
	std::vector<double> vEnergyF, vEnergyR;
	vEnergyF = SampleEnergy(FHC, Set);
	vEnergyR = SampleEnergy(RHC, Set);

	vEnergyF.insert(vEnergyF.end(), vEnergyR.begin(), vEnergyR.end());
	return vEnergyF;
}

std::vector<double> Engine::SampleEnergy(Current Horn, bool Set)
{
	std::vector<double> vEnergy;
	for (unsigned int i = 0; i < vDriver(Horn); ++i)
		vEnergy.push_back(SampleEnergy(Horn, i, Set));

	return vEnergy;
}

double Engine::SampleEnergy(Current Horn, unsigned int ID, bool Set)
{
	double Energy = -1.0;
	if (vSampleNu(Horn, ID))
	{
		Energy = vSampleNu(Horn, ID)->GetRandom();
		if (Set)
			vNeutrino(Horn, ID)->SetEnergy(Energy);
		return Energy;
	}
	else
		return -1.0;
}

double Engine::MakeSampler(Detector *Box, bool Eff)
{
	return MakeSampler(Box, Eff, FHC) +
	       MakeSampler(Box, Eff, RHC);
}

double Engine::MakeSampler(Detector *Box, bool Eff, Current Horn)
{
	double Integral = 0;
	for (unsigned int i = 0; i < vDriver(Horn); ++i)
		Integral += MakeSampler(Box, Eff, Horn, i);

	return Integral;
}

double Engine::MakeSampler(Detector *Box, bool Eff, Current Horn, unsigned int ID)
{
	delete vSampleNu(Horn, ID);

	std::stringstream ssName;
	ssName << "sample_";
	switch (Horn)
	{
		case FHC:
			ssName << "FHC";
			break;
		case RHC:
			ssName << "RHC";
			break;
	}
	ssName << "_" << ID;

	double Start, End;
	double EnStep = RangeWidth(Start, End);
	TH1D *hSampleNu = new TH1D(ssName.str().c_str(), "Neutrinos in detector", BinNumber(), Start, End);

	double Weight, Integral = 0;
	for (double Energy = Start; Energy < End; Energy += EnStep)
	{
		vNeutrino(Horn, ID)->SetEnergy(Energy);

		Weight = DecayNumber(Box, Eff, Horn, ID);
		hSampleNu->Fill(Energy + EnStep, Weight);

		Integral += EnStep * Weight;
	}

	switch (Horn)
	{
		case FHC:
			vSampleNuFHC.at(ID) = hSampleNu;
			break;
		case RHC:
			vSampleNuRHC.at(ID) = hSampleNu;
			break;
	}

	return Integral;
}

double Engine::DecayNumberIntegrated(Detector *Box, bool Eff)
{
	return MakeSampler(Box, Eff);	//this is needed
;
	//return DecayNumberIntegrated(Box, FHC) + 
	       //DecayNumberIntegrated(Box, RHC);
}

double Engine::DecayNumberIntegrated(Detector *Box, Current Horn)
{
	double Events = 0;
	for (unsigned int i = 0; i < vDriver(Horn); ++i)
	       Events += vSampleNu(Horn, i)->Integral("WIDTH");

	return Events;
}

double Engine::DecayNumber(Detector *Box, bool Eff)
{
	return DecayNumber(Box, Eff, FHC) + DecayNumber(Box, Eff, RHC);
}

double Engine::DecayNumber(Detector *Box, bool Eff, Current Horn)
{
	double Signal = 0;
	for (unsigned int i = 0; i < vDriver(Horn); ++i)
		Signal += DecayNumber(Box, Eff, Horn, i);
	return Signal;
}

double Engine::DecayNumber(Detector *Box, bool Eff, Current Horn, unsigned int ID)
{
	//std::cout << "n " << Intensity(Horn, ID) << "\t" << Box->DecayProb(vNeutrino(Horn, ID)) << std::endl;
	return (Eff ? Box->Efficiency(vNeutrino(Horn, ID)) : 1.0) * Intensity(Horn, ID) * Box->DecayProb(vNeutrino(Horn, ID));
}

void Engine::MakeFlux()
{
	MakeFlux(FHC);
	MakeFlux(RHC);
}

void Engine::MakeFlux(Current Horn)
{
	for (unsigned int i = 0; i < vDriver(Horn); ++i)
		MakeFlux(Horn, i);
}

void Engine::MakeFlux(Current Horn, unsigned int ID)
{
	vDriver(Horn, ID)->MakeFlux(vNeutrino(Horn, ID));
}

double Engine::Intensity(Current Horn, unsigned int ID)
{
	return vDriver(Horn, ID)->Intensity(vNeutrino(Horn, ID));
}

void Engine::ScaleDetector(Detector *Box)
{
	ScaleBaseline(Box);
	ScalePOT(Box);
	ScaleArea(Box);
}

void Engine::ScaleBaseline(Detector *Box)
{
	ScaleBaseline(Box->Get("Baseline"));
}

void Engine::ScalePOT(Detector *Box)
{
	ScalePOT(1.0e7 * Box->Get("Years") * Box->Get("POT/s"));
}

void Engine::ScaleArea(Detector *Box)
{
	ScaleArea(Box->Area() * 1.0e4);
}

void Engine::ScaleBaseline(double Baseline)
{
	Scale(1.0/(Baseline*Baseline));
}

void Engine::ScalePOT(double POT)
{
	Scale(POT);
}

void Engine::ScaleArea(double Area)
{
	Scale(Area);
}

void Engine::Scale(double X)
{
	Scale(X, FHC);
	Scale(X, RHC);
}

void Engine::Scale(double X, Current Horn)
{
	for (unsigned int i = 0; i < vDriver(Horn); ++i)
		Scale(X, Horn, i);
}

void Engine::Scale(double X, Current Horn, unsigned int ID)
{
	vDriver(Horn, ID)->Scale(X);
}

double Engine::BinNumber()
{
	if (vDriver(FHC))
		return vDriver(FHC, 0)->BinNumber();
	else if (vDriver(RHC))
		return vDriver(RHC, 0)->BinNumber();
}

double Engine::RangeWidth(double &Start, double &End)
{
	if (vDriver(FHC))
		return vDriver(FHC, 0)->RangeWidth(Start, End);
	else if (vDriver(RHC))
		return vDriver(RHC, 0)->RangeWidth(Start, End);
}

//access to pointers
//
unsigned int Engine::vNeutrino(Current Horn)
{
	switch (Horn)
	{
		case FHC:
			return vNeutrinoFHC.size();
		case RHC:
			return vNeutrinoRHC.size();
	}
}

Neutrino* Engine::vNeutrino(Current Horn, unsigned int i)
{
	if (i < vNeutrino(Horn))
	{
		switch (Horn)
		{
			case FHC:
				return vNeutrinoFHC.at(i);
			case RHC:
				return vNeutrinoRHC.at(i);
		}
	}
	else
		return NULL;
}

TH1D* Engine::vSampleNu(Current Horn, unsigned int i)
{
	if (i < vNeutrino(Horn))
	{
		switch (Horn)
		{
			case FHC:
				return vSampleNuFHC.at(i);
			case RHC:
				return vSampleNuRHC.at(i);
		}
	}
	else
		return NULL;
}

unsigned int Engine::vDriver(Current Horn)
{
	switch (Horn)
	{
		case FHC:
			return vDriverFHC.size();
		case RHC:
			return vDriverRHC.size();
	}
}

Driver* Engine::vDriver(Current Horn, unsigned int i)
{
	if (i < vDriver(Horn))
	{
		switch (Horn)
		{
			case FHC:
				return vDriverFHC.at(i);
			case RHC:
				return vDriverRHC.at(i);
		}
	}
	else
		return NULL;
}
