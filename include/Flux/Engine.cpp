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

void Engine::SampleEnergy(std::vector<double> &vE, std::vector<double> &vI)
{
	std::vector<double> vEnergyF, vEnergyR;
	std::vector<double> vIntensF, vIntensR;
	SampleEnergy(vEnergyF, vIntensF, FHC);
	SampleEnergy(vEnergyR, vIntensR, RHC);

	vE.clear();
	vE.insert(vE.end(), vEnergyF.begin(), vEnergyF.end());
	vE.insert(vE.end(), vEnergyR.begin(), vEnergyR.end());

	vI.clear();
	vI.insert(vI.end(), vIntensF.begin(), vIntensF.end());
	vI.insert(vI.end(), vIntensR.begin(), vIntensR.end());
}

void Engine::SampleEnergy(std::vector<double> &vE, std::vector<double> &vI, Current Horn)
{
	vE.clear();
	vI.clear();

	for (unsigned int i = 0; i < vDriver(Horn); ++i)
	{
		vE.push_back(SampleEnergy(Horn, i));
		vI.push_back(IntensitySample(Horn, i, vE.back()));
	}
}

double Engine::SampleEnergy(Current Horn, unsigned int ID)
{
	if (vSampleNu(Horn, ID))
		return vSampleNu(Horn, ID)->GetRandom();
	else
		return -1.0;
}

double Engine::MakeSampler(Detector *Box, std::vector<double> &vInt)
{
	vInt.clear();
	return MakeSampler(Box, vInt, FHC) +
	       MakeSampler(Box, vInt, RHC);
}

double Engine::MakeSampler(Detector *Box, std::vector<double> &vInt, Current Horn)
{
	double Integral = 0;
	double Ret;
	for (unsigned int i = 0; i < vDriver(Horn); ++i)
	{
		Ret = MakeSampler(Box, Horn, i);
		vInt.push_back(Ret);
		Integral += Ret;
	}

	return Integral;
}

double Engine::MakeSampler(Detector *Box)
{
	return MakeSampler(Box, FHC) +
	       MakeSampler(Box, RHC);
}

double Engine::MakeSampler(Detector *Box, Current Horn)
{
	double Integral = 0;
	double Ret;
	for (unsigned int i = 0; i < vDriver(Horn); ++i)
	{
		Ret = MakeSampler(Box, Horn, i);
		Integral += Ret;
	}

	return Integral;
}

double Engine::MakeSampler(Detector *Box, Current Horn, unsigned int ID)
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

		Weight = DecayNumber(Box, Horn, ID);
		hSampleNu->Fill(Energy + EnStep, Weight);

		Integral += EnStep * Weight;
	}

	if (Integral <= 0)
	{
		delete hSampleNu;
		hSampleNu = 0;
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

double Engine::DecayNumber(Detector *Box)
{
	return DecayNumber(Box, FHC) + DecayNumber(Box, RHC);
}

double Engine::DecayNumber(Detector *Box, Current Horn)
{
	double Signal = 0;
	for (unsigned int i = 0; i < vDriver(Horn); ++i)
		Signal += DecayNumber(Box, Horn, i);
	return Signal;
}

double Engine::DecayNumber(Detector *Box, Current Horn, unsigned int ID)
{
	//std::cout << Intensity(Horn, ID) << "\t" << Box->DecayProb(vNeutrino(Horn, ID)) << std::endl;
	return Box->Efficiency(vNeutrino(Horn, ID)) * Intensity(Horn, ID) * Box->DecayProb(vNeutrino(Horn, ID));
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

double Engine::IntensitySample(Current Horn, unsigned int ID)
{
	return IntensitySample(Horn, ID, vNeutrino(Horn, ID)->EnergyKin());
}

double Engine::IntensitySample(Current Horn, unsigned int ID, double Energy)
{
	return vDriver(Horn, ID)->InterpolateIntensity(vSampleNu(Horn, ID), Energy);
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
