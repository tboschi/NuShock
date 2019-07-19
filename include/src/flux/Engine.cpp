#include "Engine.h"
#include <omp.h>

Engine::Engine(std::string fc) :
	fluxConfig(fc)
{
}

Engine::~Engine()
{
	for (iD = mDriver.begin(); iD != mDriver.end(); ++iD)
		delete iD->second;

	for (iS = sampleNu.begin(); iS != sampleNu.end(); ++iS)
		delete iS->second;
}

void Engine::Reset()
{
	mNeutrino.clear();

	for (iD = mDriver.begin(); iD != mDriver.end(); ++iD)
		delete iD->second;
	mDriver.clear();

	for (iS = sampleNu.begin(); iS != sampleNu.end(); ++iS)
		delete iS->second;
	sampleNu.clear();
}

//load neutrino to mDriver with specific current. ID is user input and must be remembered
void Engine::BindNeutrino(std::string uuid, Neutrino &N, Current horn)
{
	uuid += "_" + HornName(horn);

	mNeutrino[uuid] = &N;

	switch (horn)
	{
		case FHC:
			mDriver[uuid] = new Driver(fluxConfig, FHC);
			break;
		case RHC:
			mDriver[uuid] = new Driver(fluxConfig, RHC);
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

void Engine::SampleEnergy(std::vector<double> &vE, std::vector<double> &vI, Current horn)
{
	vE.clear();
	vI.clear();

	double integral = 0;
	for (iD = mDriver.begin(); iD != mDriver.end(); ++iD)
	{
		if (iD->first.find(HornName(horn)) != std::string::npos)
		{
			vE.push_back(SampleEnergy(iD->first));
			vI.push_back(IntensitySample(iD->first, vE.back()));
		}
	}
}

double Engine::SampleEnergy(std::string uuid)
{
	uuid = "sample_" + uuid;
	if (sampleNu[uuid])
		return sampleNu[uuid]->GetRandom();
	else
		return -1.0;
}

//create sample for every mDriver
double Engine::MakeSampler(Detector *box, double ue, double um, double ut)
{
	double integral = 0;
	for (iD = mDriver.begin(); iD != mDriver.end(); ++iD)
		integral += MakeSampler(box, iD->first, ue, um, ut);

	return integral;
}

//create sample for every mDriver of same current
double Engine::MakeSampler(Detector *box, Current horn, double ue, double um, double ut)
{
	if (horn == both)
		return MakeSampler(box, ue, um, ut);

	double integral = 0;
	for (iD = mDriver.begin(); iD != mDriver.end(); ++iD)
		if (iD->first.find(HornName(horn)) != std::string::npos)
			integral += MakeSampler(box, iD->first, ue, um, ut);

	return integral;
}

//create sample for every mDriver and saving each integral into a vector
double Engine::MakeSampler(Detector *box, std::vector<double> &vInt, double ue, double um, double ut)
{
	vInt.clear();
	double integral = 0;
	for (iD = mDriver.begin(); iD != mDriver.end(); ++iD)
	{
		double ret = MakeSampler(box, iD->first, ue, um, ut);
		vInt.push_back(ret);
		integral += ret;
	}

	return integral;
}

//create sample for every mDriver with same current and saving each integral into a vector
double Engine::MakeSampler(Detector *box, std::vector<double> &vInt, Current horn, double ue, double um, double ut)
{
	if (horn == both)
		return MakeSampler(box, vInt, ue, um, ut);

	vInt.clear();
	double integral = 0;
	for (iD = mDriver.begin(); iD != mDriver.end(); ++iD)
	{
		if (iD->first.find(HornName(horn)) != std::string::npos)
		{
			double ret = MakeSampler(box, iD->first, ue, um, ut);
			vInt.push_back(ret);
			integral += ret;
		}
	}

	return integral;
}

//// user should never use this function
double Engine::MakeSampler(Detector *box, std::string uuid, double ue, double um, double ut)
{
	//delete old sample
	std::string name = "sample_" + uuid;

	delete sampleNu[name];

	double start, end;
	double enStep = RangeWidth(start, end);
	TH1D *hSample = new TH1D(name.c_str(), "Neutrinos in detector", BinNumber(), start, end);

	if (ue < 0)
		ue = mNeutrino[uuid]->Ue();
	if (um < 0)
		um = mNeutrino[uuid]->Um();
	if (ut < 0)
		ut = mNeutrino[uuid]->Ut();

	mNeutrino[uuid]->SetMixings(ue, um, ut);

	double integral = 0;
	for (double energy = start; energy < end; energy += enStep)
	{
		mNeutrino[uuid]->SetEnergy(energy);

		double weight = DecayNumber(box, uuid);
		hSample->Fill(energy + enStep, weight);

		integral += enStep * weight;
	}

	if (integral <= 0)
	{
		delete hSample;
		hSample = 0;
	}

	sampleNu[name] = hSample;

	return integral;
}

//make all the fluxes
double Engine::DecayNumber(Detector *box)
{
	double events = 0;
	for (iD = mDriver.begin(); iD != mDriver.end(); ++iD)
		events += DecayNumber(box, iD->first);

	return events;
}

double Engine::DecayNumber(Detector *box, Current horn)
{
	double events = 0;
	for (iD = mDriver.begin(); iD != mDriver.end(); ++iD)
		if (iD->first.find(HornName(horn)) != std::string::npos)
			events += DecayNumber(box, iD->first);

	return events;
}

double Engine::DecayNumber(Detector *box, std::string uuid)
{
	return box->Efficiency(*mNeutrino[uuid]) * box->DecayProb(*mNeutrino[uuid]) * Intensity(uuid);
}

//make all the fluxes
void Engine::MakeFlux()
{
	for (iD = mDriver.begin(); iD != mDriver.end(); ++iD)
		MakeFlux(iD->first);
}

void Engine::MakeFlux(Current horn)
{
	for (iD = mDriver.begin(); iD != mDriver.end(); ++iD)
		if (iD->first.find(HornName(horn)) != std::string::npos)
			MakeFlux(iD->first);
}

void Engine::MakeFlux(std::string uuid)
{
	mDriver[uuid]->MakeFlux(*mNeutrino[uuid]);
}

double Engine::Intensity(std::string uuid)
{
	return mDriver[uuid]->Intensity(*mNeutrino[uuid]);
}

double Engine::IntensitySample(std::string uuid)
{
	return IntensitySample(uuid, mNeutrino[uuid]->EnergyKin());
}

double Engine::IntensitySample(std::string uuid, double Energy)
{
	return mDriver[uuid]->InterpolateIntensity(sampleNu[uuid], Energy);
}

void Engine::ScaleToDetector(Detector *box)
{
	ScaleBaseline(box);
	ScalePOT(box);
	ScaleArea(box);
}

void Engine::ScaleBaseline(Detector *box)
{
	ScaleBaseline(box->Get("Baseline"));
}

void Engine::ScalePOT(Detector *box)
{
	ScalePOT(1.0e7 * box->Get("Years") * box->Get("POT/s"));
}

void Engine::ScaleArea(Detector *box)
{
	ScaleArea(box->Area() * 1.0e4);
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
	for (iD = mDriver.begin(); iD != mDriver.end(); ++iD)
		Scale(X, iD->first);
}

void Engine::Scale(double X, Current horn)
{
	for (iD = mDriver.begin(); iD != mDriver.end(); ++iD)
		if (iD->first.find(HornName(horn)) != std::string::npos)
			Scale(X, iD->first);
}

void Engine::Scale(double X, std::string uuid)
{
	mDriver[uuid]->Scale(X);
}

double Engine::BinNumber()
{
	if (mDriver.begin()->second)
		return mDriver.begin()->second->BinNumber();
	else
		return -1.0;
}

double Engine::RangeWidth(double &start, double &end)
{
	if (mDriver.begin()->second)
		return mDriver.begin()->second->RangeWidth(start, end);
	else
	{
		start = -1;
		end = -1;
	}
}


std::string Engine::HornName(Current horn)
{
	std::string name;
	switch (horn)
	{
		case FHC:
			name = "FHC";
			break;
		case RHC:
			name = "RHC";
			break;
	}

	return name;
}

