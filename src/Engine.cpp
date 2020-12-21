#include "Engine.h"

Engine::Engine(const std::string &card) : _drive(card)
{
}

// std::unordered_map<Nu::Flavour, std::shared_ptr<TH1D> > _distributions
void Engine::MakeDistribution(const Neutrino &N, const Detector &box) {
	_box = box;
	_N = N;
	_distributions = _drive.MakeFlux(_N);
	for (auto &fc : _distributions)
		fc.second->Scale(_box.Exposure());
}

// store channel for neutrino decay
void Engine::SetChannels(const std::vector<Channel::Name> &chans) {
	_channels = chans;
}

// store channel for neutrino decay
void Engine::SetChannels(Channel::Name chans) {
	_channels.push_back(chan);
}


// use in FullSensitivity like this
/*
 * auto Numbers = [&](double lu2) {
 * 			double ue = kUe ? std::pow(10.0, 0.5 * lu2) : 0.;
 *			double um = kUm ? std::pow(10.0, 0.5 * lu2) : 0.;
 *			double ut = kUt ? std::pow(10.0, 0.5 * lu2) : 0.;
 * 			return Engine.NumberEvents(ue, um, ut);
 * 			}
 * std::vector<double> BinarySearch::solve(Numbers, lu2bot, lu2top);
 */

double Engine::operator()(double x) {
	//std::cout << "computing number of events at " << lu2 << " :\t";

	double tot = 0;
	for (int i = 0; i < channels.size(); ++i)
	{
		engine->SetDecay(channels[i]);
		tot += engine->MakeSampler(box, ue, um, ut);
	}
	//std::cout << tot << " vs " << threshold << std::endl;
	return tot - (t < 0 ? threshold : t);
}

//make all the fluxes
double Engine::EventNumber(double ue, double um, double ut)
{
	// compute decay rates
	DecayRates hnl(_N);

	double T = hnl.Total(ue, um, ut);	// total decay rate

	double tot = 0.;
	for (Channel::Name chan : _channels) {
		double num = 0.;
		for (const auto &fc : _distributions) {
			double dist = 0.;
			for (int bin = 1; bin <= fc.second->GetNbinsX(); ++bin) {
				_N.SetE(fc.second->GetBinCenter(bin));
				if (_N.EKin() > 0. && _N.Beta() < 1.)
					dist += fc.second->GetBinContent(bin)
				      	* fc.second->GetBinWidth(bin)
				      	* _box.Efficiency(chan, _N.M(), _N.E())
				      	* _box.Probability(T / (_N.Beta() * _N.Gamma()));
			}
		
			switch (fc.first) {
				case Nu::E_:
				case Nu::Eb:
					num += dist * ue*ue;
					break;
				case Nu::M_:
				case Nu::Mb:
					num += dist * um*um;
					break;
				case Nu::T_:
				case Nu::Tb:
					num += dist * ut*ut;
					break;
				default:
					break;
			}
		}
		tot += num * hnl.Branch(chan, ue, um, ut);
	}

	return tot;
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

Neutrino& Engine::GetNeutrino(std::string uuid)
{
	if (mNeutrino.count(uuid))
		return mNeutrino[uuid];
	else
		Neutrino();
}

//load neutrino to mDriver with specific current. ID is user input and must be remembered
void Engine::BindNeutrino(std::string uuid, Neutrino &N, Current horn)
{
	//uuid += "_" + HornName(horn);

	if (mNeutrino.count(uuid))
	{
		std::cout << "UUID already in use! Call ReleaseNeutrino(UUID) first" << std::endl;
		return;
	}

	mNeutrino[uuid] = N;

	bool type;
	switch (horn)
	{
		case RHC:
			type = false;
			//std::cout << "RHCtype " << type << std::endl;
			break;
		case FHC:
		default:
			//std::cout << "FHCtype " << type << std::endl;
			type = true;
			break;
	}

	//std::cout << "type " << type << std::endl;

	mDriver[uuid] = new Driver(fluxConfig, type);
}

void Engine::ReleaseNeutrino(std::string uuid)
{
	mNeutrino.erase(uuid);
}

void Engine::SampleEnergy(std::map<std::string, double> &mE,
			  std::map<std::string, double> &mI)
{
	mE.clear();
	mI.clear();
	for (iD = mDriver.begin(); iD != mDriver.end(); ++iD)
	{
		mE[iD->first] = SampleEnergy(iD->first);
		mI[iD->first] = IntensitySample(iD->first, mE[iD->first]);
	}
}

void Engine::SampleEnergy(std::map<std::string, double> &mE,
			  std::map<std::string, double> &mI,
			  Current horn)
{
	if (horn == both)
		return SampleEnergy(mE, mI);

	mE.clear();
	mI.clear();
	for (iD = mDriver.begin(); iD != mDriver.end(); ++iD)
	{
		if (iD->first.find(HornName(horn)) != std::string::npos)
		{
			mE[iD->first] = SampleEnergy(iD->first);
			mI[iD->first] = IntensitySample(iD->first, mE[iD->first]);
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
double Engine::MakeSampler(Detector *box,
			   double ue, double um, double ut)
{
	double integral = 0;
	for (iD = mDriver.begin(); iD != mDriver.end(); ++iD)
		integral += MakeSampler(box, iD->first, ue, um, ut);

	return integral;
}

//create sample for every mDriver of same current
double Engine::MakeSampler(Detector *box, Current horn,
			   double ue, double um, double ut)
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
double Engine::MakeSampler(Detector *box, std::map<std::string, double> &mInt,
			   double ue, double um, double ut)
{
	mInt.clear();
	double integral = 0;
	for (iD = mDriver.begin(); iD != mDriver.end(); ++iD)
	{
		mInt[iD->first] = MakeSampler(box, iD->first, ue, um, ut);
		integral += mInt[iD->first];
		//std::cout << "making sample " << iD->first << " -> " << mInt[iD->first] << std::endl;
	}

	return integral;
}

//create sample for every mDriver with same current and saving each integral into a vector
double Engine::MakeSampler(Detector *box, std::map<std::string, double> &mInt, Current horn,
			   double ue, double um, double ut)
{
	if (horn == both)
		return MakeSampler(box, mInt, ue, um, ut);

	mInt.clear();
	double integral = 0;
	for (iD = mDriver.begin(); iD != mDriver.end(); ++iD)
	{
		if (iD->first.find(HornName(horn)) != std::string::npos)
		{
			mInt[iD->first] = MakeSampler(box, iD->first, ue, um, ut);
			integral += mInt[iD->first];
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
		ue = mNeutrino[uuid].Ue();
	if (um < 0)
		um = mNeutrino[uuid].Um();
	if (ut < 0)
		ut = mNeutrino[uuid].Ut();

	//std::cout << "setting mixings to " << ue << ", " << um << ", " << ut << std::endl;
	mNeutrino[uuid].SetMixings(ue, um, ut);

	double integral = 0;
	for (double energy = start; energy < end; energy += enStep)
	{
		mNeutrino[uuid].SetEnergy(energy + enStep/2.0);

		double weight = DecayNumber(box, uuid);
		hSample->Fill(energy + enStep/2.0, weight);

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
double Engine::DecayNumber(Detector &box)
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





void Engine::MakeFlux(Current horn)
{
	for (iD = mDriver.begin(); iD != mDriver.end(); ++iD)
		if (iD->first.find(HornName(horn)) != std::string::npos)
			MakeFlux(iD->first);
}

void Engine::MakeFlux(std::string uuid)
{
	//std::cout << "Making flux " << uuid << std::endl;
	mDriver[uuid]->MakeFlux(mNeutrino[uuid]);
}

double Engine::Intensity(std::string uuid)
{
	return mDriver[uuid]->Intensity(mNeutrino[uuid]);
}

double Engine::IntensitySample(std::string uuid)
{
	return IntensitySample(uuid, mNeutrino[uuid].EnergyKin());
}

double Engine::IntensitySample(std::string uuid, double Energy)
{
	std::string name = "sample_" + uuid;
	return mDriver[uuid]->InterpolateIntensity(sampleNu[name], Energy);
}


double Engine::BinNumber()
{
	if (mDriver.begin()->second)
		return mDriver.begin()->second->BinNumber();
	else
		return -1.0;
}

double Engine::RangeWidth()
{
	double start, end;
	return RangeWidth(start, end);
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

