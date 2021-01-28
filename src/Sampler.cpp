#include "flux/Sampler.h"

Sampler::Sampler(const Detector &box, const Driver &drive) :
		//Neutrino nu, Mixing mix) :
	_box(box),
	_drive(drive)
	//_mix(std::move(mix))	// mix is blueprint
{
}
/*
	// create neutrino and antineutrino versions of nu
	if (nu.IsDiracParticle()) {
		_N0 = std::move(nu);
		size_t optB = (nu.GetOptions() & ~Neutrino::antiparticle) | Neutrino::antiparticle;
		_NB = Neutrino(nu.M(), optB);
	}
	else { // isDiracAntiparticle
		_NB = std::move(nu);
		size_t opt0 = (nu.GetOptions() & ~Neutrino::antiparticle) | Neutrino::particle;
		_N0 = Neutrino(nu.M(), opt0);
	}

	_mix.Flatten();
}
*/

// do this check beforehand
bool Sampler::Bind(Channel::Name chan, Mixing mix, Neutrino N) {
	_N    = std::move(N);
	_chan = std::move(chan);

	_rate.SetNeutrino(_N);

	mix.Flatten();
	return _rate.IsAllowed(_chan) && _drive.MakeFlux(_N, mix);
}

/*
bool Sampler::IsAllowed(Channel::Name chan, std::initializer_list nus, const Mixing &mix) {
		// production		// decays
	for (const auto &n : nus) {
		DecayRates r(_N0);
		if (!r.IsAllowed(chan))
			return false;
	}

}
*/

// passing just one variable
/*
std::shared_ptr<TH1D> Sampler::MakeSampler(Channel::Name chan, Neutrino N, double x)
{
	//Mixing mix = _mix * x;
	return MakeSampler(chan, _mix * x);	
}
*/

std::shared_ptr<TH1D> Sampler::MakeSampler(const Mixing &mix)
{
	std::shared_ptr<TH1D> spectrum = _drive.Spectrum(_N, mix);
	if (!spectrum)
		return nullptr;

	std::cout << "using spectrum " << spectrum->GetName() << "\t" << spectrum->GetTitle() << "\n";
	spectrum->Scale(_box.Scaling());

	for (int bin = 1; bin <= spectrum->GetNbinsX(); ++bin) {
		_N.SetE(spectrum->GetBinCenter(bin)); // set energy
		double cont = spectrum->GetBinContent(bin) * spectrum->GetBinWidth(bin);
		if (_N.EKin() <= 0. || _N.Beta() > 1.) { // unphysical neutrino
			spectrum->SetBinContent(bin, 0.);
			continue;
		}
		spectrum->SetBinContent(bin, cont //* _box.Efficiency(chan, N.M(), N.E())
			* _box.Probability(_rate.Total(mix) / (_N.Beta() * _N.Gamma())) );

		//std::cout << "spec " << spectrum->GetBinLowEdge(bin)
		//	  << "\t" << cont  << "\t" << rate.Total(mix)
		//	  << "\t" << rate.Branch(chan, mix)
		//	  << "\t" << N.Beta() * N.Gamma()
		//	  << "\t" << _box.Probability(rate.Total(mix) / (N.Beta() * N.Gamma()))
		//	  //<< "\t" << _box.Efficiency(chan, N.M(), N.E())
		//	  << "\t = " << spectrum->GetBinContent(bin) << "\n";
	}
	spectrum->Scale(_rate.Branch(_chan, mix));

	return spectrum;
}


/*
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
*/
