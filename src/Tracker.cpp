#include "detector/Tracker.h"

Tracker::Tracker(std::string det, std::string track) : Detector(det)
{
	CardDealer cd(track);
	Init(track);
}

Tracker::Tracker(std::string card) : Detector(card)
{
	CardDealer cd(card);
	if (cd.Get("tracker_configuration", card))
		cd = CardDealer(card);

	Init(cd);
}

// init will set default values to make class constant
void Tracker::Init(const CardDealer &cd) {

	if (!cd.Get("threshold_", _thresholds))
		throw std::logic_error("No detection threshold specified");

	std::vector<std::string> fills = {"hadron", "muon", "pion", "gamma"};
	for (std::string f : fills)
		if (!_thresholds.count(f))
			_thresholds[f] = 0.;

	if (!cd.Get("resolution_", _resolutions))
		throw std::logic_error("No detection resolution specified");

	fills = {"hadron_angle", "hadron_energy", "hadron_enbias", "hadron_inrange", "hadron_exiting",
		 "muon_angle", "muon_energy", "muon_enbias", "muon_inrange", "muon_exiting",
		 "pion_angle", "pion_energy", "pion_enbias", "pion_inrange", "pion_exiting",
		 "gamma_angle", "gamma_energy", "gamma_enbias", "gamma_inrange", "gamma_exiting",
		 "containment", "vertex", "sagitta"};
	for (std::string f : fills)
		if (!_resolutions.count(f))
			_resolutions[f] = 0.;

	if (!cd.Get("identify_", _identifies))
		std::cerr << "No identification parameters set; cannot do background generation\n";

	fills = {"conversion_EM", "length_MIP", "pair_angle"};
	for (std::string f : fills)
		if (!_identifies.count(f))
			_identifies[f] = 0.;

	if (!cd.Get("track_step", _step))
		_step = 0.01; // 1 cm

	if (!cd.Get("verbosity", kVerbose))
		kVerbose = false;
}

// transform a particle into an event inside the detector
Tracker::Event Tracker::GenerateEvent(Particle &&P) const
{
	// pick random module
	std::uniform_real_distribution<> rndm(0, Weight());
	double w = rndm(RNG::_mt);
	auto imod = std::find_if(_modules.begin(), _modules.end(),
			[&](const std::pair<std::string, Module> &m) {
				bool ret = w < Weight(m.first);
				w -= Weight(m.first);
				return ret;
				} );
			       	
	return GenerateEvent(imod->first, std::move(P));
	// particle is in module in_mod
}

Tracker::Event Tracker::GenerateEvent(std::string mod, Particle &&P) const
{
	std::uniform_real_distribution<> rndm(-1., 1.);
	std::normal_distribution<> gaus(0., 1. / 3.);

	double x = std::min(1., std::max(-1., gaus(RNG::_mt)));
	double y = std::min(1., std::max(-1., gaus(RNG::_mt)));
	double z = rndm(RNG::_mt);
	if (_shapes.at(mod) == "tubx") {
		z *= std::sqrt(1 - y*y);
	}
	else if (_shapes.at(mod) == "tuby") {
		z *= std::sqrt(1 - x*x);
	}
	else if (_shapes.at(mod) == "sphere") {
		z *= std::sqrt(1 - x*x - y*y);
	}
	//else if (_shapes.at(mod) == "box" || _shapes.at(mod) == "tubz")
		// nothing to be done
		
	//std::cout << "random " << x << ", " << y << ", " << z << "\n";
	x *= Width(mod) / 2.;
	y *= Height(mod) / 2.;
	z = (z + 1.) * Length(mod) / 2. + Baseline(mod);

	Track T(x, y, z);
	if (!IsDetectable(T))
		//std::5out << "T is out! " << T << "\n";
		return GenerateEvent(mod, std::move(P));

	// 99.99% of particles are inside the detector
	P.SetTheta(std::atan2(std::sqrt(x*x + y*y), z));
	P.SetPhi(std::atan2(y, x));

	return std::make_pair(P, T);
}

// check if event is valid: position of the particle, threshold
bool Tracker::Reconstruct(Event &evt) const
{
	//std::cout << "reconst \t" << evt.first << "\n        \t" << evt.second << "\n";
	//std::cout << "event detect (" << std::boolalpha
	//	  << IsDetectable(evt.first) << ", " << IsDetectable(evt.second)
	//	  << ")\n";
	if (IsDetectable(evt))	//valid event
	{
		Particle &P = evt.first;
		Track &T = evt.second;

		// if track is not set..
		if (T.Length() <= 0.)
			ComputeTrack(evt);

		//double sag = ComputeSagitta(evt);
		//if (sag > _resolutions.at("sagitta"))
		if (DetectSagitta(evt))
			P.ChargeID(); // can charge ID correctly
		else
			P.ChargeID(0.);

		// do some smearing
		Smearing(evt);

		// check just threshold again
		return (IsDetectable(P));
	}

	return false;
}

void Tracker::ComputeTrack(Event &evt) const	//This should not change P
{
	double tmax, depth;

	// copy particle
	Particle P = evt.first;
	// but reference to track, as this is updated
	Track &T = evt.second;
	// copy original position
	TVector3 pos = static_cast<TVector3>(evt.second);
	std::string mod = WhichModule(T);
	Material::Name mat = MadeOf(mod);
	std::normal_distribution<> gaus;
	std::exponential_distribution<> expd;
	std::uniform_real_distribution<> rndm;

	double cover = 0;
	double radlen = 0;
	int layer = 0;

	switch (std::abs(P.Pdg()))
	{
		case 11:
			tmax = std::log( P.E() / Material::CriticalEnergy(mat) ) - 1.0;
			gaus = std::normal_distribution<>(tmax, tmax / 3.);
			depth = Material::RadiationLength(mat) * gaus(RNG::_mt) * 0.01;
			T.SetLength(mod, depth);
			T.SetEnergyDeposited(mod, P.E());
			break;
		case 22:
			depth = 9./7. * Material::RadiationLength(mat) * 0.01;
			expd = std::exponential_distribution<>(1. / depth);
			T.SetLength(mod, expd(RNG::_mt));
			T.SetEnergyDeposited(mod, P.E());
			break;
		case 13: //quit when particle decays too!
			gaus = std::normal_distribution<>(_step, _resolutions.at("vertex"));
			while (IsDetectable(evt) && !IsDecayed(P, T.Length())) {
				double dx = gaus(RNG::_mt);	// random step
				double loss = dx * Material::Bethe(MadeOf(mod), P.M(), P.Beta(), P.Gamma());
				P.SetE(P.E() - loss); // reduce energy

				T.SetLength(mod, T.Length(mod) + dx);
				T.SetEnergyDeposited(mod, T.EnergyDeposited(mod) + loss);
				T += dx * P.Vect().Unit();

				// update position of particle
				if (!IsInside(mod, T))
					mod = WhichModule(T);
			}
			break;
		case 211:
			depth = Material::NuclearRadiationLength(mat) * 0.01;
			gaus = std::normal_distribution<>(_step, _resolutions.at("vertex"));
			expd = std::exponential_distribution<>(1. / depth);
			while (IsDetectable(evt) && !IsDecayed(P, T.Length())) {
				if (cover > radlen)	//distance more than interaction length
				{	//multiplicity of hadronic interaction
					int mult = std::ceil(std::pow(P.EKin() * 1e6, 0.18) * 0.15);
					if (mult > 1)
						++layer;

					radlen = expd(RNG::_mt);
					cover = 0;
				}


				double dx = gaus(RNG::_mt);	// random step
				double loss = dx * Material::Bethe(MadeOf(mod), P.M(), P.Beta(), P.Gamma());
				P.SetE(P.E() - loss); // reduce energy

				// distance covered so far in this layer
				cover += dx;

				T.SetLength(mod, T.Length(mod) + dx);
				T.SetEnergyDeposited(mod, T.EnergyDeposited(mod) + loss);
				T += dx * P.Vect().Unit();

				// update position of particle
				if (!IsInside(mod, T))
					mod = WhichModule(T);
			}

			T.SetShower((rndm(RNG::_mt) > 1./layer));
			break;
		default:
			break;
	}
	T.SetPosition(pos);
}

// return true if sagitta can be computed in any module
bool Tracker::DetectSagitta(Event &evt) const {
	const Particle &P = evt.first;
	Track &T = evt.second;

	if (P.RealQ() == 0)
		return false;

	for (const auto &m : _modules) {
		auto B = MagneticField(m.first);
		double sag = 0.;
		std::array<double, 3> p = {{P.Px(), P.Py(), P.Pz()}};
		for (int i = 0; i < 3; ++i) {
			if (B[i] > 0) {
				// perpendicular P component wrt to B
				double pp = std::sqrt(std::pow(p[(i+1)%3], 2) + std::pow(p[(i+2)%3], 2));
				double s = 0.3 / 8. * std::pow(T.Length(m.first), 2)
					* P.RealQ() * B[i] / pp;
				sag += s * s;
			}
		}

		if (std::sqrt(sag) * _resolutions.at("sagitta") > _resolutions.at("vertex"))
			return true;
	}

	return false;
}

/*
double Tracker::ComputeSagitta(Event &evt) const
{
	const Particle &P = evt.first;
	Track &T = evt.second;

	if (P.RealQ() == 0)
		return 0;

	double sag = 0;
	double max_B = 0.;
	// get track of modules with magnetic field
	for (const auto &m : _modules) {
		//std::cout << "in " << m.first << " has " << T.Length(m.first) << "\n";
		if (MagneticField(m.first) <= 0.)
			continue;

		//double track2 = std::pow(T.Length(m.first) * std::cos(P.Theta()), 2) +
				//std::pow(T.Length(m.first) * P.Theta() * P.Phi(), 2);
		sag += T.Length(m.first);

		if (max_B < MagneticField(m.first))
			max_B = MagneticField(m.first);
	}
	
	//double rad = std::abs(0.3 * P.RealCharge()) * minB / P.Momentum());

	return std::abs(0.3 * P.RealQ() * max_B * sag / P.P())
	    * (std::pow(std::cos(P.Theta()), 2)
	     + std::pow(std::sin(P.Theta()) * std::sin(P.Phi()), 2));

	//without saggitta approximation
	//double rad = 0.1 * std::abs(P.RealCharge()) *
	//	     Get("MagneticField") / P.Momentum();
	//double track = P.TrackFGT() * sqrt(pow(cos(P.Theta()), 2) +
	//				   pow(P.Theta() * P.Phi(), 2));
	//double sag = rad * ( 1 - cos(track/rad) );
}
*/


// smear the track
void Tracker::Smearing(Event &evt) const
{
	Particle &P = evt.first;
	Track &T = evt.second;

	double sigma = 0.;
	std::string angle;

	//double perc = 0.;
	//bool hires = false;

	std::normal_distribution<> gaus;
	switch (std::abs(P.Pdg()))
	{
		case 11:
		case 22:
			angle = "gamma_angle";

			sigma = std::sqrt(std::pow(_resolutions.at("gamma_energy"), 2) / P.EKin()
					+ std::pow(_resolutions.at("gamma_enbias"), 2));
			gaus = std::normal_distribution<>(P.EKin(), sigma * P.EKin());
			P.SetEKin(gaus(RNG::_mt));
	
			break;
		case 13:
			angle = "muon_angle";

			// high resolution if most of track is in LAr
			// or module w/ magnetic field was crossed
			/*
			for (const auto &m : _modules)
				perc += T.Length(m.first) * T.DepositedEnergy(m.first);
				if (T.Length(m.first)
				if (MagneticField(m.first) > 0. && T.Length(m.first) > 0.) {
					hires = true;
					break;
				}
				if (MadeOf(m.first) == Material::LAr)
					perc += T.Length(m.first);
			}
			hires = hires || (1 - perc / T.Length()) > _resolutions.at("containment");
			*/

			if (1 - T.Importance("out") / T.Importance() > _resolutions.at("containment"))
				sigma = _resolutions.at("muon_inrange");
			else
				sigma = _resolutions.at("muon_exiting");

			gaus = std::normal_distribution<>(P.P(), sigma * P.P());
			P.SetP(gaus(RNG::_mt));
			break;
		case 211:
			angle = "pion_angle";

			// high resolution if most of track is in LAr
			// or module w/ magnetic field was crossed
			/*
			for (const auto &m : _modules) {
				if (MagneticField(m.first) > 0. && T.Length(m.first) > 0.) {
					hires = true;
					break;
				}
				if (MadeOf(m.first) == Material::LAr)
					perc += T.Length(m.first);
			}
			hires = hires || (1 - perc / T.Length()) > _resolutions.at("containment");
			*/

			//if (1 - T.Length("out") / T.Length() > _resolutions.at("containment")) {
			if (1 - T.Importance("out") / T.Importance() > _resolutions.at("containment"))
				sigma = _resolutions.at("pion_inrange");
			else
				sigma = _resolutions.at("pion_exiting");

			gaus = std::normal_distribution<>(P.P(), sigma * P.P());
			P.SetP(gaus(RNG::_mt));
			break;
		default:	//I don't care about other particles
			return;
	}

	// do angles
	gaus = std::normal_distribution<>(P.Theta(), _resolutions.at(angle) / Const::Deg);
	P.SetTheta(gaus(RNG::_mt));
	gaus = std::normal_distribution<>(P.Phi(), _resolutions.at(angle) / Const::Deg);
	P.SetPhi(gaus(RNG::_mt));
}



bool Tracker::IsDecayed(const Particle &P, double dx) const	//Threshold check
{
	std::uniform_real_distribution<> rndm;
	return rndm(RNG::_mt) > std::exp(-dx / (Const::C * P.LabSpace()));
}

// check if inside detector or above threshold
bool Tracker::IsDetectable(const Event &evt) const {
	// check threshold		// check position
	return IsDetectable(evt.first) && IsDetectable(evt.second);
}

bool Tracker::IsDetectable(const Track &T) const {
	return IsContained(T);
}

bool Tracker::IsDetectable(const Particle &P) const
{
	double threshold = 0.0;
	std::string part;
	std::uniform_real_distribution<> rndm;
	switch (std::abs(P.Pdg()))
	{
		case 11:
		case 22:
			threshold = _thresholds.at("gamma");
			break;
		case 13:
			threshold = _thresholds.at("muon");
			break;
		case 211:
			threshold = _thresholds.at("pion");
			break;
		case 311:	//neutral kaons
		case 2112:	//neutrons
			if (rndm(RNG::_mt) < 0.1)	//90% efficiency
				return false;
			// else is hadron det threshold
		case 2212:	//protons
			threshold = _thresholds.at("hadron");
			break;
		default:
			if (P.RealQ() != 0)
				threshold = _thresholds.at("hadron");
			else
				return false;
			break;
	}
	return (P.EKin() > threshold);
}

// special treatment for pi0
std::array<Tracker::Event, 2> Tracker::Pi0Decay(Tracker::Event &&pi0) const
{
	Particle &P0 = pi0.first;
	if (std::abs(P0.Pdg()) != 111)
		throw std::invalid_argument("Tracker: given particle is not a PI0");
	Track &T0 = pi0.second;

	Particle gammaA(22, P0.M()/2.0, 0, 0,  P0.M()/2.0); 
	Particle gammaB(22, P0.M()/2.0, 0, 0, -P0.M()/2.0); 

	std::uniform_real_distribution<> angle(-Const::pi, Const::pi);
	double theta = angle(RNG::_mt);
	double phi = angle(RNG::_mt);

	gammaA.SetTheta(theta);
	gammaB.SetTheta(theta + Const::pi);

	gammaA.SetPhi(phi);
	gammaB.SetPhi(phi + Const::pi);

	gammaA.Boost(P0.BoostVector());
	gammaB.Boost(P0.BoostVector());

	// gamma decay here
	double depth = 9./7. * Material::RadiationLength(MadeOf(WhichModule(T0))) * 0.01;
	std::exponential_distribution<> expd(1. / depth);
	//

	Track moveA = T0 + T0.Unit() * expd(RNG::_mt);
	Track moveB = T0 + T0.Unit() * expd(RNG::_mt);

	Tracker::Event pA{gammaA, moveA};
	Tracker::Event pB{gammaB, moveB};
	return {pA, pB};
}

// loop through events and add some chaos
std::vector<Tracker::Event> Tracker::MisIdentify(std::vector<Tracker::Event> &&evts) const
{
	if (!_identifies.size())
		throw std::logic_error("No identification parameter set\n");

	std::vector<Tracker::Event> electrons;

	for (auto ie = evts.begin(); ie != evts.end(); )
	{
		if (!IsDetectable(*ie)) {
			ie = evts.erase(ie);
			continue;
		}

		bool elec = true;
		Particle &P = ie->first;
		Track &T = ie->second;
		switch (std::abs(P.Pdg()))
		{
			case 13:	// nothing to do for muon
				break;
			case 211:
				if (T.Length() - T.Length("out") > _identifies.at("length_MIP")
				 && !T.IsShower())
				{	//it looks like a muon
					P.SetPdg(P.Q() / 3 * 13);
					P.SetM(Const::MMuon);
					--ie;		//must recheck (for thr for instance)
				}
				break;
			case 11:
				elec = true;
				for (auto ee = electrons.begin(); ee != electrons.end(); ) {
					Particle &eP = ee->first;
					Track &eT = ee->second;

					double angle = P.Vect().Angle(eP.Vect());
					if (angle >= _identifies.at("pair_angle") / Const::Deg) {
						++ee;
						continue;
					}

					//two track are too close, it can be a pair production);
					Track dir = Track(T.X(), T.Y(), T.Z());
					double T_lin = T.Length() - T.Length("out");
					double eT_lin = eT.Length() - eT.Length("out");
					double dist = std::sqrt(std::pow(T_lin, 2) + std::pow(eT_lin, 2)
						    + 2. * T_lin * eT_lin * std::cos(angle));
					dir.SetLength("in", dist);

					Particle gamma(22, P + eP);
					gamma.SetM(0);

					*ie = std::make_pair(gamma, dir);
					--ie;		//must recheck (for thr for instance)

					ee = electrons.erase(ee);
					elec = false;
					break;
				}

				if (elec)
					electrons.push_back(*ie);
				break;
			case 22:
				//short displacement > 2cm
				if((T.Length() - T.Length("out")) < _identifies.at("conversion_EM"))
				{	//conversion before 3cm 
					P.SetPdg(std::distance(ie, evts.begin()) % 2 == 0 ? 11 : -11);
					P.SetM(Const::MElectron);

					--ie;		//must recheck (for thr for instance)
				}
				break;
		}

		++ie;
	}

	return evts;
}
