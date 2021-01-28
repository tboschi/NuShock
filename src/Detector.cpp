#include "detector/Detector.h"

Detector::Detector(const std::string &card) :
	_baseline(-1.),
	_length(-1.),
	_width(-1.),
	_height(-1.),
	_volume(-1.),
	_exposure(-1.),
	_weight(-1.)
{
	CardDealer cd(card);

	if (!cd.Get("default", _default)) {
		_default.clear();
		std::cerr << "Detector: WARNING no default module set!\n";
	}

	std::vector<std::string> keys = cd.ListKeys();
	for (std::string k : keys) {
		if (k.find_first_of('_') != std::string::npos)
			k.erase(k.find_first_of('_'));
		else
			continue;

		if (_modules.count(k) || k == "ratio")
			continue;
		// k is now detector module name
		Module detect;
		if (cd.Get(k + "_", detect))
			_modules[k] = detect;
		//else
			//std::cerr << "Detector: No configuration for module \"" << k << "\"\n";

		//std::string matk = k + "_material";
		std::string type;
		if (cd.Get(k + "_material", type))
			_materials[k] = Material::fromString(type);
		else {
			_materials[k] = Material::Air; // default material
			//std::cerr << "Detector: No material for module \"" << k << "\"\n";
		}
		if (cd.Get(k + "_shape", type))
			_shapes[k] = type;
		else {
			_shapes[k] = "box";	// default shape
			//std::cerr << "Detector: No shape for module \"" << k << "\"\n";
		}
	}
	if (!_materials.count("out"))
		_materials["out"] = Material::Lead;	// make outside of lead to stop everything

	// POTs is total number of POTs
	if (!cd.Get("POT", _POTs)) {
		double years;
		if (!cd.Get("years", years))
			throw std::logic_error("Detector: no \"years\" defined in configuration!\n");
		double pots;
		if (cd.Get("POT/s", pots))
			_POTs = 1.e7 * years * pots;
		else if (cd.Get("POT/y", pots))
			_POTs = years * pots;
		else
			throw std::logic_error("Detector: no \"POT\" or \"POT/s\""
							"or \"POT/y\" defined in configuration!\n");
	}

	if (!cd.Get("ratio_", _modes)) {
		_modes["FHC"] = 1.;
		_modes["RHC"] = 1.;
	}

	double sum = std::accumulate(_modes.begin(), _modes.end(), 0.,
			[](double sum, const std::pair<std::string, double> &m) {
				return sum + m.second; } );
	std::for_each(_modes.begin(), _modes.end(),
			[&](std::pair<const std::string, double> &m) {
				m.second /= sum; } );

	if (!cd.Get("beam", _Eb))
		_Eb = 0.;	// no beam energy specified
}

std::ostream & operator<<(std::ostream &os, const Detector &box) {
	os << "<";
	std::string sep = ", ";
	std::string del = "";
	for (const auto & im : box._modules) {
		os << del << im.first << " [" << Material::toString(box.MadeOf(im.first)) << "]";
		del = sep;
	}
	return os << ">";
}

/*
double Detector::Efficiency(Channel::Name chan, double mass, double energy) const
{
	if (_mass_ener_func.count(chan))
		return _mass_ener_func.at(chan)->Interpolate(mass, energy);
	else
		return 1.;
}

//load efficiency file
//key will be a combination such as CHANNEL_MODULE_FERMION
void Detector::LoadEfficiency(Channel::Name chan, std::string file)
{
	double rate = 1.;
	for (const auto &im : _modules)
		if (file.find(im.first) != std::string::npos) {
			if (_efficiencies.find(im.first) != _efficiencies.end()) {
				std::cerr << "Efficiency file for module " << im.first
					  << " already loaded\n";
				return;
			}
			rate = Weight(im.first) / Weight();
			break;
		}


	TFile infile(file.c_str());
	std::shared_ptr<TH2D> hhfunc(static_cast<TH2D*>(infile.Get("hhfunc")));
	hhfunc->SetDirectory(0);
	hhfunc->Scale(rate);
	if (!_mass_ener_func.count(chan))
		_mass_ener_func[chan] = hhfunc;
	else
		_mass_ener_func[chan]->Add(hhfunc.get());
}
*/

Material::Name Detector::MadeOf(std::string mod) const {
	if (!_materials.count(mod))
		return _materials.at(_default);
	return _materials.at(mod);
}

// flux scaling for a flux nu/ cm² / GeV / POTs @ 1 m
double Detector::Scaling() const {
	return Section() / std::pow(Baseline(), 2) * 1.e4;
}

double Detector::Scaling(std::string mod) const {
	return Section(mod) / std::pow(Baseline(mod), 2) * 1.e4;
}

// exposure in ton * MW * year / m² per POT unit!
double Detector::Exposure() const {
	if (_exposure < 0.)
		_exposure = std::accumulate(_modules.begin(), _modules.end(), 0.,
			[&](double sum, const std::pair<std::string, Module> &m) {
				return sum + Exposure(m.first);
			} );
	return _exposure;
}

double Detector::Exposure(std::string mod) const {
	return _Eb * Weight(mod) / std::pow(Baseline(mod), 2);
}

double Detector::Weight() const {
	if (_weight < 0.)
		_weight = std::accumulate(_modules.begin(), _modules.end(), 0.,
			[&](double sum, const std::pair<std::string, Module> &m) {
				return sum + Weight(m.first);
			} );
	return _weight;
}

double Detector::Weight(std::string mod) const {
	if (!_modules.count(mod) || !_modules.at(mod).count("weight"))
		return _modules.at(_default).at("weight");
	return _modules.at(mod).at("weight");
}

// baseline is the baseline of the closest module to the target
double Detector::Baseline() const {
	if (_baseline < 0.) {
		auto b = *std::min_element(_modules.begin(), _modules.end(),
				[&](const std::pair<std::string, Module> &a,
				   const std::pair<std::string, Module> &b) {
					return Baseline(a.first) < Baseline(b.first);
				});	
		_baseline = Baseline(b.first);
	}
	return _baseline;
}

double Detector::Baseline(std::string mod) const {
	if (!_modules.count(mod) || !_modules.at(mod).count("baseline"))
		return _modules.at(_default).at("baseline");
	return _modules.at(mod).at("baseline");
}

// full length of detector, from first baseline to last baseline+length
// it's not the sum of single components
double Detector::Length() const {
	if (_length < 0.) {
		auto l = *std::max_element(_modules.begin(), _modules.end(),
				[&](const std::pair<std::string, Module> &a,
				   const std::pair<std::string, Module> &b) {
					return Baseline(a.first) < Baseline(b.first);
				});	
		_length = Baseline(l.first) + Length(l.first) - Baseline();
	}
	return _length;
}

double Detector::Length(std::string mod) const {
	if (!_modules.count(mod) || !_modules.at(mod).count("length"))
		return _modules.at(_default).at("length");
	return _modules.at(mod).at("length");
}

// y direction
double Detector::Height() const {
	if (_height < 0.) {
	auto h = *std::max_element(_modules.begin(), _modules.end(),
			[&](const std::pair<std::string, Module> &a,
			   const std::pair<std::string, Module> &b) {
				return Height(a.first) < Height(b.first);
			});	
		_height = Height(h.first);
	}
	return _height;
}

double Detector::Height(std::string mod) const {
	if (!_modules.count(mod) || !_modules.at(mod).count("height"))
		return _modules.at(_default).at("height");
	return _modules.at(mod).at("height");
}

// x direction
double Detector::Width() const {
	if (_width < 0.) {
	auto h = *std::max_element(_modules.begin(), _modules.end(),
			[&](const std::pair<std::string, Module> &a,
			   const std::pair<std::string, Module> &b) {
				return Width(a.first) < Width(b.first);
			});	
		_width = Width(h.first);
	}
	return _width;
}

double Detector::Width(std::string mod) const {
	if (!_modules.count(mod) || !_modules.at(mod).count("width"))
		return _modules.at(_default).at("width");
	return _modules.at(mod).at("width");
}

double Detector::Section() const {
	return Height() * Width();
}

double Detector::Section(std::string mod) const {
	return Height(mod) * Width(mod);
}

double Detector::Radius() const {
	return std::sqrt(Section() / Const::pi);
}

double Detector::Radius(std::string mod) const {
	return std::sqrt(Section(mod) / Const::pi);
}

// upgrade to sum of individual volumes
double Detector::Volume() const
{
	if (_volume < 0.)
		_volume = std::accumulate(_modules.begin(), _modules.end(), 0.,
			[&](double sum, const std::pair<std::string, Module> &m) {
				return sum + Volume(m.first);
			});
	return _volume;
}

double Detector::Volume(std::string mod) const
{
	return Section(mod) * Length(mod);
}

double Detector::AngularAcceptance() const
{
	return std::atan2(Radius(), Baseline());
}

double Detector::AngularAcceptance(std::string mod) const
{
	return std::atan2(Radius(), Baseline(mod));
}

// check if a particle at origin moving with given theta and phi
// will enter the detector
bool Detector::AngularAccept(const Particle &P) const
{
	return AngularAccept(P.Theta(), P.Phi());
}

bool Detector::AngularAccept(double theta, double phi) const
{
	return std::any_of(_modules.begin(), _modules.end(),
			[&](const std::pair<std::string, Module> &m) {
				return AngularAccept(m.first, theta, phi);
			});
}

bool Detector::AngularAccept(std::string mod, const Particle &P) const
{
	return AngularAccept(mod, P.Theta(), P.Phi());
}

bool Detector::AngularAccept(std::string mod, double theta, double phi) const
{
	return theta <= AngularAcceptance(mod);
	//double x = Baseline(mod) * std::tan(theta) * std::cos(phi);
	//double y = Baseline(mod) * std::tan(theta) * std::sin(phi);
	//return (std::abs(x / Width(mod)) < 1.) && (std::abs(y / Height(mod)) < 1.);
}

double Detector::Probability(double tby) const {
	return std::accumulate(_modules.begin(), _modules.end(), 0.,
			[&](double sum, const std::pair<std::string, Module> &m) {
				return sum + Probability(m.first, tby);
			});
}

double Detector::Probability(std::string mod, double tby) const {
	return std::exp(- tby * Const::M2GeV * Baseline(mod))
		* (1 - std::exp(- tby * Const::M2GeV * Length(mod)));
}

std::vector<std::string> Detector::Modules() const {
	std::vector<std::string> mods;
	mods.reserve(_modules.size());
	for (const auto &m : _modules)
		mods.push_back(m.first);
	return mods;
}

std::string Detector::WhichModule(const Track &t) const {
	auto m = std::find_if(_modules.begin(), _modules.end(),
			[&](const std::pair<std::string, Module> &m) {
				return IsInside(m.first, t);
			});
	if (m == _modules.end())
		return "out"; // it is outside
	return m->first;
}

// check if inside minimum box containing detector
bool Detector::IsContained(const Track &t) const {
	// baseline is closest to origin, length is furthest
	// width and height are largest
	return (std::abs(t.Z() - Baseline() - Length() / 2.) <= Length() / 2.
	     && std::abs(t.Y()) <= Height() / 2
	     && std::abs(t.X()) <= Width() / 2);
}


bool Detector::IsInside(const Track &t) const {
	return std::any_of(_modules.begin(), _modules.end(),
			[&](const std::pair<std::string, Module> &m) {
				return IsInside(m.first, t);
			});
}

bool Detector::IsInside(std::string mod, const Track &t) const {
	if (!_shapes.count(mod))
		return IsInside(_default, t);

	double rz = t.Z() - Baseline(mod) - Length(mod) / 2.;
	double ry = t.Y();
	double rx = t.X();

	if (_shapes.at(mod) == "box")
		return (std::abs(rz) <= Length(mod) / 2.
		     && std::abs(ry) <= Height(mod) / 2.
		     && std::abs(rx) <= Width(mod) / 2.);

	if (_shapes.at(mod) == "tubx") {
		double rr = Height(mod) * Length(mod) / Const::pi;
		return (std::abs(rx) <= Width(mod) / 2.
		     && (ry * ry + rz * rz) <= rr);
	}

	if (_shapes.at(mod) == "tuby") {
		double rr = Width(mod) * Length(mod) / Const::pi;
		return (std::abs(ry) <= Height(mod) / 2.
		     && (rx * rx + rz * rz) <= rr);
	}

	if (_shapes.at(mod) == "tubz") {
		double rr = Width(mod) * Height(mod) / Const::pi;
		return (std::abs(rz) <= Length(mod) / 2.
		     && (rx * rx + ry * ry) <= rr);
	}

	if (_shapes.at(mod) == "sphere") {
		double rrr = 0.75 * Volume(mod) / Const::pi;
		return (rx * rx + ry * ry + rz * rz <= rrr);
	}

	return false;
}

double Detector::MagneticField(std::string mod) const {
	if (!_modules.count(mod) || !_modules.at(mod).count("B_field"))
		return _modules.at(_default).at("B_field");
	return _modules.at(mod).at("B_field");
}

double Detector::BeamEnergy() const {
	return _Eb;
}

double Detector::POTs() const {
	return _POTs;
}

double Detector::POTs(std::string hc) const {
	if (!_modes.count(hc))
		return _POTs;
	return _POTs * _modes.at(hc);
}
