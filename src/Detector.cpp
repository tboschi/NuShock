#include "Detector.h"

Detector::Detector(const std::string &card) :
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
		if (_modules.count(k))
			continue;
		// k is now detector module name
		Modules detect;
		if (cd.Get(k, detect)) {
			_modules[k] = detect;
			std::string matk = k + "_material";
			std::string type;
			if (cd.Get(mat, type))
				_materials[k] = Material::fromString(type);
		}
		else if
			std::cerr << "Detector: No configuration for module \"" << k << "\"\n";
	}

	if (!cd.Get("POT", _POTs)) {
		double pots;
		if (!cd.Get("POT/s", _POTs))
			throw std::logic_error("Detector: no \"POTs\" or \"POT/s\" defined in configuration!\n");
		double years;
		if (!cd.Get("years", years))
			throw std::logic_error("Detector: no \"years\" defined in configuration!\n");
		_POTs = 1.e7 * years * pots;
	}

	if (!cd.Get("eff_", _efficiencies))
		std::cerr << "Detector: No efficiency files found\n";
}

double Detector::Efficiency(Channel::Name chan, double mass, double energy)
{
	if (_mass_ener_func.count(chan))
		return _mass_ener_func[chan]->Interpolate(mass, energy);
	else
		return 1.;
}

//load efficiency file
//key will be a combination such as CHANNEL_MODULE_FERMION
void Detector::LoadEfficiency(std::string type, Channel::Name chan)
{
	if (!_efficiencies.count(type)) {
		std::cout << "No efficiency file for \"" << type << "\"\n";
		return;
	}

	double rate = 1.;
	for (const auto &im : _modules)
		if (type.find(im.first) != std::string::npos) {
			rate = Weight(im.first) / Weight();
			break;
		}


	TFile infile(_efficiencies[type].c_str());
	TH2D* hhfunc = static_cast<TH2D*>(infile.Get("hhfunc"));
	hhfunc->SetDirectory(0);
	hhfunc->Scale(rate);
	if (!_mass_ener_func.count(chan))
		_mass_ener_func[chan] = hhfunc;
	else
		_mass_ener_func[chan]->Add(hhfunc);
	if (!mhFunc.count(channel))
	{
		mhFunc[channel] = dynamic_cast<TH2D*> (hist->Clone());
		mhFunc[channel]->SetDirectory(0);
		mhFunc[channel]->Scale(rat);
	}
	else	//adding module with rate
		mhFunc[channel]->Add(hist, rat);
}

double Detector::Material(cstd::string mod) {
	if (_modules[mod].count("material"))
		return _modules[_default]["material"];
}

double Detector::Exposure() {
	return std::accumulate(_modules.begin(), _modules.end(), 0.,
			[](double sum, const std::pair<std::string, Module> &m) {
				return sum + Exposure(m.first);
			} );
}

double Detector::Exposure(std::string mod) {
	return _POTs * Area(mod) / std::pow(Baseline(mod), 2);
}

double Detector::Weight() {
	return std::accumulate(_modules.begin(), _modules.end(), 0.,
			[](double sum, const std::pair<std::string, Module> &m) {
				return sum + Weight(m.first);
			} );
}

double Detector::Weight(std::string mod) {
	if (_modules[mod].count("weight"))
		return _modules[_default]["weight"];
	return _modules[mod]["weight"];
}

// baseline is the baseline of the closest module to the target
double Detector::Baseline() {
	return std::min_element(_modules.begin(), _modules.end(),
			[](const std::pair<std::string, Module> &a,
			   const std::pair<std::string, Module> &b) {
				return Baseline(a.first) < Baseline(b.first);
			});	
}

double Detector::Baseline(std::string mod) {
	if (_modules[mod].count("baseline"))
		return _modules[_default]["baseline"];
	return _modules[mod]["baseline"];
}

// full length of detector
// sum of single components
double Detector::Length() {
	auto l = *std::min_element(_modules.begin(), _modules.end(),
			[](const std::pair<std::string, Module> &a,
			   const std::pair<std::string, Module> &b) {
				return Baseline(a.first) < Baseline(b.first);
			});	
	return Length(l.first);
}

double Detector::Length(std::string mod) {
	if (_modules[mod].count("length"))
		return _modules[_default]["length"];
	return _modules[mod]["length"];
}

double Detector::Height() {
	auto h = *std::max_element(_modules.begin(), _modules.end(),
			[](const std::pair<std::string, Module> &a,
			   const std::pair<std::string, Module> &b) {
				return Height(a.first) < Height(b.first);
			});	
	return Height(h.first);
}

double Detector::Height(std::string mod) {
	if (_modules[mod].count("height"))
		return _modules[_default]["height"];
	return _modules[mod]["height"];
}

double Detector::Width() {
	auto h = *std::max_element(_modules.begin(), _modules.end(),
			[](const std::pair<std::string, Module> &a,
			   const std::pair<std::string, Module> &b) {
				return Width(a.first) < Width(b.first);
			});	
	return Width(h.first);
}

double Detector::Width(std::string mod) {
	if (_modules[mod].count("height"))
		return _modules[_default]["height"];
	return _modules[mod]["height"];
}

double Detector::Area() {
	return Height() * Width();
}

double Detector::Area(std::string mod) {
	return Height(mod) * Width(mod);
}

double Detector::Volume()
{
	return Area() * Length();
}

double Detector::Volume(std::string mod)
{
	return Area(mod) * Length(mod);
}

double Detector::AngularAcceptance()
{
	// area / pi is "radius"
	return std::atan2(Area() / Const::Pi, Baseline());
}

double Detector::AngularAcceptance(std::string mod)
{
	// area / pi is "radius"
	return std::atan2(Area(mod) / Const::Pi, Baseline(mod));
}

double Detector::Probability(double tby) {
	return std::accumulate(_modules.begin(), _modules.end(), 0.,
			[](double sum, const std::pair<std::string, Module> &m) {
				return sum + Probability(mod, tby);
			});
}

double Detector::Probability(std::string mod, double tby) {
	return std::exp(- tby * Baseline(mod)) * (1 - std::exp(- tby * Length(mod)));
}

bool Detector::IsInside(const Track &t) {
	return std::any_of(_modules.begin(), _modules.end(),
			[&](const std::pair<std::string, Module> &m) {
				return IsInside(m.first, t);
			});
}

bool Detector::IsInside(std::string mod, const Track &t) {
	return (std::abs(t.Z() - Baseline(mod) - Length(mod) / 2.) <= Length(mod) / 2.
	     && std::abs(t.Z() - Width(mod) / 2.) <= Width(mod) / 2.
	     && std::abs(t.Z() - Height(mod) / 2.) <= Height(mod) / 2.);
}
