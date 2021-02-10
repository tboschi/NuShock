#include "detector/Driver.h"

/*
// return spectrum based on flavour
std::shared_ptr<TH1D> Driver::Spectrum(Nu::Flavour flv, const Mixing &mix) const {

	if (!_distr.count(flv))	// skip null histograms
		return nullptr;
	if (!_distr.at(flv))
		return nullptr;

	auto spec = std::shared_ptr<TH1D>(static_cast<TH1D*>(_distr.at(flv)->Clone()));
	spec->Scale(mix(flv, 2));
	return spec;
}
*/


Driver::Driver(const std::string &card) :
	_mass(-1.),
	_helix(0)
{
	CardDealer cd(card);

	// create production object for light neutrino
	//Neutrino nu();

	Init(cd);
}

// load fluxes for neutrinos and antineutrinos
void Driver::Init(const CardDealer &cd)
{
	std::map<std::string, std::string> files;
	if (!cd.Get("flux_", files))
		throw std::invalid_argument("Driver: flux card does not contain any flux description\n");

	for (const auto &ff : files)
		_fxNu[Nu::fromString(ff.first)] = std::move(Flux::fromROOT(ff.second));

	// load tau flux modifiers
	if (!cd.Get("mod_", files))
		return;

	for (const auto &id : files) {
		Driver::Modifier mod;
		std::string line;
		std::ifstream inmod(id.second);
		while (std::getline(inmod, line)) {
			if (line.find_first_of('#') != std::string::npos)
				line.erase(line.find_first_of('#'));

			if (line.empty())
				continue;
			std::array<double, 4> entry;
			std::stringstream ssl(line);
			ssl >> entry[0] >> entry[1] >> entry[2] >> entry[3];
			mod.push_back(std::move(entry));
		}
		_modifiers[Flux::fromString(id.first)] = std::move(mod);
	}
}

/*
// return spectrum based on neutrino (majorana vs particle vs antiparticle)
std::shared_ptr<TH1D> Driver::Spectrum(const Neutrino &N, const Mixing &mix) const {
	if ( _mass != N.M() || (_helix != N.Helicity() )  // asking for a different spectrum
		throw std::invalid_argument("Driver: flux was not made for this neutrino."
				" Run MakeFlux first?");

	std::vector<Nu::Flavour> flvs;
	if (N.IsMajorana())
		flvs = Nu::All();
	else if (N.IsParticle())
		flvs = Nu::Part();
	else // is antiparticle
		flvs = Nu::Anti();

	std::shared_ptr<TH1D> spec = nullptr;
	for (const auto &flv : flvs) {
		if (!_distr.count(flv))
			continue;
		if (!_distr.at(flv))	// skip null histograms
			continue;
		if (!spec) {	// clone first non null and scale to mix
			spec = std::shared_ptr<TH1D>(static_cast<TH1D*>(_distr.at(flv)->Clone()));
			spec->Scale(mix(flv, 2));
		}
		else	// add the rest
			spec->Add(_distr.at(flv).get(), mix(flv, 2));
	}

	return spec;
}

// return spectrum based on flavour
std::shared_ptr<TH1D> Driver::Spectrum(Nu::Flavour flv, const Mixing &mix) const {

	if (!_distr.count(flv))	// skip null histograms
		return nullptr;
	if (!_distr.at(flv))
		return nullptr;

	auto spec = std::shared_ptr<TH1D>(static_cast<TH1D*>(_distr.at(flv)->Clone()));
	spec->Scale(mix(flv, 2));
	return spec;
}
*/

// mix defaults to {1, 1, 1}
/*
bool Driver::MakeFlux(std::initializer_list<Neutrino> nus, const Mixing &mix) {
	bool ret = false;
	for (const Neutrino &N : nus)
		ret |= MakeFlux(N, mix);
	return ret;
}
*/

// make all components depending on neutrino type
// production is independent of majorana/dirac but depends on helicity
Spectrum Driver::MakeSpectrum(const Neutrino &N, const Mixing &mix) const
{	
	ProductionRate hnl(N);

	std::vector<Nu::Flavour> flvs;
	if (N.IsMajorana())
		flvs = Nu::All();
	else if (N.IsParticle())
		flvs = Nu::Part();
	else // is antiparticle
		flvs = Nu::Anti();

	Spectrum::Distribution distr;

	for (const auto &flv : Nu::All())
		if (mix(flv))
			distr[flv] = MakeComponent(hnl, flv, N.M()); 

	// spectrum object
	return Spectrum{std::move(N), std::move(distr)};
}

std::shared_ptr<TH1D> Driver::MakeComponent(ProductionRate &hnl, Nu::Flavour nu, double mass) const
{
	if (!_fxNu.count(nu))
		return nullptr;

	switch (nu) {
		case Nu::E0:
		case Nu::EB:
			return MakeElectron(hnl, _fxNu.at(nu));
		case Nu::M0:
		case Nu::MB:
			return MakeMuon(hnl, _fxNu.at(nu));
		case Nu::T0:
		case Nu::TB: // only tau component requires neutrino mass
			return MakeTau(hnl, _fxNu.at(nu), mass);
		default:
			return nullptr;
	}
}

std::shared_ptr<TH1D> Driver::MakeElectron(ProductionRate &heavy, const Flux::Component &fxNu) const
{
	//auto fxHNL = fxNu;
	// start with empty map
	// only allowed processes will be added
	// if nothing is allowed a null pointer is returned
	Flux::Component fxHNL;

	//pi+ -> e+ nu_e
	if (fxNu.count(Flux::Parent::Pion)
	&& heavy.IsAllowed(Production::Channel::PionE)) {
		fxHNL[Flux::Parent::Pion] = std::shared_ptr<TH1D>(static_cast<TH1D*>
					(fxNu.at(Flux::Parent::Pion)->Clone()));
		fxHNL[Flux::Parent::Pion]->Scale(heavy.Scale(Production::Channel::PionE));
	}

	//K+ -> pi0 e+ nu_e	(5.07 %)	//I just assume that both decays are equally probable 
	//K+ -> e+ nu_e		(1.582e-3)	//for each different energy, so it is a linear combination
	if (fxNu.count(Flux::Parent::Kaon)
	&& (heavy.IsAllowed(Production::Channel::KaonE) || heavy.IsAllowed(Production::Channel::KaonCE))) {
		fxHNL[Flux::Parent::Kaon] = std::shared_ptr<TH1D>(static_cast<TH1D*>
					(fxNu.at(Flux::Parent::Kaon)->Clone()));
		double kf = (1.582e-3 * heavy.Scale(Production::Channel::KaonE) 
			       + 5.07 * heavy.Scale(Production::Channel::KaonCE))
			  / (1.582e-3 + 5.07);
		fxHNL[Flux::Parent::Kaon]->Scale(kf);
	}

	//K0 -> pi+ e+ nu_e
	if (fxNu.count(Flux::Parent::Kaon0)
	&& heavy.IsAllowed(Production::Channel::Kaon0E)) {
		fxHNL[Flux::Parent::Kaon0] = std::shared_ptr<TH1D>(static_cast<TH1D*>
					(fxNu.at(Flux::Parent::Kaon0)->Clone()));
		fxHNL[Flux::Parent::Kaon0]->Scale(heavy.Scale(Production::Channel::Kaon0E));
	}

	//mu+ -> nu_mu_bar e+ nu_e
	if (fxNu.count(Flux::Parent::Muon)
	&& heavy.IsAllowed(Production::Channel::MuonE)) {
		fxHNL[Flux::Parent::Muon] = std::shared_ptr<TH1D>(static_cast<TH1D*>
					(fxNu.at(Flux::Parent::Muon)->Clone()));
		fxHNL[Flux::Parent::Muon]->Scale(heavy.Scale(Production::Channel::MuonE));
	}

	//Ds+ -> e+ nu_e
	if (fxNu.count(Flux::Parent::Charm)
	&& heavy.IsAllowed(Production::Channel::CharmE)) {
		fxHNL[Flux::Parent::Charm] = std::shared_ptr<TH1D>(static_cast<TH1D*>
					(fxNu.at(Flux::Parent::Charm)->Clone()));
		fxHNL[Flux::Parent::Charm]->Scale(heavy.Scale(Production::Channel::CharmE));
	}

	if (!fxHNL.size())
		return nullptr;

	// combine into total
	return std::accumulate(std::next(fxHNL.begin()), fxHNL.end(), fxHNL.begin()->second,
			[](std::shared_ptr<TH1D> total,
				const std::pair<Flux::Parent, std::shared_ptr<TH1D> > &hist) {
				total->Add(hist.second.get());
				return total; } );
}

std::shared_ptr<TH1D> Driver::MakeMuon(ProductionRate &heavy, const Flux::Component &fxNu) const
{
	// start with empty map
	// only allowed processes will be added
	// if nothing is allowed a null pointer is returned
	Flux::Component fxHNL;
	//auto fxHNL = fxNu;

	//pi+ -> mu+ nu_mu
	if (fxNu.count(Flux::Parent::Pion)
	&& heavy.IsAllowed(Production::Channel::PionM)) {
		fxHNL[Flux::Parent::Pion] = std::shared_ptr<TH1D>(static_cast<TH1D*>
					(fxNu.at(Flux::Parent::Pion)->Clone()));
		fxHNL[Flux::Parent::Pion]->Scale(heavy.Scale(Production::Channel::PionM));
	}

	//K+ -> mu+ nu_mu	(63.56%)
	//K+ -> pi0 mu+ nu_mu	(3.53%)
	if (fxNu.count(Flux::Parent::Kaon)
	&& (heavy.IsAllowed(Production::Channel::KaonM) || heavy.IsAllowed(Production::Channel::KaonCM))) {
		fxHNL[Flux::Parent::Kaon] = std::shared_ptr<TH1D>(static_cast<TH1D*>
					(fxNu.at(Flux::Parent::Kaon)->Clone()));
		double kf = (63.56 * heavy.Scale(Production::Channel::KaonM)
			    + 3.35 * heavy.Scale(Production::Channel::KaonCM))
			  / (63.56 + 3.35);
		fxHNL[Flux::Parent::Kaon]->Scale(kf);
	}

	//K0 -> pi- mu+ nu_mu
	if (fxNu.count(Flux::Parent::Kaon0)
	&& heavy.IsAllowed(Production::Channel::Kaon0M)) {
		fxHNL[Flux::Parent::Kaon0] = std::shared_ptr<TH1D>(static_cast<TH1D*>
					(fxNu.at(Flux::Parent::Kaon0)->Clone()));
		fxHNL[Flux::Parent::Kaon0]->Scale(heavy.Scale(Production::Channel::Kaon0M));
	}

	//mu- -> nu_mu e- nu_e_bar
	if (fxNu.count(Flux::Parent::Muon)
	&& heavy.IsAllowed(Production::Channel::MuonM)) {
		fxHNL[Flux::Parent::Muon] = std::shared_ptr<TH1D>(static_cast<TH1D*>
					(fxNu.at(Flux::Parent::Muon)->Clone()));
		fxHNL[Flux::Parent::Muon]->Scale(heavy.Scale(Production::Channel::MuonM));
	}

	//Ds+ -> mu+ nu_mu
	if (fxNu.count(Flux::Parent::Charm)
	&& heavy.IsAllowed(Production::Channel::CharmM)) {
		fxHNL[Flux::Parent::Charm] = std::shared_ptr<TH1D>(static_cast<TH1D*>
					(fxNu.at(Flux::Parent::Charm)->Clone()));
		fxHNL[Flux::Parent::Charm]->Scale(heavy.Scale(Production::Channel::CharmM));
	}

	if (!fxHNL.size())
		return nullptr;

	return std::accumulate(std::next(fxHNL.begin()), fxHNL.end(), fxHNL.begin()->second,
			[](std::shared_ptr<TH1D> total,
				const std::pair<Flux::Parent, std::shared_ptr<TH1D> > &hist) {
				total->Add(hist.second.get());
				return total; } );
}

//Make tauonic components, requires neutrino mass
//void Driver::MakeTauComponent(Flux *fxFlux, Neutrino &N)
std::shared_ptr<TH1D> Driver::MakeTau(ProductionRate &heavy, const Flux::Component &fxNu, double mass) const
{
	//auto fxHNL = fxNu;
	// start with empty map
	// only allowed processes will be added
	// if nothing is allowed a null pointer is returned
	Flux::Component fxHNL;

	//Ds+ -> tau+ nu_tau
	if (fxNu.count(Flux::Parent::Charm)
	&& heavy.IsAllowed(Production::Channel::CharmT)) {
		fxHNL[Flux::Parent::Charm] = std::shared_ptr<TH1D>(static_cast<TH1D*>
					(fxNu.at(Flux::Parent::Charm)->Clone()));
		double mul = 1.;
		if (_modifiers.count(Flux::Parent::Charm))
			mul = Stretch(fxHNL[Flux::Parent::Charm], _modifiers.at(Flux::Parent::Charm), mass);
		fxHNL[Flux::Parent::Charm]->Scale(heavy.Scale(Production::Channel::CharmT) * mul);
	}

	//tau+ -> pi+ nu_tau
	if (fxNu.count(Flux::Parent::Pion)
	&& heavy.IsAllowed(Production::Channel::TauPI)) {
		fxHNL[Flux::Parent::Pion] = std::shared_ptr<TH1D>(static_cast<TH1D*>
					(fxNu.at(Flux::Parent::Pion)->Clone()));
		double mul = 1.;
		if (_modifiers.count(Flux::Parent::Pion))
			mul = Stretch(fxHNL[Flux::Parent::Pion], _modifiers.at(Flux::Parent::Pion), mass);
		fxHNL[Flux::Parent::Pion]->Scale(heavy.Scale(Production::Channel::TauPI) * mul);
	}

	//tau+ -> pi+ pi0 nu_tau	//crossing simmetries
	if (fxNu.count(Flux::Parent::PPion)
	&& heavy.IsAllowed(Production::Channel::Tau2PI)) {
		fxHNL[Flux::Parent::PPion] = std::shared_ptr<TH1D>(static_cast<TH1D*>
					(fxNu.at(Flux::Parent::PPion)->Clone()));
		double mul = 1.;
		if (_modifiers.count(Flux::Parent::PPion))
			mul = Stretch(fxHNL[Flux::Parent::PPion], _modifiers.at(Flux::Parent::PPion), mass);
		fxHNL[Flux::Parent::PPion]->Scale(heavy.Scale(Production::Channel::Tau2PI) * mul);
	}

	//tau+ -> nu_tau_bar e+ nu_e
	if (fxNu.count(Flux::Parent::TauE)
	&& heavy.IsAllowed(Production::Channel::TauET)) {
		fxHNL[Flux::Parent::TauE] = std::shared_ptr<TH1D>(static_cast<TH1D*>
					(fxNu.at(Flux::Parent::TauE)->Clone()));
		double mul = 1.;
		if (_modifiers.count(Flux::Parent::TauE))
			mul = Stretch(fxHNL[Flux::Parent::TauE], _modifiers.at(Flux::Parent::TauE), mass);
		fxHNL[Flux::Parent::TauE]->Scale(heavy.Scale(Production::Channel::TauET) * mul);
	}

	//tau+ -> nu_tau_bar mu+ nu_mu
	if (fxNu.count(Flux::Parent::TauM)
	&& heavy.IsAllowed(Production::Channel::TauMT)) {
		fxHNL[Flux::Parent::TauM] = std::shared_ptr<TH1D>(static_cast<TH1D*>
					(fxNu.at(Flux::Parent::TauM)->Clone()));
		double mul = 1.;
		if (_modifiers.count(Flux::Parent::TauM))
			mul = Stretch(fxHNL[Flux::Parent::TauM], _modifiers.at(Flux::Parent::TauM), mass);
		fxHNL[Flux::Parent::TauM]->Scale(heavy.Scale(Production::Channel::TauMT) * mul);
	}

	if (!fxHNL.size())
		return nullptr;

	return std::accumulate(std::next(fxHNL.begin()), fxHNL.end(), fxHNL.begin()->second,
			[](std::shared_ptr<TH1D> total,
				const std::pair<Flux::Parent, std::shared_ptr<TH1D> > &hist) {
				total->Add(hist.second.get());
				return total; } );
}


//Modificator for charm to tau flux
double Driver::Stretch(std::shared_ptr<TH1D> hist, const Driver::Modifier &mod, double mass) const
{
	double sx = 0., ex = 0., vp = 0.;
	for (const auto &m : mod)
		if (m[0] > mass || std::abs(m[0] - mass) < 1.e-9) {
			sx = m[1]; ex = m[2]; vp = 1./m[3];
			break;
		}

	if (sx >= hist->GetXaxis()->GetXmin() && ex <= hist->GetXaxis()->GetXmax())
	{
		TH1D *htmp = static_cast<TH1D*>(hist->Clone());
		hist->Reset("ICES");

		int ia = htmp->FindFirstBinAbove();
		int ib = htmp->FindLastBinAbove();
		double eA = htmp->GetBinContent(ia);
		double eB = htmp->GetBinContent(ib);
		for (int i = ia; i <= ib; ++i) {
			double shift = sx + (htmp->GetBinContent(i) - eA) * (ex - sx);
			hist->SetBinContent(htmp->FindBin(shift),
					    htmp->GetBinContent(i) * (ex - sx) / (eB - eA));
		}
		return vp;
	}

	return 0.;
}
