#include "Driver.h"

Driver::Driver(const std::string &card)
{
	CardDealer cd(fluxConfig);

	// create production object for light neutrino
	Neutrino nu();
	_light = Production(nu);

	Init(cd);
}

// load fluxes for neutrinos and antineutrinos
void Init(const CardDealer &cd)
{
	std::map<std::string, std::string> > files;
	if (!cd.Get("flux_", files))
		throw std::invalid_argument("Driver: flux card does not contain any flux description\n");

	for (const auto &ff : files)
		_fxNu[Flavour::fromString(ff.first)] = std::move(Flux(ff.second));

	// load tau flux modifiers
	if (!fc.Get("mod_", files))
		return;

	for (const auto &id : files) {
		Modifier mod;
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
		_modifiers[Flux::fromString(id.first)] = std::move(modifier);
	}
}

//return the intensity of the flux at given energy
//
double Driver::Intensity(Neutrino &N, bool neut) const
{
	//double Energy = N.EnergyKin();
	double en = N.Momentum();	//as p ~ E for light neutrinos
	if (neut)
		for (auto ff = std::begin(Nu::Neut); ff != std::end(Nu::Neut); ++f)
			MakeComponent(*ff, fxHNL[*ff], N);
	else
		for (auto ff = std::begin(Nu::Anti); ff != std::end(Nu::Anti); ++f)
			MakeComponent(*ff, fxHNL[*ff], N);

	std::accumulate(std::begin(
	//std::cout << fxHeavyElectron << ", " << fxHeavyMuon << ", " << fxHeavyTau << std::endl;
	double Intensity = 0;
	if (fxHeavyElectron)
		Intensity += N.Ue(2) * InterpolateIntensity(fxHeavyElectron->Get(Flux:::Total), Energy);
	if (fxHeavyMuon)
		Intensity += N.Um(2) * InterpolateIntensity(fxHeavyMuon->Get(Flux::Total), Energy);
	if (fxHeavyTau)
		Intensity += N.Ut(2) * InterpolateIntensity(fxHeavyTau->Get(Flux::Total), Energy);

	return Intensity;
}

// prepare flux by scaling each component by the neutrino
// production rate. HNL flux is just each component times
// time the production rate
//
// opt can be Neutrino::particle or Neutrino::antiparticle
std::unordered_map<Nu::Flavour, std::shared_ptr<TH1D> > Driver::MakeFlux(const Neutrino &N)
{
	// production for HNL
	Production hnl(N);

	std::map<Nu::Flavour, std::shared_ptr<TH1D> > dist;
	if (N.IsParticle())
		for (const auto flv : Nu::Neut)
			dist[flv] = MakeComponent(hnl, flv, N.M());
	else if (N.IsAntiParticle())
		for (const auto flv : Nu::Anti)
			dist[flv] = MakeComponent(hnl, flv, N.M());

	return dist;
}

std::shared_ptr<TH1D> MakeComponent(Production &hnl, Nu::Flavour nu, double mass) {
	switch (nu) {
		case Nu::E_:
		case Nu::Eb:
			return MakeElectron(hnl, _fxNu[nu]);
		case Nu::M_:
		case Nu::Mb:
			return MakeMuon(hnl, _fxNu[nu]);
		case Nu::T_:
		case Nu::Tb: // only tau component requires neutrino mass
			return MakeTau(hnl, _fxNu[nu], mass);
		default:
			throw std::invalid_argument("Driver: unknwon neutrino flavour\n");
	}
}

//Make electronic components, requires Neutrino (T if neutrino, F if antineutrino), the Flux object, the mass
//void Driver::MakeElecComponent(Flux *fxFlux, Neutrino &N)
std::shared_ptr<TH1D> Driver::MakeElectron(Production &heavy, const Flux::Component &fxNu)
{
	auto fxHNL = fxNu;

	//pi+ -> e+ nu_e
	if (fxHNL.count(Flux::Parent::Pion))
		fxHNL[Flux::Parent::Pion]->Scale(heavy.Scale(Channel::PionE)
					      / _light.Scale(Channel::PionE));

	//K+ -> pi0 e+ nu_e	(5.07 %)	//I just assume that both decays are equally probable 
	//K+ -> e+ nu_e		(1.582e-3)	//for each different energy, so it is a linear combination
	if (fxHNL.count(Flux::Parent::Kaon)) {
		double kf = (1.582e-3 * heavy.Scale(Channel::KaonE) / _light.Scale(Channel::KaonE)
			       + 5.07 * heavy.Scale(Channel::KaonCE) / _light.Scale(Channel::KaonE))
			  / (1.582e-3 + 5.07);
		fxHNL[Flux::Parent::Kaon]->Scale(kf);
	}

	//K0 -> pi+ e+ nu_e
	if (fxHNL.count(Flux::Parent::Kaon0))
		fxHNL[Flux::Parent::Kaon0]->Scale(heavy.Scale(Channel::Kaon0E)
					       / _light.Scale(Channel::Kaon0E));

	//mu+ -> nu_mu_bar e+ nu_e
	if (fxHNL.count(Flux::Parent::Muon))
		fxHNL[Flux::Parent::Muon]->Scale(heavy.Scale(Channel::MuonE)
					      / _light.Scale(Channel::MuonE));

	//Ds+ -> e+ nu_e
	if (fxHNL.count(Flux::Parent::Charm))
		fxHNL[Flux::Parent::Charm]->Scale(heavy.Scale(Channel::CharmE)
					       / _light.Scale(Channel::CharmE));

	// combine into total
	return std::accumulate(fxHNL.begin()+1, fxHNL.end(), fxHNL.begin(),
			[](std::shared_ptr<TH1D> total, std::shared_ptr<TH1D> hist) {
				total->Add(hist);
				return total; } );
}

//Make muonic components, requires Neutrino (T if neutrino, F if antineutrino), the Flux object, the mass
//void Driver::MakeMuonComponent(Flux *fxFlux, Neutrino &N)
std::shared_ptr<TH1D> Driver::MakeMuon(Production &heavy, const Flux::Component &fxNu)
{
	auto fxHNL = fxNu;

	//pi+ -> mu+ nu_mu
	if (fxHNL.count(Flux::Parent::Pion))
		fxHNL[Flux::Parent::Pion]->Scale(heavy.Scale(Channel::PionM)
					      / _light.Scale(Channel::PionM));

	//K+ -> mu+ nu_mu	(63.56%)
	//K+ -> pi0 mu+ nu_mu	(3.53%)
	if (fxHNL.count(Flux::Parent::Kaon)) {
		double kf = (63.56 * heavy.Scale(Channel::KaonM) / _light.Scale(Channel::KaonM)
			    + 3.35 * heavy.Scale(Channel::KaonCM) / _light.Scale(Channel::KaonM))
			  / (63.56 + 3.35);
		fxHNL[Flux::Parent::Kaon]->Scale(kf);
	}

	//K0 -> pi- mu+ nu_mu
	if (fxHNL.count(Flux::Parent::Kaon0))
		fxHNL[Flux::Parent::Kaon0]->Scale(heavy.Scale(Channel::Kaon0M)
					       / _light.Scale(Channel::Kaon0M));

	//mu- -> nu_mu e- nu_e_bar
	if (fxHNL.count(Flux::Parent::Muon))
	fxHNL[Flux::Parent::Muon]->Scale(heavy.Scale(Channel::MuonM)
				      / _light.Scale(Channel::MuonM));

	//Ds+ -> mu+ nu_mu
	if (fxHNL.count(Flux::Parent::Charm))
		fxHNL[Flux::Parent::Charm]->Scale(heavy.Scale(Channel::CharmM)
					       / _light.Scale(Channel::CharmM));

	return std::accumulate(fxHNL.begin()+1, fxHNL.end(), fxHNL.begin(),
			[](std::shared_ptr<TH1D> total, std::shared_ptr<TH1D> hist) {
				total->Add(hist);
				return total; } );
}

//Make tauonic components, requires Neutrino (T if neutrino, F if antineutrino), the Flux object, the mass
//void Driver::MakeTauComponent(Flux *fxFlux, Neutrino &N)
std::shared_ptr<TH1D> Driver::MakeTau(Production &heavy, const Flux::Component &fxNu, double mass)
{
	auto fxHNL = fxNu;

	//Ds+ -> tau+ nu_tau
	if (fxHNL.count(Flux::Parent::Charm)) {
		double mul = 1.;
		if (modifiers.count(Flux::Parent::Charm))
			mul = Stretch(fxHNL[Flux::Parent::Charm], modifiers[Flux::Parent::Charm], mass);
		fxHNL[Flux::Parent::Charm]->Scale(heavy.Scale(Channel::CharmT) * mul
					       / _light.Scale(Channel::CharmT));
	}

	//tau+ -> pi+ nu_tau
	if (fxHNL.count(Flux::Parent::Pion)) {
		double mul = 1.;
		if (modifiers.count(Flux::Parent::Pion))
			mul = Stretch(fxHNL[Flux::Parent::Pion], modifiers[Flux::Parent::Pion], mass)
		fxHNL[Flux::Parent::Pion]->Scale(heavy.Scale(Channel::TauPI) * mul
					      / _light.Scale(Channel::TauPI));
	}

	//tau+ -> pi+ pi0 nu_tau	//crossing simmetries
	if (fxHNL.count(Flux::Parent::PPion)) {
		double mul = 1.;
		if (modifiers.count(Flux::Parent::PPion))
			mul = Stretch(fxHNL[Flux::Parent::PPion], modifiers[Flux::Parent::PPion], mass);
		fxHNL[Flux::Parent::PPion]->Scale(heavy.Scale(Channel::Tau2PI) * mul
					       / _light.Scale(Channel::Tau2PI));
	}

	//tau+ -> nu_tau_bar e+ nu_e
	if (fxHNL.count(Flux::Parent::TauE)) {
		double mul = 1.;
		if (modifiers.count(Flux::Parent::TauE))
			mul = Stretch(fxHNL[Flux::Parent::TauE], modifiers[Flux::Parent::TauE], mass);
		fxHNL[Flux::Parent::TauE]->Scale(heavy.Scale(Channel::TauET) * mul
					      / _light.Scale(Channel::TauET));
	}

	//tau+ -> nu_tau_bar mu+ nu_mu
	if (fxHNL.count(Flux::Parent::TauM)) {
		double mul = 1.;
		if (modifiers.count(Flux::Parent::TauM))
			mul = Stretch(fxHNL[Flux::Parent::TauM], modifiers[Flux::Parent::TauM], mass)
		fxHNL[Flux::Parent::TauM]->Scale(heavy.Scale(Channel::TauMT) * mul
					      / _light.Scale(Channel::TauMT));
	}

	fxHNL[Flux::Parent::Total]->Reset("ICES");
	for (const auto &ih : fxHNL)
		if (ih.first != Flux::Parent::Total)
			fxHNL[Flux::Parent::Total]->Add(ih.second.get());

	return std::accumulate(fxHNL.begin()+1, fxHNL.end(), fxHNL.begin(),
			[](std::shared_ptr<TH1D> total, std::shared_ptr<TH1D> hist) {
				total->Add(hist);
				return total; } );
}


//Modificator for charm to tau flux
double Driver::Stretch(std::shared_ptr<TH1D> hist, const Modifier &mod, double mass)
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
			double shif = sx + (htpm->GetBinContent(i) - eA) * (ex - sx);
			hist->SetBinContent(htmp->FindBin(shif),
					    htmp->GetBinContent(i) * (ex - sx) / (eB - eA));
		}
		return 1.;
	}

	return 0.;
}
