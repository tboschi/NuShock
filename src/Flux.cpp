#include "Flux.h"

//ctor
Flux::Flux(const std::string &histFile) :
{
	TFile infile(histFile.c_str(), "READ");
	if (infile.IsZombie())
		throw std::invalid_argument("Flux: file does not exist\n");


	TIter next(infile.GetListOfKeys());
	TKey* k;
	while (k = static_cast<TKey*>(next())) {
		std::string name = k->GetName();
		std::shared_ptr<TH1D> hflux(static_cast<TH1D*>(infile.Get(name.c_str())));
		hflux->SetDirectory(0);
		_hists[name] = std::move(hflux);
	}

	if (!_hists.size())
		throw std::logic_error("No histogram collected! Very bad\n")

	if (!_hists.count("htotal")) {
		_hists["htotal"] = std::shared_ptr<TH1D>(*_hists.begin());
		Combine();
	}
}

	//CloneCopy(hTotal, InFile->Get("htotal"));
	//CloneCopy(hPion,  InFile->Get("hpion"));
	//CloneCopy(hPPion, InFile->Get("h2pion"));
	//CloneCopy(hKaon,  InFile->Get("hkaon"));
	//CloneCopy(hKaon0, InFile->Get("hkaon0"));
	//CloneCopy(hCharm, InFile->Get("hcharm"));
	//CloneCopy(hMuon,  InFile->Get("hmuon"));
	//CloneCopy(hTauE,  InFile->Get("htaue"));
	//CloneCopy(hTauM,  InFile->Get("htaum"));


/*
//Clone functions, so that an object from this class owns valid copies of the histograms
void Flux::CloneCopy(TH1D*& T, TObject* X)
{
	if (X)
	{
		T = dynamic_cast<TH1D*> (X->Clone());
		T->SetDirectory(0);
	}
	else
		T = NULL;
}

void Flux::CloneCopy(TH1D*& T, TH1D* X)
{
	if (X)
	{
		T = dynamic_cast<TH1D*> (X->Clone());
		T->SetDirectory(0);
	}
	else
		T = NULL;
}
*/


//Get functions


std::shared_ptr<TH1D> Flux::Get(const std::string &name) const
{
	if (_hists.count(name))
		return _hists[name];
	return nullptr;
}

std::shared_ptr<TH1D> Flux::Get(const Component &name) const
{
	switch (name)
	{
		case Total:
			return _hists["htotal"];
		case Pion:
			return _hists["hpion"];
		case PPion:
			return _hists["hppion"];
		case Kaon:
			return _hists["hkaon"];
		case Kaon0:
			return _hists["hkaon0"];
		case Charm:
			return _hists["hcharm"];
		case Muon:
			return _hists["hmuon"];
		case TauE:
			return _hists["htaue"];
		case TauM:
			return _hists["htaum"];
		default:
			return nullptr;
	}
}

void Flux::Combine()
{
	_hists["htotal"]->Reset("ICES");
	for (const auto &ih : _hists)
		_hists["htotal"]->Add(ih.second.get());
}

void Flux::Scale(double X)
{
	for (const auto &ih : _hists)
		ih.second->Scale(x);
}

void Flux::Scale(const Component &name, double x)
{
	if (_hists.count(name))
		_hists[name]->Scale(x);
}

// stretch the flux such that it has starting and ending energies at sx and ex
bool Flux::Stretch(const Component &name, double sx, double ex)
{
	if (!_hists.count(name))
		return false;

	if (sx >= RangeStart() && ex <= RangeEnd())
	{
		TH1D *htmp = static_cast<TH1D*>(_hists[name]->Clone());
		htmp->Reset("ICES");

		int ia = htmp->FindFirstBinAbove();
		int ib = htmp->FindLastBinAbove();
		double eA = htmp->GetBinContent(ia);
		double eB = htmp->GetBinContent(ib);
		for (int i = ia; i <= ib; ++i) {
			double shifted = sx + (htpm->GetBinContent(i) - eA) * (ex - sx);
			_hists[name]->SetBinContent(htmp->FindBin(shifted),
						    htmp->GetBinContent(i) * (ex - sx) / (eB - eA));
		}

		delete htmp;
		return true;
	}

	return false;
}

double Flux::RangeStart()
{
	return _hists["htotal"]->GetBinLowEdge(1);
}

double Flux::RangeEnd()
{
	return _hists["htotal"]->GetBinLowEdge(BinNumber() + 1);
}

double Flux::BinNumber()
{
	return _hists["htotal"]->GetNbinsX();
}

double Flux::BinWidth()
{
	return _hists["htotal"]->GetBinWidth(1);
}
