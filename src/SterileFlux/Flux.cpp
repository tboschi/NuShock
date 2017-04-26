#include "Flux.h"

//ctor
Flux::Flux(TH1D* Total, TH1D* Pion, TH1D* Kaon, TH1D* Kaon0, TH1D* Muon)
{
	CloneTotal(Total);
	ClonePion(Pion);
	CloneKaon(Kaon);
	CloneKaon0(Kaon0);
	CloneMuon(Muon);
}

Flux::Flux(std::string HistFile)
{
	TFile* InFile = new TFile(HistFile.c_str(), "READ");

	CloneTotal((TH1D*) InFile->Get("htotal"));
	ClonePion((TH1D*) InFile->Get("hpion"));
	CloneKaon((TH1D*) InFile->Get("hkaon"));
	CloneKaon0((TH1D*) InFile->Get("hkaon0"));
	CloneMuon((TH1D*) InFile->Get("hmuon"));

	InFile->Close();
}

//copy ctor
Flux::Flux(const Flux & f)
{
	CloneTotal(f.GetTotal());
	ClonePion(f.GetPion());
	CloneKaon(f.GetKaon());
	CloneKaon0(f.GetKaon0());
	CloneMuon(f.GetMuon());
}

//dector
/*
Flux::~Flux()
{
	delete hTotal;
	delete hPion;
	delete hKaon;
	delete hKaon0;
	delete hMuon;
}*/

//Clone functions, so that an object from this class owns valid copies of the histograms
void Flux::CloneTotal(TH1D* X)
{
	hTotal = (TH1D*) X->Clone();
	hTotal->SetDirectory(0);
}

void Flux::ClonePion(TH1D* X)
{
	hPion = (TH1D*) X->Clone();
	hPion->SetDirectory(0);
}

void Flux::CloneKaon(TH1D* X)
{
	hKaon = (TH1D*) X->Clone();
	hKaon->SetDirectory(0);
}

void Flux::CloneKaon0(TH1D* X)
{
	hKaon0 = (TH1D*) X->Clone();
	hKaon0->SetDirectory(0);
}

void Flux::CloneMuon(TH1D* X)
{
	hMuon = (TH1D*) X->Clone();
	hMuon->SetDirectory(0);
}

//Get functions

TH1D* Flux::GetTotal() const
{
	return hTotal;
}

TH1D* Flux::GetPion() const
{
	return hPion;
}

TH1D* Flux::GetKaon() const
{
	return hKaon;
}

TH1D* Flux::GetKaon0() const
{
	return hKaon0;
}

TH1D* Flux::GetMuon() const
{
	return hMuon;
}
