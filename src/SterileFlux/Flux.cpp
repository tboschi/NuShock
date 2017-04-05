#include "Flux.h"

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
	TFile* InFile = new TFile(HistFile, "READ");

	CloneTotal((TH1D*) InFile->Get("htotal"));
	CloneTotal((TH1D*) InFile->Get("hpion"));
	CloneTotal((TH1D*) InFile->Get("hkaon"));
	CloneTotal((TH1D*) InFile->Get("hkaon0"));
	CloneTotal((TH1D*) InFile->Get("hmuon"));

	InFile->Close();
}

//Clone functions

void Flux::CloneTotal(TH1D* X)	//Rescale each component
{
	hTotal = (TH1D*) X->Clone();
}

void Flux::ClonePion(TH1D* X);
{
	hPion = (TH1D*) X->Clone();
}

void Flux::CloneKaon(TH1D* X);
{
	hKaon = (TH1D*) X->Clone();
}

void Flux::CloneKaon0(TH1D* X);
{
	hKaon0 = (TH1D*) X->Clone();
}

void Flux::CloneMuon(TH1D* X);
{
	hMuon = (TH1D*) X->Clone();
}

//Get functions

TH1D* Flux::GetTotal()
{
	return hTotal;
}

TH1D* Flux::GetPion()
{
	return hPion;
}

TH1D* Flux::GetKaon()
{
	return hKaon;
}

TH1D* Flux::GetKaon0()
{
	return hKaon0;
}

TH1D* Flux::GetMuon()
{
	return hMuon;
}
