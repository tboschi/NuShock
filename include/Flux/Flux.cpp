#include "Flux/Flux.h"

//ctor
Flux::Flux(std::string HistFile) :
	hTotal(0),
        hPion(0),
        h2Pion(0),
        hKaon(0),
        hKaon0(0),
        hCharm(0),
        hMuon(0),
        hTauE(0),
        hTauM(0)
{
	TFile* InFile = new TFile(HistFile.c_str(), "READ");

	CloneCopy(hTotal, InFile->Get("htotal"));
	CloneCopy(hPion,  InFile->Get("hpion"));
	CloneCopy(h2Pion, InFile->Get("h2pion"));
	CloneCopy(hKaon,  InFile->Get("hkaon"));
	CloneCopy(hKaon0, InFile->Get("hkaon0"));
	CloneCopy(hCharm, InFile->Get("hcharm"));
	CloneCopy(hMuon,  InFile->Get("hmuon"));
	CloneCopy(hTauE,  InFile->Get("htaue"));
	CloneCopy(hTauM,  InFile->Get("htaum"));

	InFile->Close();

}

//copy ctor
Flux::Flux(const Flux & f)
{
	CloneCopy(hTotal, f.Get(Total));
	CloneCopy(hPion,  f.Get(Pion));
	CloneCopy(h2Pion, f.Get(2Pion));
	CloneCopy(hKaon,  f.Get(Kaon));
	CloneCopy(hKaon0, f.Get(Kaon0));
	CloneCopy(hCharm, f.Get(Charm));
	CloneCopy(hMuon,  f.Get(Muon));
	CloneCopy(hTauE,  f.Get(TauE));
	CloneCopy(hTauM,  f.Get(TauM));
}

//detor
Flux::~Flux()
{
	delete hTotal;
	delete hPion;
	delete h2Pion;
	delete hKaon;
	delete hKaon0;
	delete hCharm;
	delete hMuon;
	delete hTauE;
	delete hTauM;
}

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


//Get functions

TH1D* Flux::Get(Hist KeyName) const
{
	switch (KeyName)
	{
		case Total:
			return hTotal;
		case Pion:
			return hPion;
		case 2Pion:
			return h2Pion;
		case Kaon:
			return hKaon;
		case Kaon0:
			return hKaon0;
		case Charm:
			return hCharm;
		case Muon:
			return hMuon;
		case TauE:
			return hTauE;
		case TauM:
			return hTauM;
		default:
			return NULL;
	}
}

void Flux::Scale(double X, Hist KeyName)
{
	switch (KeyName)
	{
		case Total:
			Scale(X, Hist::Pion);
			Scale(X, Hist::2Pion);
			Scale(X, Hist::Kaon);
			Scale(X, Hist::Kaon0);
			Scale(X, Hist::Charm);
			Scale(X, Hist::Muon);
			Scale(X, Hist::TauE);
			Scale(X, Hist::TauM);
			break;
		case Pion:
			hPion->Scale(X);
			break;
		case 2Pion:
			h2Pion->Scale(X);
			break;
		case Kaon:
			hKaon->Scale(X);
			break;
		case Kaon0:
			hKaon0->Scale(X);
			break;
		case Muon:
			hMuon->Scale(X);
			break;
		case TauE:
			hTauE->Scale(X);
			break;
		case TauM:
			hTauM->Scale(X);
			break;
		deafult:
			break;
	}
}

double Flux::RangeStart()
{
	return Get(Hist::Total)->GetBinCenter(0);
}

double Flux::RangeEnd()
{
	return Get(Hist::Total)->GetBinCenter(BinNumber());
}

double Flux::BinNumber()
{
	return Get(Hist::Total)->GetNbinsX(0);
}

double Flux::BinWidth()
{
	return (RangeEnd() - RangeStart()) / BinNumber();
}
