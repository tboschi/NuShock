#include "Flux.h"

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
	CloneCopy(hTotal, f.GetTotal());
	CloneCopy(hPion,  f.GetPion());
	CloneCopy(h2Pion, f.Get2Pion());
	CloneCopy(hKaon,  f.GetKaon());
	CloneCopy(hKaon0, f.GetKaon0());
	CloneCopy(hCharm, f.GetCharm());
	CloneCopy(hMuon,  f.GetMuon());
	CloneCopy(hTauE,  f.GetTauE());
	CloneCopy(hTauM,  f.GetTauM());
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

TH1D* Flux::GetTotal() const
{
	return hTotal;
}

TH1D* Flux::GetPion() const
{
	return hPion;
}

TH1D* Flux::Get2Pion() const
{
	return h2Pion;
}

TH1D* Flux::GetKaon() const
{
	return hKaon;
}

TH1D* Flux::GetKaon0() const
{
	return hKaon0;
}

TH1D* Flux::GetCharm() const
{
	return hCharm;
}

TH1D* Flux::GetMuon() const
{
	return hMuon;
}

TH1D* Flux::GetTauE() const
{
	return hTauE;
}

TH1D* Flux::GetTauM() const
{
	return hTauM;
}
