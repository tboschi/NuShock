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
	CloneCopy(hTotal, f.Get(Hist::Total));
	CloneCopy(hPion,  f.Get(Hist::Pion));
	CloneCopy(h2Pion, f.Get(Hist::2Pion));
	CloneCopy(hKaon,  f.Get(Hist::Kaon));
	CloneCopy(hKaon0, f.Get(Hist::Kaon0));
	CloneCopy(hCharm, f.Get(Hist::Charm));
	CloneCopy(hMuon,  f.Get(Hist::Muon));
	CloneCopy(hTauE,  f.Get(Hist::TauE));
	CloneCopy(hTauM,  f.Get(Hist::TauM));
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
		case Hist::Total:
			return hTotal;
		case Hist::Pion:
			return hPion;
		case Hist::2Pion:
			return h2Pion;
		case Hist::Kaon:
			return hKaon;
		case Hist::Kaon0:
			return hKaon0;
		case Hist::Charm:
			return hCharm;
		case Hist::Muon:
			return hMuon;
		case Hist::TauE:
			return hTauE;
		case Hist::TauM:
			return hTauM;
		default:
			return NULL;
	}
}
