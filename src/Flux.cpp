#include "Flux.h"

//ctor
Flux::Flux(std::string HistFile) :
	hTotal(0),
        hPion(0),
        hPPion(0),
        hKaon(0),
        hKaon0(0),
        hCharm(0),
        hMuon(0),
        hTauE(0),
        hTauM(0)
{
	TFile* InFile = new TFile(HistFile.c_str(), "READ");
	if (InFile->IsZombie())
		std::cout << "File " << HistFile << " does not exist" << std::endl;

	CloneCopy(hTotal, InFile->Get("htotal"));
	CloneCopy(hPion,  InFile->Get("hpion"));
	CloneCopy(hPPion, InFile->Get("h2pion"));
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
	CloneCopy(hPPion, f.Get(PPion));
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
	delete hPPion;
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

TH1D* Flux::Get(Hist Name) const
{
	switch (Name)
	{
		case Total:
			return hTotal;
		case Pion:
			return hPion;
		case PPion:
			return hPPion;
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

void Flux::Add()
{
	Get(Total)->Reset("ICES");

	Add(Pion);
	Add(PPion);
	Add(Kaon);
	Add(Kaon0);
	Add(Charm);
	Add(Muon);
	Add(TauE);
	Add(TauM);
}

void Flux::Add(Hist Name)
{
	TH1D* hComponent = Get(Name);

	if (hComponent)
		Get(Total)->Add(hComponent);
}

void Flux::Scale(double X)
{
	Scale(Pion, X);
	Scale(PPion, X);
	Scale(Kaon, X);
	Scale(Kaon0, X);
	Scale(Charm, X);
	Scale(Muon, X);
	Scale(TauE, X);
	Scale(TauM, X);
}

void Flux::Scale(Hist Name, double X)
{
	TH1D* hComponent = Get(Name);

	if (hComponent)
		hComponent->Scale(X);
}

bool Flux::Stretch(Hist Name, double Sx, double Ex)
{
	TH1D* hComponent = Get(Name);

	if (hComponent && Sx >= RangeStart() && Ex <= RangeEnd())
	{
		TH1D *hTemp = dynamic_cast<TH1D*> (hComponent->Clone());
		hComponent->Reset("ICES");

		double A = hTemp->GetXaxis()->GetBinCenter(hTemp->FindFirstBinAbove(0));
		double B = hTemp->GetXaxis()->GetBinCenter(hTemp->FindLastBinAbove(0));

		double m = (Ex - Sx) / (B - A);
		double Start = RangeStart();
		double End   = RangeEnd();
		double EnStep = (End - Start) / BinNumber();
		for (double Energy = Start; Energy < End; Energy += EnStep)
		{
			double Flux = hTemp->GetBinContent(hTemp->FindBin(Energy));
			hComponent->Fill( Sx + (Energy - A) * (Ex - Sx) / (B - A), Flux * (Ex - Sx) / (B - A) );
		}

		delete hTemp;
		return true;
	}
	else
		return false;
}

double Flux::RangeStart()
{
	return Get(Total)->GetBinLowEdge(1);
}

double Flux::RangeEnd()
{
	return Get(Total)->GetBinLowEdge(BinNumber() + 1);
}

double Flux::BinNumber()
{
	return Get(Total)->GetNbinsX();
}

double Flux::BinWidth()
{
	return Get(Total)->GetBinWidth(1);
}
