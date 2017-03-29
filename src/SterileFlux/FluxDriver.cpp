#include <FluxDriver.h>

FluxDriver::FluxDriver(std::string SourceName)
{
	SourceFile = new TFile(SourceName.c_str(), "OPEN");

	hTotalFlux = (TH1D*) SourceFile->Get("htotalflux");
	hMuonPion = (TH1D*) SourceFile->Get("hmuonpiona");
	hMuonKaon = (TH1D*) SourceFile->Get("hmuonkaon");
	hElectronPion = (TH1D*) SourceFile->Get("helectronpion");
	hElectronKaon = (TH1D*) SourceFile->Get("helectronkaon");
	hElectronKaon3 = (TH1D*) SourceFile->Get("helectronkaon3");
	hMuonKaonOther = (TH1D*) SourceFile->Get("hmuonkaonother");
}

FluxDriver::~FluxDriver()
{
	SourceFile->Close();
}

TH1D *GetHist()
{
	return hTotalFlux;
}

void FluxDriver::MakeSterileFlux(double M_Sterile)
{
	//Clone from original fluxes
	sMuonPion = (TH1D*) hMuonPion->Clone();
	sMuonKaon = (TH1D*) hMuonKaon->Clone();
	sElectronPion = (TH1D*) hElectronPion->Clone();
	sElectronKaon = (TH1D*) hElectronKaon->Clone();
	sElectronKaon3 = (TH1D*) hElectronKaon3->Clone();
	sMuonKaonOther = (TH1D*) hMuonKaonOther->Clone();

	//Scale accordingly
	sMuonPion->Scale(Um*Um * ShrockFactor(M_Pion, M_Muon, M_Sterile));
	sMuonKaon->Scale(Um*Um * ShrockFactor(M_Kaon, M_Muon, M_Sterile));
	sElectronPion->Scale(Ue*Ue * ShrockFactor(M_Pion, M_Electron, M_Sterile));
	sElectronKaon->Scale(Ue*Ue * ShrockFactor(M_Kaon, M_Electron, M_Sterile));
	sElectronKaon3->Scale(Ue*Ue);
	sMuonKaonOther->Scale(Um*Um * ShrockFactor(M_Kaon, M_Muon, M_Sterile));

	//Add fluxes to total
	sTotalFlux->Add(sMuonPion);
	sTotalFlux->Add(sMuonKaon);
	sTotalFlux->Add(sElectronPion);
	sTotalFlux->Add(sElectronKaon);
	sTotalFlux->Add(sElectronKaon3);
	sTotalFlux->Add(sMuonKaonOther);
}

//Sample energy and sets values to given pointers. Detector smearing must me implemented!!
double FluxDriver::SampleEnergy(Flux *StdFlux, Flux *HeavyFlux)
{
	double RanEnergy = sTotalFlux->GetRandom();

	StdFlux->SetEnergy(RanEnergy);
	StdFlux->SetMuonPion(hMuonPion->GetBinContent(hMuonPion->FindBin(RanEnergy)));
	StdFlux->SetMuonKaon(hMuonKaon->GetBinContent(hMuonKaon->FindBin(RanEnergy)));
	StdFlux->SetElectronPion(hElectronPion->GetBinContent(hElectronPion->FindBin(RanEnergy)));
	StdFlux->SetElectronKaon(hElectronKaon->GetBinContent(hElectronKaon->FindBin(RanEnergy)));
	StdFlux->SetElectronKaon3(hElectronKaon3->GetBinContent(hElectronKaon3->FindBin(RanEnergy)));
	StdFlux->SetMuonKaonOther(hMuonKaonOther->GetBinContent(hMuonKaonOther->FindBin(RanEnergy)));

	HeavyFlux->SetEnergy(RanEnergy);
	HeavyFlux->SetMuonPion(sMuonPion->GetBinContent(sMuonPion->FindBin(RanEnergy)));
	HeavyFlux->SetMuonKaon(sMuonKaon->GetBinContent(sMuonKaon->FindBin(RanEnergy)));
	HeavyFlux->SetElectronPion(sElectronPion->GetBinContent(sElectronPion->FindBin(RanEnergy)));
	HeavyFlux->SetElectronKaon(sElectronKaon->GetBinContent(sElectronKaon->FindBin(RanEnergy)));
	HeavyFlux->SetElectronKaon3(sElectronKaon3->GetBinContent(sElectronKaon3->FindBin(RanEnergy)));
	HeavyFlux->SetMuonKaonOther(sMuonKaonOther->GetBinContent(sMuonKaonOther->FindBin(RanEnergy)));

	return RanEnergy;
}
