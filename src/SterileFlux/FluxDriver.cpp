#include <FluxDriver.h>

FluxDriver::FluxDriver(std::string ConfigFlux)
{
	fxNuMuon = 0;
	fxNuMuonBar = 0;
	fxNuElectron = 0;
	fxNuElectronBar = 0;

	std::string Line, Key, Element;
	std::stringstream ssL;
	
	ConfigFile.open(FluxConfig.c_str())
	while (std::getline(ConfigFile, Line))
	{
		if (Line[0] == '#') continue;

		ssL.str("");
		ssL.clear();
		ssL << Line;

		ssL >> Key >> Name;
		if (Key.find("Muon_") != std::string::npos) fxNuMuon = new Flux(Name);
		if (Key.find("MuonBar_") != std::string::npos) fxNuMuonBar = new Flux(Name);
		if (Key.find("Electron_") != std::string::npos) fxNuElectron = new Flux(Name);
		if (Key.find("ElectronBar_") != std::string::npos) fxNuElectronBar = new Flux(Name);
	}
	ConfigFile.close();
}

FluxDriver::~FluxDriver()
{
	SourceFile->Close();
}

TH1D *FluxDriver::GetHist()
{
	return sTotalFlux;
}

void FluxDriver::MakeSterileFlux(double M_Sterile, double M_Lepton, double U_x, double U_l)
{
	//Clone from original fluxes
	sPion = (TH1D*) hPion->Clone();
	sKaon = (TH1D*) hKaon->Clone();
	sKaon0 = (TH1D*) hKaon0->Clone();
	sMuon = (TH1D*) hMuon->Clone();

	//Scale accordingly
	sPion->Scale(U_x*U_x * Kine::ShrockFactor(M_Pion, M_Lepton, M_Sterile));
	sKaon->Scale(U_x*U_x * Kine::ShrockFactor(M_Kaon, M_Lepton, M_Sterile));
	sKaon0->Scale(U_x*U_x);
	sMuon->Scale(U_x*U_x);		//Not correct

	//Add fluxes to total
	sTotalFlux->Add(sPion);
	sTotalFlux->Add(sKaon);
	sTotalFlux->Add(sKaon0);
	sTotalFlux->Add(sMuon);
}

//Need to think about this
//void FluxDriver::DetectorModel(TH1D * hHist, double Probabilty)
//{
//	for (int i = 1; i < sPion->GetNbinsX()+1; ++i)
//		hHist->SetBinContent(i, sPion->GetBinContent(i)*Probability);
//}

//Sample energy and sets values to given pointers.
/*
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
*/
