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

	SterileFlux = new Flux(0);

	M_Pion = genie::PDGLibrary::Instance()->Find(211)->Mass();
	M_Kaon = genie::PDGLibrary::Instance()->Find(321)->Mass();
	M_Muon = genie::PDGLibrary::Instance()->Find(13)->Mass();
	M_Electron = genie::PDGLibrary::Instance()->Find(11)->Mass();
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

//Sample energy and return a flux. Detector smearing must me implemented!!
Flux *FluxDriver::SampleEnergy()
{
	double RanEnergy = sTotalFlux->GetRandom();
	double fMuonPion = sMuonPion->GetBinContent(sMuonPion->FindBin(RanEnergy));
        double fMuonKaon = sMuonKaon->GetBinContent(sMuonKaon->FindBin(RanEnergy));
        double fElectronPion = sElectronPion->GetBinContent(sElectronPion->FindBin(RanEnergy));
        double fElectronKaon = sElectronKaon->GetBinContent(sElectronKaon->FindBin(RanEnergy));
        double fElectronKaon3 = sElectronKaon3->GetBinContent(sElectronKaon3->FindBin(RanEnergy));
        double fMuonKaonOther = sMuonKaonOther->GetBinContent(sMuonKaonOther->FindBin(RanEnergy));

	//Return a pointer to Flux object with components
	SterileFlux->SetAll(RanEnergy, fMuonPion, fMuonKaon, fElectronPion, fElectronKaon, fElectronKaon3, fMuonKaonOther);
	return SterileFlux;
}

//Kinematic factors
double FluxDriver::ShrockFactor(double M_Meson, double M_Lepton, double M_Sterile) //Correct scaling of flux
{
	double dM_a = pow(M_Lepton/M_Meson, 2.0);
	double dM_i = pow(M_Sterile/M_Meson, 2.0);	

	if (M_Meson >= M_Lepton + M_Sterile)
		return ShrockRho(dM_a, dM_i)/(dM_a * pow(1-dM_a, 2.0));
	else return 0;
}

double FluxDriver::ShrockRho(double X double Y)
{
	return ShrockFM(X, Y)*sqrt(ShrockLambda(1, X, Y));
}

double FluxDriver::ShrockFM(double X, double Y)
{
	return X+Y - (X-Y)*(X-Y);
}

double FluxDriver::ShrockLambda(double X, double Y, double Z)
{
	return X*X + Y*Y + Z*Z - 2*(X*Y + X*Z + Y*Z);
}
