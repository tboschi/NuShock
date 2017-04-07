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
	delete fxNuMuon;
	delete fxNuMuonBar;
	delete fxNuElectron;
	delete fxNuElectronBar;
}

void FluxDriver::MakeSterileFlux(double M_Sterile, double U_e, double U_m, double U_t)
{
	hTotalSterile->Reset("ICES");	//Reset before generating
	hPionlSterile->Reset("ICES");
	hKaonSterile->Reset("ICES"); 
	hKaon0Sterile->Reset("ICES");
	hMuonSterile->Reset("ICES"); 

	if (!fxNuMuon)		//Model on nu mu flux
	{
		Flux sxFlux = *fxNuMuon;

		sxFlux.GetPion()->Scale(U_m*U_m * Kine::ShrockFactor(M_Pion, M_Muon, M_Sterile));
		hPionSterile->Add(sxFlux.GetPion());

		sxFlux.GetKaon()->Scale(U_m*U_m * Kine::ShrockFactor(M_Kaon, M_Muon, M_Sterile));
		hKaonSterile->Add(sxFlux.GetKaon());

		sxFlux.GetKaon0()->Scale(U_m*U_m);
		hKaon0Sterile->Add(sxFlux.GetKaon0());

		sxFlux.GetMuon()->Scale(U_m*U_m);
		hMuonSterile->Add(sxFlux.GetMuon());
	}

	if (!fxNuMuonBar)	//Model on nu mu bar flux
	{
		Flux sxFlux = *fxNuMuonBar;

		sxFlux.GetPion()->Scale(U_m*U_m * Kine::ShrockFactor(M_Pion, M_Muon, M_Sterile));
		hPionSterile->Add(sxFlux.GetPion());

		sxFlux.GetKaon()->Scale(U_m*U_m * Kine::ShrockFactor(M_Kaon, M_Muon, M_Sterile));
		hKaonSterile->Add(sxFlux.GetKaon());

		sxFlux.GetKaon0()->Scale(U_m*U_m);
		hKaon0Sterile->Add(sxFlux.GetKaon0());

		sxFlux.GetMuon()->Scale(U_m*U_m);
		hMuonSterile->Add(sxFlux.GetMuon());
	}
	
	if (!fxNuElectron)	//Model on nu electron flux
	{
		Flux sxFlux = *fxNuElectron;

		sxFlux.GetPion()->Scale(U_e*U_e * Kine::ShrockFactor(M_Pion, M_Muon, M_Sterile));
		hPionSterile->Add(sxFlux.GetPion());

		sxFlux.GetKaon()->Scale(U_e*U_e * Kine::ShrockFactor(M_Kaon, M_Muon, M_Sterile));
		hKaonSterile->Add(sxFlux.GetKaon());

		sxFlux.GetKaon0()->Scale(U_e*U_e);
		hKaon0Sterile->Add(sxFlux.GetKaon0());

		sxFlux.GetMuon()->Scale(U_e*U_e);
		hMuonSterile->Add(sxFlux.GetMuon());
	}

	if (!fxNuElectronBar)	//Model on nu electron bar flux
	{
		Flux sxFlux = *fxNuElectronBar;

		sxFlux.GetPion()->Scale(U_e*U_e * Kine::ShrockFactor(M_Pion, M_Muon, M_Sterile));
		hPionSterile->Add(sxFlux.GetPion());

		sxFlux.GetKaon()->Scale(U_e*U_e * Kine::ShrockFactor(M_Kaon, M_Muon, M_Sterile));
		hKaonSterile->Add(sxFlux.GetKaon());

		sxFlux.GetKaon0()->Scale(U_e*U_e);
		hKaon0Sterile->Add(sxFlux.GetKaon0());

		sxFlux.GetMuon()->Scale(U_e*U_e);
		hMuonSterile->Add(sxFlux.GetMuon());
	}

	hTotalSterile->Add(hPionSterile);
	hTotalSterile->Add(hKaonSterile);
	hTotalSterile->Add(hKaon0Sterile);
	hTotalSterile->Add(hMuonSterile);
}

double FluxDriver::SampleEnergy()		//Sample energy
{
	return hTotalSterile->GetRandom();
} 

//Get individual histograms

TH1D* FluxDriver::GetSterile()
{
	return hTotalSterile;
} 

TH1D* FluxDriver::GetPion()
{
	return hPionSterile;
} 

TH1D* FluxDriver::GetKaon()
{
	return hKaonSterile;
} 

TH1D* FluxDriver::GetKaon0()
{
	return hKaon0Sterile;
} 

TH1D* FluxDriver::GetMuon()
{
	return hMuonSterile;
} 
