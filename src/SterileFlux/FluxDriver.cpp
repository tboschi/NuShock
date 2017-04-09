#include <FluxDriver.h>

FluxDriver::FluxDriver(std::string FluxConfig) : 
	M_Electron(Const::fMElectron),
	M_Muon(Const::fMMuon),
	M_Pion(Const::fMPion),
	M_Kaon(Const::fMKaon)
{
	fxNuMuon = 0;
	fxNuMuonBar = 0;
	fxNuElectron = 0;
	fxNuElectronBar = 0;

	std::string Line, Key, Name;
	std::stringstream ssL;
	
	std::ifstream ConfigFile(FluxConfig.c_str());
	while (std::getline(ConfigFile, Line))
	{
		if (Line[0] == '#') continue;

		ssL.str("");
		ssL.clear();
		ssL << Line;

		if (ssL >> Key >> Name)
		{
			if (Key.find("Muon_") != std::string::npos) fxNuMuon = new Flux(Name);
			if (Key.find("MuonBar_") != std::string::npos) fxNuMuonBar = new Flux(Name);
			if (Key.find("Electron_") != std::string::npos) fxNuElectron = new Flux(Name);
			if (Key.find("ElectronBar_") != std::string::npos) fxNuElectronBar = new Flux(Name);
		}
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
	hPionSterile->Reset("ICES");
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

		sxFlux.GetPion()->Scale(U_e*U_e * Kine::ShrockFactor(M_Pion, M_Electron, M_Sterile));
		hPionSterile->Add(sxFlux.GetPion());

		sxFlux.GetKaon()->Scale(U_e*U_e * Kine::ShrockFactor(M_Kaon, M_Electron, M_Sterile));
		hKaonSterile->Add(sxFlux.GetKaon());

		sxFlux.GetKaon0()->Scale(U_e*U_e);
		hKaon0Sterile->Add(sxFlux.GetKaon0());

		sxFlux.GetMuon()->Scale(U_e*U_e);
		hMuonSterile->Add(sxFlux.GetMuon());
	}

	if (!fxNuElectronBar)	//Model on nu electron bar flux
	{
		Flux sxFlux = *fxNuElectronBar;

		sxFlux.GetPion()->Scale(U_e*U_e * Kine::ShrockFactor(M_Pion, M_Electron, M_Sterile));
		hPionSterile->Add(sxFlux.GetPion());

		sxFlux.GetKaon()->Scale(U_e*U_e * Kine::ShrockFactor(M_Kaon, M_Electron, M_Sterile));
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

void FluxDriver::SetBaseline(double Baseline)
{
	hTotalSterile->Scale(1e6/(Baseline*Baseline));
	hPionSterile->Scale(1e6/(Baseline*Baseline));
	hKaonSterile->Scale(1e6/(Baseline*Baseline));
	hKaon0Sterile->Scale(1e6/(Baseline*Baseline));
	hPionSterile->Scale(1e6/(Baseline*Baseline));
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
