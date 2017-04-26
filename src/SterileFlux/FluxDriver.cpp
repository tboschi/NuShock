#include <FluxDriver.h>

FluxDriver::FluxDriver(std::string FluxConfig) : 
	M_Electron(Const::fMElectron),
	M_Muon(Const::fMMuon),
	M_Pion(Const::fMPion),
	M_Pion0(Const::fMPion0),
	M_Kaon(Const::fMKaon),
	M_Kaon0(Const::fMKaon0)
{
	hTotalSterile = new TH1D("htotal", "Total flux", 100,0,20);
        hPionSterile = new TH1D("hpion", "Pion flux", 100,0,20);
        hKaonSterile = new TH1D("hkaon", "Kaon flux", 100,0,20);
        hKaon0Sterile = new TH1D("hkaon0", "Kaon0 flux", 100,0,20);
        hMuonSterile = new TH1D("hmuon", "Muon flux", 100,0,20);

	hTotalStandard = new TH1D("itotal", "Total flux", 100,0,20);
        hPionStandard = new TH1D("ipion", "Pion flux", 100,0,20);
        hKaonStandard = new TH1D("ikaon", "Kaon flux", 100,0,20);
        hKaon0Standard = new TH1D("ikaon0", "Kaon0 flux", 100,0,20);
        hMuonStandard = new TH1D("imuon", "Muon flux", 100,0,20);

	fxNuMuon = 0;
	fxNuMuonBar = 0;
	fxNuElectron = 0;
	fxNuElectronBar = 0;

	std::string Line, Key, Name;
	std::stringstream ssL;
	Kine = false;
	
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

			if (Key.find("Kinematics") != std::string::npos)
			{
				if (Name.find("on") != std::string::npos)
					Kine = true;
				else Kine = false;
			}
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

	delete hTotalSterile;
	delete hPionSterile;
	delete hKaonSterile;
	delete hKaon0Sterile;
	delete hMuonSterile;
}

void FluxDriver::MakeSterileFlux(double M_Sterile, double U_e, double U_m, double U_t)
{
	hTotalSterile->Reset("ICES");	//Reset before generating
	hPionSterile->Reset("ICES");
	hKaonSterile->Reset("ICES"); 
	hKaon0Sterile->Reset("ICES");
	hMuonSterile->Reset("ICES"); 

	if (fxNuMuon)		//Model on nu mu flux
	{
		Flux sxFlux(*fxNuMuon);

		//pi+ -> mu+ nu_mu
		sxFlux.GetPion()->Scale(U_m*U_m * Kine::ShrockFactor(M_Pion, M_Muon, M_Sterile));
		hPionSterile->Add(sxFlux.GetPion());

		//K+ -> mu+ nu_mu	&&	pi0 mu+ nu_mu
		double KaonFactor = 63.44/(63.44+3.32) * U_m*U_m * Kine::ShrockFactor(M_Kaon, M_Muon, M_Sterile);	//Two body
		if (Kine)
			KaonFactor += 3.32/(63.44+3.32) * U_m*U_m * Kine::ShrockFactor(M_Kaon - M_Pion0, M_Muon, M_Sterile);	//Three body
		else 	KaonFactor += 3.32/(63.44+3.32) * U_m*U_m;	//Simple scaling
		sxFlux.GetKaon()->Scale(KaonFactor);
		hKaonSterile->Add(sxFlux.GetKaon());

		//K0 -> pi- mu+ nu_mu
		if (Kine)
			sxFlux.GetKaon0()->Scale(U_m*U_m * Kine::ShrockFactor(M_Kaon0 - M_Pion, M_Muon, M_Sterile));	//Three body
		else 	sxFlux.GetKaon0()->Scale(U_m*U_m);	//simple scaling
		hKaon0Sterile->Add(sxFlux.GetKaon0());

		//mu -> nu_mu e nu_e
		sxFlux.GetMuon()->Scale(U_m*U_m);
		hMuonSterile->Add(sxFlux.GetMuon());
	}

	if (fxNuMuonBar)	//Model on nu mu bar flux
	{
		Flux sxFlux(*fxNuMuonBar);

		//pi- -> mu- nu_mu_bar
		sxFlux.GetPion()->Scale(U_m*U_m * Kine::ShrockFactor(M_Pion, M_Muon, M_Sterile));
		hPionSterile->Add(sxFlux.GetPion());

		//K- -> mu- nu_mu_bar	&&	pi0 mu- nu_mu_bar
		double KaonFactor = 63.44/(63.44+3.32) * U_m*U_m * Kine::ShrockFactor(M_Kaon, M_Muon, M_Sterile);	//Two body
		if (Kine)
			KaonFactor += 3.32/(63.44+3.32) * U_m*U_m * Kine::ShrockFactor(M_Kaon - M_Pion0, M_Muon, M_Sterile);	//Three body
		else 	KaonFactor += 3.32/(63.44+3.32) * U_m*U_m;	//simple scaling
		sxFlux.GetKaon()->Scale(KaonFactor);
		hKaonSterile->Add(sxFlux.GetKaon());

		//K0 -> pi- mu+ nu_mu
		if (Kine)
			sxFlux.GetKaon0()->Scale(U_m*U_m * Kine::ShrockFactor(M_Kaon0 - M_Pion, M_Muon, M_Sterile));	//Three body
		else 	sxFlux.GetKaon0()->Scale(U_m*U_m);	//simple scaling
		hKaon0Sterile->Add(sxFlux.GetKaon0());

		//mu -> nu_mu e nu_e
		sxFlux.GetMuon()->Scale(U_m*U_m);
		hMuonSterile->Add(sxFlux.GetMuon());
	}
	
	if (fxNuElectron)	//Model on nu electron flux
	{
		Flux sxFlux(*fxNuElectron);

		//pi+ -> e+ nu_e
		sxFlux.GetPion()->Scale(U_e*U_e * Kine::ShrockFactor(M_Pion, M_Electron, M_Sterile));
		hPionSterile->Add(sxFlux.GetPion());

		//K+ -> pi0 e+ nu_e
		if (Kine)
			sxFlux.GetKaon()->Scale(U_e*U_e * Kine::ShrockFactor(M_Kaon - M_Pion0, M_Electron, M_Sterile));	//Three body
		else 	sxFlux.GetKaon()->Scale(U_e*U_e);	//simple scaling
		hKaonSterile->Add(sxFlux.GetKaon());

		//K0 -> pi- e+ nu_e
		if (Kine)
			sxFlux.GetKaon0()->Scale(U_e*U_e * Kine::ShrockFactor(M_Kaon0 - M_Pion, M_Electron, M_Sterile));	//Three body
		else 	sxFlux.GetKaon0()->Scale(U_e*U_e);	//simple scaling
		hKaon0Sterile->Add(sxFlux.GetKaon0());

		//mu -> nu_mu e nu_e
		sxFlux.GetMuon()->Scale(U_e*U_e);
		hMuonSterile->Add(sxFlux.GetMuon());
	}

	if (fxNuElectronBar)	//Model on nu electron bar flux
	{
		Flux sxFlux(*fxNuElectronBar);

		//pi- -> e- nu_e_bar
		sxFlux.GetPion()->Scale(U_e*U_e * Kine::ShrockFactor(M_Pion, M_Electron, M_Sterile));
		hPionSterile->Add(sxFlux.GetPion());

		//K- -> pi0 e- nu_e_bar
		if (Kine)
			sxFlux.GetKaon()->Scale(U_e*U_e * Kine::ShrockFactor(M_Kaon - M_Pion0, M_Electron, M_Sterile));	//Three body
		else 	sxFlux.GetKaon()->Scale(U_e*U_e);	//simple scaling
		hKaonSterile->Add(sxFlux.GetKaon());

		//K0 -> pi- e- nu_e_bar
		if (Kine)
			sxFlux.GetKaon0()->Scale(U_e*U_e * Kine::ShrockFactor(M_Kaon0 - M_Pion, M_Electron, M_Sterile));	//Three body
		else 	sxFlux.GetKaon0()->Scale(U_e*U_e);	//simple scaling
		hKaon0Sterile->Add(sxFlux.GetKaon0());

		//mu -> nu_mu e nu_e
		sxFlux.GetMuon()->Scale(U_e*U_e);
		hMuonSterile->Add(sxFlux.GetMuon());
	}

	hTotalSterile->Add(hPionSterile);
	hTotalSterile->Add(hKaonSterile);
	hTotalSterile->Add(hKaon0Sterile);
	hTotalSterile->Add(hMuonSterile);
}

/*
void FluxDriver::MakeSterilePDF(double Energy, double Prob)
{
	int Bin = hTotalSterile->FindBin(Energy);
	double Entry = hTotalSterile->GetBinContent(Bin) * Prob;
	hTotalSterile->SetBinContent(Bin, Entry);
}
*/

void FluxDriver::MakeStandardFlux()
{
	hTotalStandard->Reset("ICES");	//Reset before generating
	hPionStandard->Reset("ICES");
	hKaonStandard->Reset("ICES"); 
	hKaon0Standard->Reset("ICES");
	hMuonStandard->Reset("ICES"); 

	if (fxNuMuon)		//Model on nu mu flux
	{
		Flux sxFlux = *fxNuMuon;

		hPionStandard->Add(sxFlux.GetPion());
		hKaonStandard->Add(sxFlux.GetKaon());
		hKaon0Standard->Add(sxFlux.GetKaon0());
		hMuonStandard->Add(sxFlux.GetMuon());
	}

	if (fxNuMuonBar)	//Model on nu mu bar flux
	{
		Flux sxFlux = *fxNuMuonBar;

		hPionStandard->Add(sxFlux.GetPion());
		hKaonStandard->Add(sxFlux.GetKaon());
		hKaon0Standard->Add(sxFlux.GetKaon0());
		hMuonStandard->Add(sxFlux.GetMuon());
	}
	
	if (fxNuElectron)	//Model on nu electron flux
	{
		Flux sxFlux = *fxNuElectron;

		hPionStandard->Add(sxFlux.GetPion());
		hKaonStandard->Add(sxFlux.GetKaon());
		hKaon0Standard->Add(sxFlux.GetKaon0());
		hMuonStandard->Add(sxFlux.GetMuon());
	}

	if (fxNuElectronBar)	//Model on nu electron bar flux
	{
		Flux sxFlux = *fxNuElectronBar;

		hPionStandard->Add(sxFlux.GetPion());
		hKaonStandard->Add(sxFlux.GetKaon());
		hKaon0Standard->Add(sxFlux.GetKaon0());
		hMuonStandard->Add(sxFlux.GetMuon());
	}

	hTotalStandard->Add(hPionStandard);
	hTotalStandard->Add(hKaonStandard);
	hTotalStandard->Add(hKaon0Standard);
	hTotalStandard->Add(hMuonStandard);
}

double FluxDriver::SampleEnergy()		//Sample PDF
{
	return hTotalSterile->GetRandom();
} 

double FluxDriver::GetIntensity(double Energy)	//Return flux intensity, given energy
{
	return hTotalSterile->GetBinContent(hTotalSterile->FindBin(Energy));
} 

void FluxDriver::SetBaseline(double Baseline)
{
	hTotalSterile->Scale(1e4/(Baseline*Baseline));
	hPionSterile->Scale(1e4/(Baseline*Baseline));
	hKaonSterile->Scale(1e4/(Baseline*Baseline));
	hKaon0Sterile->Scale(1e4/(Baseline*Baseline));
	hPionSterile->Scale(1e4/(Baseline*Baseline));

	hTotalStandard->Scale(1e4/(Baseline*Baseline));
	hPionStandard->Scale(1e4/(Baseline*Baseline));
	hKaonStandard->Scale(1e4/(Baseline*Baseline));
	hKaon0Standard->Scale(1e4/(Baseline*Baseline));
	hPionStandard->Scale(1e4/(Baseline*Baseline));
}

void FluxDriver::SetPOT(double POT)
{
	hTotalSterile->Scale(POT);
	hPionSterile->Scale(POT);
	hKaonSterile->Scale(POT);
	hKaon0Sterile->Scale(POT);
	hPionSterile->Scale(POT);

	hTotalStandard->Scale(POT);
	hPionStandard->Scale(POT);
	hKaonStandard->Scale(POT);
	hKaon0Standard->Scale(POT);
	hPionStandard->Scale(POT);
}

//Get individual histograms

TH1D* FluxDriver::GetTotal()
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

TH1D* FluxDriver::GetTotalOriginal()
{
	return hTotalStandard;
} 

TH1D* FluxDriver::GetPionOriginal()
{
	return hPionStandard;
} 

TH1D* FluxDriver::GetKaonOriginal()
{
	return hKaonStandard;
} 

TH1D* FluxDriver::GetKaon0Original()
{
	return hKaon0Standard;
} 

TH1D* FluxDriver::GetMuonOriginal()
{
	return hMuonStandard;
} 
