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
	
	//hTotalSterile = new TH1D("htotal", "Total flux", 100,0,5);	//BNB
        //hPionSterile = new TH1D("hpion", "Pion flux", 100,0,5);	//BNB
        //hKaonSterile = new TH1D("hkaon", "Kaon flux", 100,0,5);	//BNB 
        //hKaon0Sterile = new TH1D("hkaon0", "Kaon0 flux", 100,0,5);	//BNB
        //hMuonSterile = new TH1D("hmuon", "Muon flux", 100,0,5);	//BNB

	hTotalStandard = new TH1D("itotal", "Total flux", 100,0,20);
        hPionStandard = new TH1D("ipion", "Pion flux", 100,0,20);
        hKaonStandard = new TH1D("ikaon", "Kaon flux", 100,0,20);
        hKaon0Standard = new TH1D("ikaon0", "Kaon0 flux", 100,0,20);
        hMuonStandard = new TH1D("imuon", "Muon flux", 100,0,20);

	fxNuMuon = 0;
	fxNuMuonBar = 0;
	fxNuElectron = 0;
	fxNuElectronBar = 0;

	std::stringstream ssL;
	std::string Line, Key, Name;
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
				KineFile = new TFile(Name.c_str(), "OPEN");
				Kine = true;

				hTemp = (TH1D*) KineFile->Get("muonmuon");
				hMuonMuon  = (TH1D*) hTemp->Clone();
				hMuonMuon->SetDirectory(0);

				hTemp = (TH1D*) KineFile->Get("muonelec");
				hMuonElec  = (TH1D*) hTemp->Clone();
				hMuonElec->SetDirectory(0);

				hTemp = (TH1D*) KineFile->Get("kaonmuon");
				hKaonMuon  = (TH1D*) hTemp->Clone();
				hKaonMuon->SetDirectory(0);

				hTemp = (TH1D*) KineFile->Get("kaonelec");
				hKaonElec  = (TH1D*) hTemp->Clone();
				hKaonElec->SetDirectory(0);
				
				hTemp = (TH1D*) KineFile->Get("kaon0muon");
				hKaon0Muon = (TH1D*) hTemp->Clone();
				hKaon0Muon->SetDirectory(0);
				
				hTemp = (TH1D*) KineFile->Get("kaon0elec");
				hKaon0Elec = (TH1D*) hTemp->Clone();
				hKaon0Elec->SetDirectory(0);

				KineFile->Close();
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

bool FluxDriver::MakeSterileFlux(double M_Sterile, double U_e, double U_m, double U_t)
{
	bool Ret = false;
	if (IsChanged(M_Sterile, U_e, U_m, U_t))
	{
		hTotalSterile->Reset("ICES");	//Reset before generating
		hPionSterile->Reset("ICES");
		hKaonSterile->Reset("ICES"); 
		hKaon0Sterile->Reset("ICES");
		hMuonSterile->Reset("ICES"); 

		if (fxNuMuon)		//Model on nu mu flux
		{
			Flux sxFlux(*fxNuMuon);
			MakeMuonComponent(sxFlux, M_Sterile, U_e, U_m, U_t);
		}
	
		if (fxNuMuonBar)	//Model on nu mu bar flux
		{
			Flux sxFlux(*fxNuMuonBar);
			MakeMuonComponent(sxFlux, M_Sterile, U_e, U_m, U_t);
		}
	
		if (fxNuElectron)	//Model on nu electron flux
		{
			Flux sxFlux(*fxNuElectron);
			MakeElecComponent(sxFlux, M_Sterile, U_e, U_m, U_t);
		}
	
		if (fxNuElectronBar)	//Model on nu electron bar flux
		{
			Flux sxFlux(*fxNuElectronBar);
			MakeElecComponent(sxFlux, M_Sterile, U_e, U_m, U_t);
		}
	
		Ret = hTotalSterile->Add(hPionSterile);
		Ret *= hTotalSterile->Add(hKaonSterile);
		Ret *= hTotalSterile->Add(hKaon0Sterile);
		Ret *= hTotalSterile->Add(hMuonSterile);
	}

	return Ret;
}

void FluxDriver::MakeMuonComponent(Flux &sxFlux, double M_Sterile, double U_e, double U_m, double U_t)
{
	//pi+ -> mu+ nu_mu
	sxFlux.GetPion()->Scale(U_m*U_m * Kine::ShrockFactor(M_Pion, M_Muon, M_Sterile));
	hPionSterile->Add(sxFlux.GetPion());

	//K+ -> mu+ nu_mu	&&	pi0 mu+ nu_mu
	double KaonFactor = 63.56/(63.56+3.35) * U_m*U_m * Kine::ShrockFactor(M_Kaon, M_Muon, M_Sterile);	//Two body
	if (Kine)
	{
		double KM = hKaonMuon->GetBinContent(hKaonMuon->FindBin(M_Sterile+1e-9)+1);	//1e-9 to prevent bin error
		KaonFactor += 3.35/(63.56+3.35) * U_m*U_m * KM;	//Three body
	}
	else if (M_Sterile < M_Kaon - M_Pion0 - M_Muon)
		KaonFactor += 3.35/(63.56+3.35) * U_m*U_m;	//simple scaling
	else KaonFactor += 0;	//simple scaling
	sxFlux.GetKaon()->Scale(KaonFactor);
	hKaonSterile->Add(sxFlux.GetKaon());

	//K0 -> pi- mu+ nu_mu
	if (Kine)
	{
		double K0M = hKaon0Muon->GetBinContent(hKaon0Muon->FindBin(M_Sterile+1e-9));	//1e-9 to prevent bin error
		sxFlux.GetKaon0()->Scale(U_m*U_m * K0M);	//Three body
	}
	else if (M_Sterile < M_Kaon0 - M_Pion - M_Muon)
		sxFlux.GetKaon0()->Scale(U_m*U_m);	//simple scaling
	else 
		sxFlux.GetKaon0()->Scale(0);	//simple scaling
	hKaon0Sterile->Add(sxFlux.GetKaon0());

	//mu -> nu_mu e nu_e
	if (Kine)
	{
		double MM = hMuonMuon->GetBinContent(hMuonMuon->FindBin(M_Sterile+1e-9));	//1e-9 to prevent bin error
		sxFlux.GetMuon()->Scale(U_m*U_m * MM);	//Three body
	}
	else if (M_Sterile < M_Muon - M_Electron)
		sxFlux.GetMuon()->Scale(U_m*U_m);
	else sxFlux.GetMuon()->Scale(0);
	hMuonSterile->Add(sxFlux.GetMuon());
}

void FluxDriver::MakeElecComponent(Flux &sxFlux, double M_Sterile, double U_e, double U_m, double U_t)
{
	//pi- -> e- nu_e_bar
	sxFlux.GetPion()->Scale(U_e*U_e * Kine::ShrockFactor(M_Pion, M_Electron, M_Sterile));
	hPionSterile->Add(sxFlux.GetPion());

	//K- -> pi0 e- nu_e_bar
	double KaonFactor = 1.582e-3/(1.582e-3+5.07) * U_e*U_e * Kine::ShrockFactor(M_Kaon, M_Electron, M_Sterile);	//Two body
	if (Kine)
	{
		double KE = hKaonElec->GetBinContent(hKaonElec->FindBin(M_Sterile+1e-9)+1);	//1e-9 to prevent bin error
		KaonFactor += 5.07/(1.582e-3+5.07) * U_e*U_e * KE;	//Three body
	}
	else if (M_Sterile < M_Kaon - M_Pion0 - M_Electron)
		KaonFactor += 5.07/(1.582e-3+5.07) * U_e*U_e;	//simple scaling
	else KaonFactor += 0;	//simple scaling
	sxFlux.GetKaon()->Scale(KaonFactor);
	hKaonSterile->Add(sxFlux.GetKaon());

	//K0 -> pi- e- nu_e_bar
	if (Kine)
	{
		double K0E = hKaon0Elec->GetBinContent(hKaon0Elec->FindBin(M_Sterile+1e-9));	//1e-9 to prevent bin error
		sxFlux.GetKaon0()->Scale(U_e*U_e * K0E);
	}
	else if (M_Sterile < M_Kaon0 - M_Pion - M_Electron)
		sxFlux.GetKaon0()->Scale(U_e*U_e);	//simple scaling
	else sxFlux.GetKaon0()->Scale(0);	//simple scaling
	hKaon0Sterile->Add(sxFlux.GetKaon0());

	//mu -> nu_mu e nu_e
	if (Kine)
	{
		double ME = hMuonElec->GetBinContent(hMuonElec->FindBin(M_Sterile+1e-9));	//1e-9 to prevent bin error
		sxFlux.GetMuon()->Scale(U_e*U_e * ME);
	}
	else if (M_Sterile < M_Muon - M_Electron)
		sxFlux.GetMuon()->Scale(U_e*U_e);
	else sxFlux.GetMuon()->Scale(0);
	hMuonSterile->Add(sxFlux.GetMuon());
}

/*
void FluxDriver::MakeSterilePDF(double Energy, double Prob)
{
	int Bin = hTotalSterile->FindBin(Energy);
	double Entry = hTotalSterile->GetBinContent(Bin) * Prob;
	hTotalSterile->SetBinContent(Bin, Entry);
}
*/

/*
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
*/

double FluxDriver::SampleEnergy()		//Sample PDF
{
	return hTotalSterile->GetRandom();
} 

double FluxDriver::GetRange()
{
	return hTotalSterile->GetBinCenter(GetBinNumber()) - hTotalSterile->GetBinCenter(0);
}

double FluxDriver::GetStartRange()
{
	return hTotalSterile->GetBinCenter(1);
}

double FluxDriver::GetEndRange()
{
	return hTotalSterile->GetBinCenter(hTotalSterile->GetNbinsX());
}

double FluxDriver::GetBinNumber()
{
	return hTotalSterile->GetNbinsX();
}

long double FluxDriver::GetIntensity(double Energy)	//Return flux intensity, given energy, simple linear interpolation
{
	int Bin = hTotalSterile->FindBin(Energy);
	double f1 = hTotalSterile->GetBinContent(Bin);
	double E1 = hTotalSterile->GetBinCenter(Bin);
	double f2;
	double E2;
	if (Energy < hTotalSterile->GetBinCenter(Bin))
	{
		f2 = hTotalSterile->GetBinContent(Bin-1);
		E2 = hTotalSterile->GetBinCenter(Bin-1);
	}
	else
	{
		f2 = hTotalSterile->GetBinContent(Bin+1);
		E2 = hTotalSterile->GetBinCenter(Bin+1);
	}

	return (Energy-E1)*(f2-f1)/(E2-E1) + f1;
	//return hTotalSterile->GetBinContent(hTotalSterile->FindBin(Energy));	//1e-6 to prevent bin error
} 

void FluxDriver::SetBaseline(double Baseline)
{
	hTotalSterile->Scale(1.0/(Baseline*Baseline));
	hPionSterile->Scale(1.0/(Baseline*Baseline));
	hKaonSterile->Scale(1.0/(Baseline*Baseline));
	hKaon0Sterile->Scale(1.0/(Baseline*Baseline));
	hMuonSterile->Scale(1.0/(Baseline*Baseline));
}

void FluxDriver::SetPOT(double POT)
{
	hTotalSterile->Scale(POT);
	hPionSterile->Scale(POT);
	hKaonSterile->Scale(POT);
	hKaon0Sterile->Scale(POT);
	hMuonSterile->Scale(POT);

	//hTotalStandard->Scale(POT);
	//hPionStandard->Scale(POT);
	//hKaonStandard->Scale(POT);
	//hKaon0Standard->Scale(POT);
	//hMuonStandard->Scale(POT);
}

void FluxDriver::SetArea(double Area)
{
	hTotalSterile->Scale(Area);
	hPionSterile->Scale(Area);
	hKaonSterile->Scale(Area);
	hKaon0Sterile->Scale(Area);
	hMuonSterile->Scale(Area);

	//hTotalStandard->Scale(Area);
	//hPionStandard->Scale(Area);
	//hKaonStandard->Scale(Area);
	//hKaon0Standard->Scale(Area);
	//hMuonStandard->Scale(Area);
}


bool FluxDriver::IsChanged(double M_Sterile, double U_e, double U_m, double U_t)
{
	bool Ret = ( M_Sterile != M_Sterile_prev || 
		     U_e != U_e_prev ||
    		     U_m != U_m_prev ||
    		     U_t != U_t_prev );

	M_Sterile_prev = M_Sterile;
	U_e_prev = U_e;
    	U_m_prev = U_m;
    	U_t_prev = U_t;

	return Ret;
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
