#include <FluxDriver.h>

FluxDriver::FluxDriver(std::string FluxConfig) : 
	M_Electron(Const::fMElectron),
	M_Muon(Const::fMMuon),
	M_Pion(Const::fMPion),
	M_Pion0(Const::fMPion0),
	M_Kaon(Const::fMKaon),
	M_Kaon0(Const::fMKaon0)
{
	BinNumber = 100;
	RangeStart = 0;
	RangeEnd = 20;
	M_Sterile_prev = -1.0;

	//muonic components
	hTotalMn = new TH1D("htotalmn", "Total flux", 100,0,20);
        hPionMn  = new TH1D("hpionmn",  "Pion flux",  100,0,20);
        hKaonMn  = new TH1D("hkaonmn",  "Kaon flux",  100,0,20);
        hKaon0Mn = new TH1D("hkaon0mn", "Kaon0 flux", 100,0,20);
        hMuonMn  = new TH1D("hmuonmn",  "Muon flux",  100,0,20);

	hTotalMa = new TH1D("htotalma", "Total flux", 100,0,20);
        hPionMa  = new TH1D("hpionma",  "Pion flux",  100,0,20);
        hKaonMa  = new TH1D("hkaonma",  "Kaon flux",  100,0,20);
        hKaon0Ma = new TH1D("hkaon0ma", "Kaon0 flux", 100,0,20);
        hMuonMa  = new TH1D("hmuonma",  "Muon flux",  100,0,20);

	//electronic components
	hTotalEn = new TH1D("htotalen", "Total flux", 100,0,20);
        hPionEn  = new TH1D("hpionen",  "Pion flux",  100,0,20);
        hKaonEn  = new TH1D("hkaonen",  "Kaon flux",  100,0,20);
        hKaon0En = new TH1D("hkaon0en", "Kaon0 flux", 100,0,20);
        hMuonEn  = new TH1D("hmuonen",  "Muon flux",  100,0,20);
	
	hTotalEa = new TH1D("htotalea", "Total flux", 100,0,20);
        hPionEa  = new TH1D("hpionea",  "Pion flux",  100,0,20);
        hKaonEa  = new TH1D("hkaonea",  "Kaon flux",  100,0,20);
        hKaon0Ea = new TH1D("hkaon0ea", "Kaon0 flux", 100,0,20);
        hMuonEa  = new TH1D("hmuonea",  "Muon flux",  100,0,20);
	
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

	delete hTotalMn;
	delete hTotalMa;
	delete hPionMn;
	delete hPionMa;
	delete hKaonMn;
	delete hKaonMa;
	delete hKaon0Mn;
	delete hKaon0Ma;
	delete hMuonMn;
	delete hMuonMa;

	delete hTotalEn;
	delete hTotalEa;
	delete hPionEn;
	delete hPionEa;
	delete hKaonEn;
	delete hKaonEa;
	delete hKaon0En;
	delete hKaon0Ea;
	delete hMuonEn;
	delete hMuonEa;
}

bool FluxDriver::MakeFlux(double M_Sterile)
{
	if (IsChanged(M_Sterile))
	{
		hTotalMn->Reset("ICES");	//Reset before generating
		hPionMn->Reset("ICES");
		hKaonMn->Reset("ICES"); 
		hKaon0Mn->Reset("ICES");
		hMuonMn->Reset("ICES"); 

		hTotalEn->Reset("ICES");	//Reset before generating
		hPionEn->Reset("ICES");
		hKaonEn->Reset("ICES"); 
		hKaon0En->Reset("ICES");
		hMuonEn->Reset("ICES"); 

		hTotalMa->Reset("ICES");	//Reset before generating
		hPionMa->Reset("ICES");
		hKaonMa->Reset("ICES"); 
		hKaon0Ma->Reset("ICES");
		hMuonMa->Reset("ICES"); 

		hTotalEa->Reset("ICES");	//Reset before generating
		hPionEa->Reset("ICES");
		hKaonEa->Reset("ICES"); 
		hKaon0Ea->Reset("ICES");
		hMuonEa->Reset("ICES"); 

		if (fxNuMuon)		//Model on nu mu flux
		{
			Flux sxFlux(*fxNuMuon);
			MakeMuonComponent(1, sxFlux, M_Sterile);

			hTotalMn->Add(hPionMn);
			hTotalMn->Add(hKaonMn);
			hTotalMn->Add(hKaon0Mn);
			hTotalMn->Add(hMuonMn);
		}
	
		if (fxNuMuonBar)	//Model on nu mu bar flux
		{
			Flux sxFlux(*fxNuMuonBar);
			MakeMuonComponent(0, sxFlux, M_Sterile);

			hTotalMa->Add(hPionMa);
			hTotalMa->Add(hKaonMa);
			hTotalMa->Add(hKaon0Ma);
			hTotalMa->Add(hMuonMa);
		}
	
		if (fxNuElectron)	//Model on nu electron flux
		{
			Flux sxFlux(*fxNuElectron);
			MakeElecComponent(1, sxFlux, M_Sterile);

			hTotalEn->Add(hPionEn);
			hTotalEn->Add(hKaonEn);
		 	hTotalEn->Add(hKaon0En);
			hTotalEn->Add(hMuonEn);
		}
	
		if (fxNuElectronBar)	//Model on nu electron bar flux
		{
			Flux sxFlux(*fxNuElectronBar);
			MakeElecComponent(0, sxFlux, M_Sterile);

			hTotalEa->Add(hPionEa);
			hTotalEa->Add(hKaonEa);
		 	hTotalEa->Add(hKaon0Ea);
			hTotalEa->Add(hMuonEa);
		}

		return true;
	}
	else 
		return false;
}

void FluxDriver::MakeMuonComponent(bool Neutrino, Flux &sxFlux, double M_Sterile)
{
	//pi+ -> mu+ nu_mu
	sxFlux.GetPion()->Scale(Kine::ShrockFactor(M_Pion, M_Muon, M_Sterile));
	if (Neutrino)
		hPionMn->Add(sxFlux.GetPion());
	else 
		hPionMa->Add(sxFlux.GetPion());

	//K+ -> mu+ nu_mu	&&	pi0 mu+ nu_mu
	double KaonFactor = 63.56/(63.56+3.35) * Kine::ShrockFactor(M_Kaon, M_Muon, M_Sterile);	//Two body
	if (Kine)
	{
		double KM = hKaonMuon->GetBinContent(hKaonMuon->FindBin(M_Sterile+1e-9)+1);	//1e-9 to prevent bin error
		KaonFactor += 3.35/(63.56+3.35) * KM;	//Three body
	}
	else if (M_Sterile < M_Kaon - M_Pion0 - M_Muon)
		KaonFactor += 3.35/(63.56+3.35);	//simple scaling
	else KaonFactor += 0;	//simple scaling
	sxFlux.GetKaon()->Scale(KaonFactor);
	if (Neutrino)
		hKaonMn->Add(sxFlux.GetKaon());
	else
		hKaonMa->Add(sxFlux.GetKaon());

	//K0 -> pi- mu+ nu_mu
	if (Kine)
	{
		double K0M = hKaon0Muon->GetBinContent(hKaon0Muon->FindBin(M_Sterile+1e-9));	//1e-9 to prevent bin error
		sxFlux.GetKaon0()->Scale(K0M);	//Three body
	}
	else if (M_Sterile > M_Kaon0 - M_Pion - M_Muon)
		sxFlux.GetKaon0()->Scale(0);	//simple scaling
	if (Neutrino)
		hKaon0Mn->Add(sxFlux.GetKaon0());
	else
		hKaon0Ma->Add(sxFlux.GetKaon0());

	//mu -> nu_mu e nu_e
	if (Kine)
	{
		double MM = hMuonMuon->GetBinContent(hMuonMuon->FindBin(M_Sterile+1e-9));	//1e-9 to prevent bin error
		sxFlux.GetMuon()->Scale(MM);	//Three body
	}
	else if (M_Sterile > M_Muon - M_Electron)
		sxFlux.GetMuon()->Scale(0);
	if (Neutrino)
		hMuonMn->Add(sxFlux.GetMuon());
	else
		hMuonMa->Add(sxFlux.GetMuon());
}

void FluxDriver::MakeElecComponent(bool Neutrino, Flux &sxFlux, double M_Sterile)
{
	//pi- -> e- nu_e_bar
	sxFlux.GetPion()->Scale(Kine::ShrockFactor(M_Pion, M_Electron, M_Sterile));
	if (Neutrino)
		hPionEn->Add(sxFlux.GetPion());
	else
		hPionEa->Add(sxFlux.GetPion());

	//K- -> pi0 e- nu_e_bar
	double KaonFactor = 1.582e-3/(1.582e-3+5.07) * Kine::ShrockFactor(M_Kaon, M_Electron, M_Sterile);	//Two body
	if (Kine)
	{
		double KE = hKaonElec->GetBinContent(hKaonElec->FindBin(M_Sterile+1e-9)+1);	//1e-9 to prevent bin error
		KaonFactor += 5.07/(1.582e-3+5.07) * KE;	//Three body
	}
	else if (M_Sterile < M_Kaon - M_Pion0 - M_Electron)
		KaonFactor += 5.07/(1.582e-3+5.07);	//simple scaling
	else KaonFactor += 0;	//simple scaling
	sxFlux.GetKaon()->Scale(KaonFactor);
	if (Neutrino)
		hKaonEn->Add(sxFlux.GetKaon());
	else
		hKaonEa->Add(sxFlux.GetKaon());

	//K0 -> pi- e- nu_e_bar
	if (Kine)
	{
		double K0E = hKaon0Elec->GetBinContent(hKaon0Elec->FindBin(M_Sterile+1e-9));	//1e-9 to prevent bin error
		sxFlux.GetKaon0()->Scale(K0E);
	}
	else if (M_Sterile > M_Kaon0 - M_Pion - M_Electron)
		sxFlux.GetKaon0()->Scale(0);	//simple scaling
	if (Neutrino)
		hKaon0En->Add(sxFlux.GetKaon0());
	else
		hKaon0Ea->Add(sxFlux.GetKaon0());

	//mu -> nu_mu e nu_e
	if (Kine)
	{
		double ME = hMuonElec->GetBinContent(hMuonElec->FindBin(M_Sterile+1e-9));	//1e-9 to prevent bin error
		sxFlux.GetMuon()->Scale(ME);
	}
	else if (M_Sterile > M_Muon - M_Electron)
		sxFlux.GetMuon()->Scale(0);
	if (Neutrino)
		hMuonEn->Add(sxFlux.GetMuon());
	else
		hMuonEa->Add(sxFlux.GetMuon());
}

/*
double FluxDriver::SampleEnergy()		//Sample PDF
{
	return hTotalM->GetRandom();
}
*/

double FluxDriver::GetRangeStart()
{
	return RangeStart;
}

double FluxDriver::GetRangeEnd()
{
	return RangeEnd;
}

int FluxDriver::GetBinNumber()
{
	return BinNumber;
}

double FluxDriver::GetIntensityNeut(double Energy, double Ue, double Um, double Ut)	//Return flux intensity, given energy, simple linear interpolation
{
	int Bin = hTotalMn->FindBin(Energy);
	double I1 = hTotalMn->GetBinContent(Bin);
	double E1 = hTotalMn->GetBinCenter(Bin);
	double I2, E2;
	if (Energy < hTotalMn->GetBinCenter(Bin))
	{
		I2 = hTotalMn->GetBinContent(Bin-1);
		E2 = hTotalMn->GetBinCenter(Bin-1);
	}
	else
	{
		I2 = hTotalMn->GetBinContent(Bin+1);
		E2 = hTotalMn->GetBinCenter(Bin+1);
	}

	double IntenseM = (Energy-E1)*(I2-I1)/(E2-E1) + I1;

	Bin = hTotalEn->FindBin(Energy);
	I1 = hTotalEn->GetBinContent(Bin);
	E1 = hTotalEn->GetBinCenter(Bin);
	I2, E2;
	if (Energy < hTotalEn->GetBinCenter(Bin))
	{
		I2 = hTotalEn->GetBinContent(Bin-1);
		E2 = hTotalEn->GetBinCenter(Bin-1);
	}
	else
	{
		I2 = hTotalEn->GetBinContent(Bin+1);
		E2 = hTotalEn->GetBinCenter(Bin+1);
	}

	double IntenseE = (Energy-E1)*(I2-I1)/(E2-E1) + I1;

	return Um*Um*IntenseM + Ue*Ue*IntenseE;
} 

double FluxDriver::GetIntensityAnti(double Energy, double Ue, double Um, double Ut)	//Return flux intensity, given energy, simple linear interpolation
{
	int Bin = hTotalMa->FindBin(Energy);
	double I1 = hTotalMa->GetBinContent(Bin);
	double E1 = hTotalMa->GetBinCenter(Bin);
	double I2, E2;
	if (Energy < hTotalMa->GetBinCenter(Bin))
	{
		I2 = hTotalMa->GetBinContent(Bin-1);
		E2 = hTotalMa->GetBinCenter(Bin-1);
	}
	else
	{
		I2 = hTotalMa->GetBinContent(Bin+1);
		E2 = hTotalMa->GetBinCenter(Bin+1);
	}

	double IntenseM = (Energy-E1)*(I2-I1)/(E2-E1) + I1;

	Bin = hTotalEa->FindBin(Energy);
	I1 = hTotalEa->GetBinContent(Bin);
	E1 = hTotalEa->GetBinCenter(Bin);
	I2, E2;
	if (Energy < hTotalEa->GetBinCenter(Bin))
	{
		I2 = hTotalEa->GetBinContent(Bin-1);
		E2 = hTotalEa->GetBinCenter(Bin-1);
	}
	else
	{
		I2 = hTotalEa->GetBinContent(Bin+1);
		E2 = hTotalEa->GetBinCenter(Bin+1);
	}

	double IntenseE = (Energy-E1)*(I2-I1)/(E2-E1) + I1;

	return Um*Um*IntenseM + Ue*Ue*IntenseE;
} 

void FluxDriver::SetBaseline(double Baseline)
{
	hTotalMn->Scale(1.0/(Baseline*Baseline));
	hPionMn->Scale(1.0/(Baseline*Baseline));
	hKaonMn->Scale(1.0/(Baseline*Baseline));
	hKaon0Mn->Scale(1.0/(Baseline*Baseline));
	hMuonMn->Scale(1.0/(Baseline*Baseline));

	hTotalMa->Scale(1.0/(Baseline*Baseline));
	hPionMa->Scale(1.0/(Baseline*Baseline));
	hKaonMa->Scale(1.0/(Baseline*Baseline));
	hKaon0Ma->Scale(1.0/(Baseline*Baseline));
	hMuonMa->Scale(1.0/(Baseline*Baseline));

	hTotalEn->Scale(1.0/(Baseline*Baseline));
	hPionEn->Scale(1.0/(Baseline*Baseline));
	hKaonEn->Scale(1.0/(Baseline*Baseline));
	hKaon0En->Scale(1.0/(Baseline*Baseline));
	hMuonEn->Scale(1.0/(Baseline*Baseline));

	hTotalEa->Scale(1.0/(Baseline*Baseline));
	hPionEa->Scale(1.0/(Baseline*Baseline));
	hKaonEa->Scale(1.0/(Baseline*Baseline));
	hKaon0Ea->Scale(1.0/(Baseline*Baseline));
	hMuonEa->Scale(1.0/(Baseline*Baseline));
}

void FluxDriver::SetPOT(double POT)
{
	hTotalMn->Scale(POT);
	hPionMn->Scale(POT);
	hKaonMn->Scale(POT);
	hKaon0Mn->Scale(POT);
	hMuonMn->Scale(POT);

	hTotalMa->Scale(POT);
	hPionMa->Scale(POT);
	hKaonMa->Scale(POT);
	hKaon0Ma->Scale(POT);
	hMuonMa->Scale(POT);

	hTotalEn->Scale(POT);
	hPionEn->Scale(POT);
	hKaonEn->Scale(POT);
	hKaon0En->Scale(POT);
	hMuonEn->Scale(POT);

	hTotalEa->Scale(POT);
	hPionEa->Scale(POT);
	hKaonEa->Scale(POT);
	hKaon0Ea->Scale(POT);
	hMuonEa->Scale(POT);
}

void FluxDriver::SetArea(double Area)
{
	hTotalMn->Scale(Area);
	hPionMn->Scale(Area);
	hKaonMn->Scale(Area);
	hKaon0Mn->Scale(Area);
	hMuonMn->Scale(Area);

	hTotalMa->Scale(Area);
	hPionMa->Scale(Area);
	hKaonMa->Scale(Area);
	hKaon0Ma->Scale(Area);
	hMuonMa->Scale(Area);

	hTotalEn->Scale(Area);
	hPionEn->Scale(Area);
	hKaonEn->Scale(Area);
	hKaon0En->Scale(Area);
	hMuonEn->Scale(Area);

	hTotalEa->Scale(Area);
	hPionEa->Scale(Area);
	hKaonEa->Scale(Area);
	hKaon0Ea->Scale(Area);
	hMuonEa->Scale(Area);
}


bool FluxDriver::IsChanged(double M_Sterile)
{
	bool Ret = (fabs(M_Sterile - M_Sterile_prev) > 1e-9);
	M_Sterile_prev = M_Sterile;

	return Ret;
}

//Get individual histograms
/*
TH1D* FluxDriver::GetTotalMn()
{
	return hTotalMn;
} 

TH1D* FluxDriver::GetPionM()
{
	return hPionM;
} 

TH1D* FluxDriver::GetKaonM()
{
	return hKaonM;
} 

TH1D* FluxDriver::GetKaon0M()
{
	return hKaon0M;
} 

TH1D* FluxDriver::GetMuonM()
{
	return hMuonM;
} 

TH1D* FluxDriver::GetTotalE()
{
	return hTotalE;
} 

TH1D* FluxDriver::GetPionE()
{
	return hPionE;
} 

TH1D* FluxDriver::GetKaonE()
{
	return hKaonE;
} 

TH1D* FluxDriver::GetKaon0E()
{
	return hKaon0E;
} 

TH1D* FluxDriver::GetMuonE()
{
	return hMuonE;
} 
*/
