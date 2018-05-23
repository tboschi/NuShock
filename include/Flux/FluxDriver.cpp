#include <Flux/FluxDriver.h>

FluxDriver::FluxDriver(std::string FluxConfig) : 
	M_Electron(Const::fMElectron),
	M_Muon(Const::fMMuon),
	M_Tau(Const::fMTau),
	M_Pion(Const::fMPion),
	M_Pion0(Const::fMPion0),
	M_Kaon(Const::fMKaon),
	M_Kaon0(Const::fMKaon0),
	M_Charm(Const::fMDs)
{
	fxNuElectron = 0;
	fxNuElectronBar = 0;
	fxNuMuon = 0;
	fxNuMuonBar = 0;
	fxNuTau = 0;
	fxNuTauBar = 0;

	std::stringstream ssL;
	std::string Line, Key, Name;
	Kine = false;
	Mod = false;
	
	std::ifstream ModFile;
	std::ifstream ConfigFile(FluxConfig.c_str());
	while (std::getline(ConfigFile, Line))
	{
		if (Line[0] == '#') continue;

		ssL.str("");
		ssL.clear();
		ssL << Line;

		//checks if the keyword is in the config file first, then initialiase the correct object
		//
		if (ssL >> Key >> Name)
		{
			if (Key.find("Electron_")    != std::string::npos) fxNuElectron = new Flux(Name);
			if (Key.find("ElectronBar_") != std::string::npos) fxNuElectronBar = new Flux(Name);
			if (Key.find("Muon_")	     != std::string::npos) fxNuMuon = new Flux(Name);
			if (Key.find("MuonBar_")     != std::string::npos) fxNuMuonBar = new Flux(Name);
			if (Key.find("Tau_")	     != std::string::npos) fxNuTau = new Flux(Name);
			if (Key.find("TauBar_")	     != std::string::npos) fxNuTauBar = new Flux(Name);

			//histogram details (number of bins, start of range, end of range)
			//
			if (Key.find("BinNumber")    != std::string::npos) BinNumber = std::strtoul(Name.c_str(), NULL, 10);
			if (Key.find("RangeStart")   != std::string::npos) RangeStart = std::strtod(Name.c_str(), NULL);
			if (Key.find("RangeEnd")     != std::string::npos) RangeEnd = std::strtod(Name.c_str(), NULL);

			//phase space functions for 3body decays are loaded here
			//
			if (Key.find("Kinematics") != std::string::npos) 
			{
				KineFile = new TFile(Name.c_str(), "OPEN");
				Kine = true;

				CloneCopy(hMuonElec,  KineFile->Get("muonelec"));
				CloneCopy(hMuonMuon,  KineFile->Get("muonmuon"));
				CloneCopy(hKaonElec,  KineFile->Get("kaonelec"));
				CloneCopy(hKaonMuon,  KineFile->Get("kaonmuon"));
				CloneCopy(hKaon0Elec, KineFile->Get("kaon0elec"));
				CloneCopy(hKaon0Muon, KineFile->Get("kaon0muon"));
				CloneCopy(hTauEElec,  KineFile->Get("taueelec"));
				CloneCopy(hTauETau,   KineFile->Get("tauetau"));
				CloneCopy(hTauMMuon,  KineFile->Get("taummuon"));
				CloneCopy(hTauMTau,   KineFile->Get("taumtau"));

				KineFile->Close();
			}

			//correction for charm to tau decay
			//
			if (Key.find("Modifier") != std::string::npos)
			{
				Mod = true;
				ModFile.open(Name.c_str());
			}
		}
	}
	ConfigFile.close();

	//initialising histograms
	//electronic components
	hTotalE = new TH1D("htotale", "Total flux", BinNumber, RangeStart, RangeEnd);
        hPionE  = new TH1D("hpione",  "Pion flux",  BinNumber, RangeStart, RangeEnd);
        hKaonE  = new TH1D("hkaone",  "Kaon flux",  BinNumber, RangeStart, RangeEnd);
        hKaon0E = new TH1D("hkaon0e", "Kaon0 flux", BinNumber, RangeStart, RangeEnd);
        hMuonE  = new TH1D("hmuone",  "Muon flux",  BinNumber, RangeStart, RangeEnd);
        hCharmE = new TH1D("hcharme", "Charm flux", BinNumber, RangeStart, RangeEnd);
	
	//muonic components
	hTotalM = new TH1D("htotalm", "Total flux", BinNumber, RangeStart, RangeEnd);
        hPionM  = new TH1D("hpionm",  "Pion flux",  BinNumber, RangeStart, RangeEnd);
        hKaonM  = new TH1D("hkaonm",  "Kaon flux",  BinNumber, RangeStart, RangeEnd);
        hKaon0M = new TH1D("hkaon0m", "Kaon0 flux", BinNumber, RangeStart, RangeEnd);
        hMuonM  = new TH1D("hmuonm",  "Muon flux",  BinNumber, RangeStart, RangeEnd);
        hCharmM = new TH1D("hcharmm", "Charm flux", BinNumber, RangeStart, RangeEnd);

	//tauonic components
	hTotalT = new TH1D("htotalt", "Total flux", BinNumber, RangeStart, RangeEnd);
        hPionT  = new TH1D("hpiont",  "Pion flux",  BinNumber, RangeStart, RangeEnd);
        h2PionT = new TH1D("h2piont", "2Pion flux", BinNumber, RangeStart, RangeEnd);
        hCharmT = new TH1D("hcharmt", "Charm flux", BinNumber, RangeStart, RangeEnd);
        hTauET  = new TH1D("htauet",  "TauE flux",  BinNumber, RangeStart, RangeEnd);
        hTauMT  = new TH1D("htaumt",  "TauM flux",  BinNumber, RangeStart, RangeEnd);

	M_Sterile_prev = -1.0;

	//load charm to tau modificator to vectors
	//
	if (ModFile.is_open())
	{
		double Mass, X, Y;
		while (getline(ModFile, Line))
		{
			ssL.str("");
			ssL.clear();
			ssL << Line;

			ssL >> Mass >> X >> Y;
			vMdir.push_back(Mass);
			vXdir.push_back(X);
			vYdir.push_back(Y);
		}
		ModFile.close();
	}
}

//deconstructor
FluxDriver::~FluxDriver()
{
	delete fxNuElectron;
	delete fxNuElectronBar;
	delete fxNuMuon;
	delete fxNuMuonBar;
	delete fxNuTau;
	delete fxNuTauBar;

	delete hTotalE;
	delete hPionE;
	delete hKaonE;
	delete hKaon0E;
	delete hMuonE;
	delete hCharmE;

	delete hTotalM;
	delete hPionM;
	delete hKaonM;
	delete hKaon0M;
	delete hMuonM;
	delete hCharmM;

	delete hTotalT;
	delete hPionT;
	delete hPionTa;
	delete hCharmT;
	delete hTauET;
	delete hTauET;
}

//clones the histograms in the kinematic files
void FluxDriver::CloneCopy(TH1D*& T, TObject* X)
{
	if (X)
	{
		T = dynamic_cast<TH1D*> (X->Clone());
		T->SetDirectory(0);
	}
	else
		T = NULL;
}

//make flux, only input needed is the mass of the neutrino
//return true if successful
//
bool FluxDriver::MakeFlux(Neutrino * N)
{
	if (!IsChanged(N))	//compute only if particle is changed have changed 
	{
		std::cerr << "WARNING: recopmuting flux with same neutrino. This should not be happening!" << std::endl;
		return false;
	}
	else
	{
		delete sxHeavyElectron;
		delete sxHeavyMuon;
		delete sxHeavyTau;
		
		if (sxHeavyElectron = new Flux(*fxNuElectron))
			MakeElecComponent(sxHeavyElectron);

		if (sxHeavyMuon = new Flux(*fxNuMuon))
			MakeElecComponent(sxHeavyMuon);

		if (sxHeavyTau = new Flux(*fxNuTau))
			MakeElecComponent(sxHeavyTau);

		return true;
	}
}

//Make electronic components, requires Neutrino (T if neutrino, F if antineutrino), the Flux object, the mass
void FluxDriver::MakeElecComponent(Flux *sxFlux)
{
	TH1D *hPoint, *hTotal;

	//pi+ -> e+ nu_e
	if (hPoint = sxFlux->Get(Hist::Pion))
	{
		hPoint->Scale(Kine::Unhelicity(M_Pion, M_Electron, GetMass(), GetHelicity()));

		hTotal->Add(hPoint);
	}

	//K+ -> pi0 e+ nu_e	(5.07 %)		//I just assume that both decays are equally probable 
	//K+ -> e+ nu_e		(1.582e-3)		//for each different energy, so it is a linear combination of the two
	if (hPoint = sxFlux->Get(Hist::Kaon))
	{
				    //branching ratios
		double KaonFactor = 1.582e-3/(1.582e-3+5.07) * Kine::Unhelicity(M_Kaon, M_Electron, GetMass(), GetHelicity());	//Two body
		if (Kine)
		{
			double KE = hKaonElec->GetBinContent(hKaonElec->FindBin(GetMass()+1e-9)+1);	//1e-9 to prevent bin error
			KaonFactor += 5.07/(1.582e-3+5.07) * KE;	//Three body
		}
		else if (GetMass() < M_Kaon - M_Pion0 - M_Electron)
			KaonFactor += 5.07/(1.582e-3+5.07);	//simple scaling
		else KaonFactor += 0;	//simple scaling
		hPoint->Scale(KaonFactor);

		hTotal->Add(hPoint);
	}

	//K0 -> pi+ e+ nu_e
	if (hPoint = sxFlux->Get(Hist::Kaon0))
	{
		if (Kine)
		{
			double K0E = hKaon0Elec->GetBinContent(hKaon0Elec->FindBin(GetMass()+1e-9));	//1e-9 to prevent bin error
			hPoint->Scale(K0E);
		}
		else if (GetMass() > M_Kaon0 - M_Pion - M_Electron)
			hPoint->Scale(0);	//simple scaling

		hTotal->Add(hPoint);
	}

	//mu+ -> nu_mu_bar e+ nu_e
	if (hPoint = sxFlux->Get(Hist::Muon))
	{
		if (Kine)
		{
			double ME = hMuonElec->GetBinContent(hMuonElec->FindBin(GetMass()+1e-9));	//1e-9 to prevent bin error
			hPoint->Scale(ME);
		}
		else if (GetMass() > M_Muon - M_Electron)
			hPoint->Scale(0);

		hTotal->Add(hPoint);
	}

	//Ds+ -> e+ nu_e
	if (hPoint = sxFlux->Get(Hist::Charm))
	{
		hPoint->Scale(Kine::Unhelicity(M_Charm, M_Electron, GetMass(), GetHelicity()));

		hTotal->Add(hPoint);
	}
}

//Make muonic components, requires Neutrino (T if neutrino, F if antineutrino), the Flux object, the mass
void FluxDriver::MakeMuonComponent(Flux *sxFlux)
{
	TH1D *hPoint, *hTotal;

	//pi+ -> mu+ nu_mu
	if (hPoint = sxFlux->Get(Hist::Pion))
	{	
		hPoint->Scale(Kine::Unhelicity(M_Pion, M_Muon, GetMass(), GetHelicity()));

		hTotal->Add(hPoint);
	}

	//K+ -> mu+ nu_mu	(63.56%)
	//K+ -> pi0 mu+ nu_mu	(3.53%)
	if (hPoint = sxFlux->Get(Hist::Kaon))
	{
		double KaonFactor = 63.56/(63.56+3.35) * Kine::Unhelicity(M_Kaon, M_Muon, GetMass(), GetHelicity());	//Two body
		if (Kine)
		{
			double KM = hKaonMuon->GetBinContent(hKaonMuon->FindBin(GetMass()+1e-9)+1);	//1e-9 to prevent bin error
			KaonFactor += 3.35/(63.56+3.35) * KM;	//Three body
		}
		else if (GetMass() < M_Kaon - M_Pion0 - M_Muon)
			KaonFactor += 3.35/(63.56+3.35);	//simple scaling
		else KaonFactor += 0;	//simple scaling
		hPoint->Scale(KaonFactor);

		hTotal->Add(hPoint);
	}

	//K0 -> pi- mu+ nu_mu
	if (hPoint = sxFlux->Get(Hist::Kaon0))
	{
		if (Kine)
		{
			double K0M = hKaon0Muon->GetBinContent(hKaon0Muon->FindBin(GetMass()+1e-9));	//1e-9 to prevent bin error
				hPoint->Scale(K0M);	//Three body
		}
		else if (GetMass() > M_Kaon0 - M_Pion - M_Muon)
			hPoint->Scale(0);	//simple scaling

		hTotal->Add(hPoint);
	}

	//mu- -> nu_mu e- nu_e_bar
	if (hPoint = sxFlux->Get(Hist::Muon))
	{
		if (Kine)
		{
			double MM = hMuonMuon->GetBinContent(hMuonMuon->FindBin(GetMass()+1e-9));	//1e-9 to prevent bin error
			hPoint->Scale(MM);	//Three body
		}
		else if (GetMass() > M_Muon - M_Electron)
			hPoint->Scale(0);

		hTotal->Add(hPoint);
	}

	//Ds+ -> mu+ nu_mu
	if (hPoint = sxFlux->Get(Hist::Charm))
	{
		hPoint->Scale(Kine::Unhelicity(M_Charm, M_Muon, GetMass(), GetHelicity()));

		hTotal->Add(hPoint);
	}	
}

//Make tauonic components, requires Neutrino (T if neutrino, F if antineutrino), the Flux object, the mass
void FluxDriver::MakeTauComponent(Flux *sxFlux)
{
	TH1D *hPoint, *hTotal;

	//Ds+ -> tau+ nu_tau
	if (hPoint = sxFlux->Get(Hist::Charm))
	{
		//needs special scaling
		hPoint->Scale(Kine::Unhelicity(M_Charm, M_Tau, GetMass(), GetHelicity()));

		//"manual" modifier from empirical observation
		//the histogram is stretched and pulled to match the MC simulation spectrum
		//
		if (Mod && GetMass() < M_Charm - M_Tau)
		{
			TH1D *hTemp = dynamic_cast<TH1D*> (hPoint->Clone());
			hPoint->Reset("ICES");

			double xdir, ydir;
			Modify(xdir ,ydir, GetMass());
			double EnStep = (RangeEnd-RangeStart)/5000.0;
			for (double Energy = RangeStart; Energy < RangeEnd; Energy += EnStep)
			{
				double Flux = hTemp->GetBinContent(hTemp->FindBin(Energy));
				hPoint->Fill(Energy*xdir, Flux*BinNumber/5000.0);	//fix end point
			}

			hPoint->Scale(xdir);		//fix peak
			hPoint->Scale(ydir);		//fix peak
			hTemp->Delete();
		}

		hTotal->Add(hPoint);
	}

	//tau+ -> pi+ nu_tau
	if (hPoint = sxFlux->Get(Hist::Pion))
	{
		hPoint->Scale(Kine::Unhelicity(M_Tau, M_Pion, GetMass(), GetHelicity()));

		hTotal->Add(hPoint);
	}

	//tau+ -> pi+ pi0 nu_tau
	if (hPoint = sxFlux->Get(Hist::2Pion))
	{
		if (false)	//not implemented
		{
			//double P2 = hTau2Pion->GetBinContent(hTau2Pion->FindBin(GetMass()+1e-9));	//1e-9 to prevent bin error
			//hPoint->Scale(P2);	//Three body
		}
		else if (GetMass() > M_Tau - M_Pion - M_Pion0)	//only hard cut threshold
			hPoint->Scale(0);

		hTotal->Add(hPoint);
	}

	//tau+ -> nu_tau_bar e+ nu_e
	if (hPoint = sxFlux->Get(Hist::TauE))
	{
		if (Kine)
		{
			double TE = hTauETau->GetBinContent(hTauETau->FindBin(GetMass()+1e-9));	//1e-9 to prevent bin error
			hPoint->Scale(TE);	//Three body
		}
		else if (GetMass() > M_Tau - M_Electron)
			hPoint->Scale(0);

		hTotal->Add(hPoint);
	}

	//tau+ -> nu_tau_bar mu+ nu_mu
	if (hPoint = sxFlux->Get(Hist::TauM))
	{
		if (Kine)
		{
			double TM = hTauMTau->GetBinContent(hTauMTau->FindBin(GetMass()+1e-9));	//1e-9 to prevent bin error
			hPoint->Scale(TM);	//Three body
		}
		else if (GetMass() > M_Tau - M_Muon)
			hPoint->Scale(0);

		hTotal->Add(hPoint);
	}
}

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

//return the intensity of the flux at given energy
//
double FluxDriver::GetIntensity(Neutrino *N)
{
	return GetIntensity(N->EnergyKin());
}

double FluxDriver::GetIntensity(double Energy)	//Return flux intensity, given energy, simple linear interpolation
{									//NvA == 1 -> Neutrino, NvA == 0 -> Antineutrino
	double Ue, Um, Ut;						//Uu == 1 -> mixing , Uu == 0 -> maximal (light neutrino)

	if (Uu)
	{
		Ue = GetUe();
		Um = GetUm();
		Ut = GetUt();
	}
	else
	{
		Ue = 1.0;
		Um = 1.0;
		Ut = 1.0;
	}

	TH1D *hTotElec, *hTotMuon, *hTotTau;

	if (NvA)
	{
		hTotElec = hTotalEn;
		hTotMuon = hTotalMn;
		hTotTau  = hTotalTn;
	}
	else
	{
		hTotElec = hTotalEa;
		hTotMuon = hTotalMa;
		hTotTau  = hTotalTa;
	}

	int Bin;
	double I1, I2, E1, E2;

	Bin = hTotElec->FindBin(Energy);
	I1 = hTotElec->GetBinContent(Bin);
	E1 = hTotElec->GetBinCenter(Bin);
	if (Energy < hTotElec->GetBinCenter(Bin))
	{
		I2 = hTotElec->GetBinContent(Bin-1);
		E2 = hTotElec->GetBinCenter(Bin-1);
	}
	else
	{
		I2 = hTotElec->GetBinContent(Bin+1);
		E2 = hTotElec->GetBinCenter(Bin+1);
	}

	double IntElec = (Energy-E1)*(I2-I1)/(E2-E1) + I1;

	Bin = hTotMuon->FindBin(Energy);
	I1 = hTotMuon->GetBinContent(Bin);
	E1 = hTotMuon->GetBinCenter(Bin);
	if (Energy < hTotMuon->GetBinCenter(Bin))
	{
		I2 = hTotMuon->GetBinContent(Bin-1);
		E2 = hTotMuon->GetBinCenter(Bin-1);
	}
	else
	{
		I2 = hTotMuon->GetBinContent(Bin+1);
		E2 = hTotMuon->GetBinCenter(Bin+1);
	}

	double IntMuon = (Energy-E1)*(I2-I1)/(E2-E1) + I1;

	Bin = hTotTau->FindBin(Energy);
	I1 = hTotTau->GetBinContent(Bin);
	E1 = hTotTau->GetBinCenter(Bin);
	if (Energy < hTotTau->GetBinCenter(Bin))
	{
		I2 = hTotTau->GetBinContent(Bin-1);
		E2 = hTotTau->GetBinCenter(Bin-1);
	}
	else
	{
		I2 = hTotTau->GetBinContent(Bin+1);
		E2 = hTotTau->GetBinCenter(Bin+1);
	}

	double IntTau = (Energy-E1)*(I2-I1)/(E2-E1) + I1;

	return Ue*Ue*IntElec + Um*Um*IntMuon + Ut*Ut*IntTau;
} 

double FluxDriver::InterpolateIntensity(TH1D* Hist, double Energy)
{
	int Bin = Hist->FindBin(Energy);
	double I1 = Hist->GetBinContent(Bin);
	double E1 = Hist->GetBinCenter(Bin);
	double I2, E2;

	double Ret = 0.0;
	if (Bin > 2 && Bin < GetBinNumber())
	{
		if (Energy < Hist->GetBinCenter(Bin) && )
		{
			I2 = Hist->GetBinContent(Bin-1);
			E2 = Hist->GetBinCenter(Bin-1);
		}
		else
		{
			I2 = Hist->GetBinContent(Bin+1);
			E2 = Hist->GetBinCenter(Bin+1);
		}

		return (Energy-E1)*(I2-I1)/(E2-E1) + I1;
	}
	else
		return 0.0;
}

void FluxDriver::SetBaseline(double Baseline)
{
	ScaleAll(1.0/(Baseline*Baseline));
}

void FluxDriver::SetPOT(double POT)
{
	ScaleAll(POT);
}

void FluxDriver::SetArea(double Area)
{
	ScaleAll(Area);
}

void FluxDriver::ScaleAll(double X)
{
	hTotalE->Scale(X);
	hPionE->Scale(X);
	hKaonE->Scale(X);
	hKaon0E->Scale(X);
	hMuonE->Scale(X);
	hCharmE->Scale(X);

	hTotalM->Scale(X);
	hPionM->Scale(X);
	hKaonM->Scale(X);
	hKaon0M->Scale(X);
	hMuonM->Scale(X);
	hCharmM->Scale(X);

	hTotalT->Scale(X);
	hPionT->Scale(X);
	h2PionT->Scale(X);
	hCharmT->Scale(X);
	hTauET->Scale(X);
	hTauMT->Scale(X);
}

bool FluxDriver::IsChanged(Neutrino *N)
{
	bool Ret = (fabs(N->GetMass() - GetMass()) > 1e-9) ||
	           (N->GetHelicity() != GetHelicity()) ||
		   (N->IsParticle() != IsParticle());

	SetMass(N->GetMass());
	SetHelicity(N->GetHelicity());
	SetParticle(N->IsParticle());

	return Ret;
}

void FluxDriver::SetMass(double Mass)
{
	fMass = Mass;
}

void FluxDriver::SetHelicity(int Helicity)
{
	iHel = Helicity;
}

void FluxDriver::SetParticle(bool Particle)
{
	bParticle = Particle;
}

void FluxDriver::SetMixings(double *Mixings)
{
	fUe = Mixings[0];
	fUm = Mixings[1];
	fUt = Mixings[2];
}

double FluxDriver::GetMass()
{
	return fMass;
}

int FluxDriver::GetHelicity()
{
	return iHel;
}

bool FluxDriver::IsParticle()
{
	return bParticle;
}

double FluxDriver::GetUe()
{
	return fUe;
}

double FluxDriver::GetUm()
{
	return fUm;
}

double FluxDriver::GetUt()
{
	return fUt;
}

//Modificator for charm to tau flux
void FluxDriver::Modify(double &xdir, double &ydir, double M_Sterile)
{
	for (unsigned int i = 0; i < vMdir.size(); ++i)
	{
		if (vMdir.at(i) > M_Sterile || fabs(M_Sterile-vMdir.at(i)) < 1e-9)
		{
			xdir = vXdir.at(i);
			ydir = vYdir.at(i);
			break;
		}
	}
}
