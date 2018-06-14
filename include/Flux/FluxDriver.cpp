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
	//fxNuElectronBar = 0;
	fxNuMuon = 0;
	//fxNuMuonBar = 0;
	fxNuTau = 0;
	//fxNuTauBar = 0;

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
			//if (Key.find("ElectronBar_") != std::string::npos) fxNuElectronBar = new Flux(Name);
			if (Key.find("Muon_")	     != std::string::npos) fxNuMuon = new Flux(Name);
			//if (Key.find("MuonBar_")     != std::string::npos) fxNuMuonBar = new Flux(Name);
			if (Key.find("Tau_")	     != std::string::npos) fxNuTau = new Flux(Name);
			//if (Key.find("TauBar_")	     != std::string::npos) fxNuTauBar = new Flux(Name);

			//histogram details (number of bins, start of range, end of range)
			//
			//if (Key.find("BinNumber")    != std::string::npos) BinNumber = std::strtoul(Name.c_str(), NULL, 10);
			//if (Key.find("RangeStart")   != std::string::npos) RangeStart = std::strtod(Name.c_str(), NULL);
			//if (Key.find("RangeEnd")     != std::string::npos) RangeEnd = std::strtod(Name.c_str(), NULL);

			//phase space functions for 3body decays are loaded here
			//hopefully this is not needed
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
	}
	ConfigFile.close();

	Mass_prev = -1.0;
}

//deconstructor
FluxDriver::~FluxDriver()
{
	delete fxNuElectron;
	//delete fxNuElectronBar;
	delete fxNuMuon;
	//delete fxNuMuonBar;
	delete fxNuTau;
	//delete fxNuTauBar;

	delete fxHeavyElectron;
	//delete fxHeavyElectronBar;
	delete fxHeavyMuon;
	//delete fxHeavyMuonBar;
	delete fxHeavyTau;
	//delete fxHeavyTauBar;
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
		delete fxHeavyElectron;
		delete fxHeavyMuon;
		delete fxHeavyTau;
		
		fxHeavyElectron = 0;
                fxHeavyMuon = 0;
                fxHeavyTau = 0;

		if (fxHeavyElectron = new Flux(*fxNuElectron))
			MakeElecComponent(fxHeavyElectron, N);

		if (fxHeavyMuon = new Flux(*fxNuMuon))
			MakeElecComponent(fxHeavyMuon, N);

		if (fxHeavyTau = new Flux(*fxNuTau))
			MakeElecComponent(fxHeavyTau, N);

		return true;
	}
}

//Make electronic components, requires Neutrino (T if neutrino, F if antineutrino), the Flux object, the mass
void FluxDriver::MakeElecComponent(Flux *fxFlux, Neutrino *N)
{
	TH1D *hPoint, *hTotal;

	//pi+ -> e+ nu_e
	if (hPoint = fxFlux->Get(Flux::Pion))
	{
		hPoint->Scale(N->ProductionScale(Amplitude::_PionE));
		//Kine::Unhelicity(M_Pion, M_Electron, GetMass(), GetHelicity()));

		hTotal->Add(hPoint);
	}

	//K+ -> pi0 e+ nu_e	(5.07 %)		//I just assume that both decays are equally probable 
	//K+ -> e+ nu_e		(1.582e-3)		//for each different energy, so it is a linear combination of the two
	if (hPoint = fxFlux->Get(Flux::Kaon))
	{
				    //branching ratios
		double KaonFactor = 1.582e-3/(1.582e-3+5.07) * N->ProductionScale(Amplitude::_KaonE);	//Two body
		if (Kine)
		{
			//double KE = hKaonElec->GetBinContent(hKaonElec->FindBin(GetMass()+1e-9)+1);	//1e-9 to prevent bin error
			KaonFactor += 5.07/(1.582e-3+5.07) * N->ProductionScale(Amplitude::_KaonCE);	//Three body
		}
		else if (N->Mass() < M_Kaon - M_Pion0 - M_Electron)
			KaonFactor += 5.07/(1.582e-3+5.07);	//simple scaling
		else KaonFactor += 0;	//simple scaling
		hPoint->Scale(KaonFactor);

		hTotal->Add(hPoint);
	}

	//K0 -> pi+ e+ nu_e
	if (hPoint = fxFlux->Get(Flux::Kaon0))
	{
		if (Kine)
		{
			//double K0E = hKaon0Elec->GetBinContent(hKaon0Elec->FindBin(GetMass()+1e-9));	//1e-9 to prevent bin error
			hPoint->Scale(N->ProductionScale(Amplitude::_Kaon0E));
		}
		else if (N->Mass() > M_Kaon0 - M_Pion - M_Electron)
			hPoint->Scale(0);	//simple scaling

		hTotal->Add(hPoint);
	}

	//mu+ -> nu_mu_bar e+ nu_e
	if (hPoint = fxFlux->Get(Flux::Muon))
	{
		if (Kine)
		{
			//double ME = hMuonElec->GetBinContent(hMuonElec->FindBin(GetMass()+1e-9));	//1e-9 to prevent bin error
			hPoint->Scale(N->ProductionScale(Amplitude::_MuonE));
		}
		else if (N->Mass() > M_Muon - M_Electron)
			hPoint->Scale(0);

		hTotal->Add(hPoint);
	}

	//Ds+ -> e+ nu_e
	if (hPoint = fxFlux->Get(Flux::Charm))
	{
		//hPoint->Scale(Kine::Unhelicity(M_Charm, M_Electron, Mass(), GetHelicity()));
		hPoint->Scale(N->ProductionScale(Amplitude::_CharmE));

		hTotal->Add(hPoint);
	}
}

//Make muonic components, requires Neutrino (T if neutrino, F if antineutrino), the Flux object, the mass
void FluxDriver::MakeMuonComponent(Flux *fxFlux, Neutrino* N)
{
	TH1D *hPoint, *hTotal;

	//pi+ -> mu+ nu_mu
	if (hPoint = fxFlux->Get(Flux::Pion))
	{	
		//hPoint->Scale(Kine::Unhelicity(M_Pion, M_Muon, Mass(), GetHelicity()));
		hPoint->Scale(N->ProductionScale(Amplitude::_PionM));

		hTotal->Add(hPoint);
	}

	//K+ -> mu+ nu_mu	(63.56%)
	//K+ -> pi0 mu+ nu_mu	(3.53%)
	if (hPoint = fxFlux->Get(Flux::Kaon))
	{
		//double KaonFactor = 63.56/(63.56+3.35) * Kine::Unhelicity(M_Kaon, M_Muon, Mass(), GetHelicity());	//Two body
		double KaonFactor = 63.56/(63.56+3.35) * N->ProductionScale(Amplitude::_KaonM);
		if (Kine)
		{
			//double KM = hKaonMuon->GetBinContent(hKaonMuon->FindBin(Mass()+1e-9)+1);	//1e-9 to prevent bin error
			//KaonFactor += 3.35/(63.56+3.35) * KM;	//Three body
			KaonFactor += 3.35/(63.56+3.35) * N->ProductionScale(Amplitude::_KaonCM);	//Three body
		}
		else if (N->Mass() < M_Kaon - M_Pion0 - M_Muon)
			KaonFactor += 3.35/(63.56+3.35);	//simple scaling
		else KaonFactor += 0;	//simple scaling
		hPoint->Scale(KaonFactor);

		hTotal->Add(hPoint);
	}

	//K0 -> pi- mu+ nu_mu
	if (hPoint = fxFlux->Get(Flux::Kaon0))
	{
		if (Kine)
		{
			//double K0M = hKaon0Muon->GetBinContent(hKaon0Muon->FindBin(GetMass()+1e-9));	//1e-9 to prevent bin error
			//hPoint->Scale(K0M);	//Three body
			hPoint->Scale(N->ProductionScale(Amplitude::_Kaon0M));	//Three body
		}
		else if (N->Mass() > M_Kaon0 - M_Pion - M_Muon)
			hPoint->Scale(0);	//simple scaling

		hTotal->Add(hPoint);
	}

	//mu- -> nu_mu e- nu_e_bar
	if (hPoint = fxFlux->Get(Flux::Muon))
	{
		if (Kine)
		{
			double MM = hMuonMuon->GetBinContent(hMuonMuon->FindBin(N->Mass()+1e-9));	//1e-9 to prevent bin error
			hPoint->Scale(MM);	//Three body
			hPoint->Scale(N->ProductionScale(Amplitude::_MuonM));	//Three body
		}
		else if (N->Mass() > M_Muon - M_Electron)
			hPoint->Scale(0);

		hTotal->Add(hPoint);
	}

	//Ds+ -> mu+ nu_mu
	if (hPoint = fxFlux->Get(Flux::Charm))
	{
		//hPoint->Scale(Kine::Unhelicity(M_Charm, M_Muon, GetMass(), GetHelicity()));
		hPoint->Scale(N->ProductionScale(Amplitude::_CharmM));

		hTotal->Add(hPoint);
	}	
}

//Make tauonic components, requires Neutrino (T if neutrino, F if antineutrino), the Flux object, the mass
void FluxDriver::MakeTauComponent(Flux *fxFlux, Neutrino *N)
{
	TH1D *hPoint, *hTotal;

	//Ds+ -> tau+ nu_tau
	if (hPoint = fxFlux->Get(Flux::Charm))
	{
		//needs special scaling
		//hPoint->Scale(Kine::Unhelicity(M_Charm, M_Tau, GetMass(), GetHelicity()));
		hPoint->Scale(N->ProductionScale(Amplitude::_CharmT));

		//"manual" modifier from empirical observation
		//the histogram is stretched and pulled to match the MC simulation spectrum
		//
		if (Mod && N->Mass() < M_Charm - M_Tau)
		{
			TH1D *hTemp = dynamic_cast<TH1D*> (hPoint->Clone());
			hPoint->Reset("ICES");

			double xdir, ydir;
			Modify(xdir ,ydir, N->Mass());
			double EnStep = (RangeEnd() - RangeStart())/5000.0;
			for (double Energy = RangeStart(); Energy < RangeEnd(); Energy += EnStep)
			{
				double Flux = hTemp->GetBinContent(hTemp->FindBin(Energy));
				hPoint->Fill(Energy*xdir, Flux*BinNumber()/5000.0);	//fix end point
			}

			hPoint->Scale(xdir);		//fix peak
			hPoint->Scale(ydir);		//fix peak
			hTemp->Delete();
		}

		hTotal->Add(hPoint);
	}

	//tau+ -> pi+ nu_tau
	if (hPoint = fxFlux->Get(Flux::Pion))
	{
		//hPoint->Scale(Kine::Unhelicity(M_Tau, M_Pion, GetMass(), GetHelicity()));
		hPoint->Scale(N->ProductionScale(Amplitude::_TauPion));

		hTotal->Add(hPoint);
	}

	//tau+ -> pi+ pi0 nu_tau	//crossing simmetries
	if (hPoint = fxFlux->Get(Flux::PPion))
	{
		if (false)	//not implemented
		{
			//double P2 = hTauPPion->GetBinContent(hTauPPion->FindBin(GetMass()+1e-9));	//1e-9 to prevent bin error
			//hPoint->Scale(P2);	//Three body
		}
		else if (N->Mass() > M_Tau - M_Pion - M_Pion0)	//only hard cut threshold
			hPoint->Scale(0);

		hTotal->Add(hPoint);
	}

	//tau+ -> nu_tau_bar e+ nu_e
	if (hPoint = fxFlux->Get(Flux::TauE))
	{
		if (Kine)
		{
			//double TE = hTauETau->GetBinContent(hTauETau->FindBin(GetMass()+1e-9));	//1e-9 to prevent bin error
			//hPoint->Scale(TE);	//Three body
			hPoint->Scale(N->ProductionScale(Amplitude::_TauEE));	//Three body
		}
		else if (N->Mass() > M_Tau - M_Electron)
			hPoint->Scale(0);

		hTotal->Add(hPoint);
	}

	//tau+ -> nu_tau_bar mu+ nu_mu
	if (hPoint = fxFlux->Get(Flux::TauM))
	{
		if (Kine)
		{
			//double TM = hTauMTau->GetBinContent(hTauMTau->FindBin(GetMass()+1e-9));	//1e-9 to prevent bin error
			//hPoint->Scale(TM);	//Three body
			hPoint->Scale(N->ProductionScale(Amplitude::_TauMT));	//Three body
		}
		else if (N->Mass() > M_Tau - M_Muon)
			hPoint->Scale(0);

		hTotal->Add(hPoint);
	}
}

double FluxDriver::Range()
{
	double Start, End;
	return Range(Start, End);
}

double FluxDriver::Range(double &Start, double &End)
{
	Start = RangeStart();
	End = RangeEnd();
	return Start - End;
}

double FluxDriver::RangeStart()
{
	if (fxHeavyElectron)
		return fxHeavyElectron->RangeStart();
	else if (fxHeavyMuon)
		return fxHeavyMuon->RangeStart();
	else if (fxHeavyTau)
		return fxHeavyTau->RangeStart();
	else
		return -1.0;
}

double FluxDriver::RangeEnd()
{
	if (fxHeavyElectron)
		return fxHeavyElectron->RangeEnd();
	else if (fxHeavyMuon)
		return fxHeavyMuon->RangeEnd();
	else if (fxHeavyTau)
		return fxHeavyTau->RangeEnd();
	else
		return -1.0;
}

int FluxDriver::BinNumber()
{
	if (fxHeavyElectron)
		return fxHeavyElectron->BinNumber();
	else if (fxHeavyMuon)
		return fxHeavyMuon->BinNumber();
	else if (fxHeavyTau)
		return fxHeavyTau->BinNumber();
	else
		return -1.0;
}

//return the intensity of the flux at given energy
//

double FluxDriver::Intensity(Neutrino *N)	//Return flux intensity, given energy, simple linear interpolation
{
	double Energy = N->EnergyKin();
	return pow(N->Ue(), 2) * InterpolateIntensity(fxHeavyElectron->Get(Flux::Total), Energy) +
	       pow(N->Um(), 2) * InterpolateIntensity(fxHeavyMuon->Get(Flux::Total), Energy) +
	       pow(N->Ut(), 2) * InterpolateIntensity(fxHeavyTau->Get(Flux::Total), Energy);
}

double FluxDriver::InterpolateIntensity(TH1D* Hist, double Energy)
{
	int Bin = Hist->FindBin(Energy);
	double I1 = Hist->GetBinContent(Bin);
	double E1 = Hist->GetBinCenter(Bin);
	double I2, E2;

	double Ret = 0.0;
	if (Bin > 2 && Bin < BinNumber())
	{
		if (Energy < Hist->GetBinCenter(Bin))
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
	fxHeavyElectron->Scale(X);
	fxHeavyMuon->Scale(X);
	fxHeavyTau->Scale(X);
}

bool FluxDriver::IsChanged(Neutrino *N)
{
	bool Ret = (fabs(N->Mass() - Mass_prev) > 1e-9) ||
	           (N->Helicity() != Helicity_prev) ||
		   (N->IsParticle() != Particle_prev);

	return Ret;
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
