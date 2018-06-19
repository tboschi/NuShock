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
	fxNuMuon = 0;
	fxNuTau = 0;

	fxHeavyElectron = 0;
	fxHeavyMuon = 0;
	fxHeavyTau = 0;

	std::stringstream ssL;
	std::string Line, Key, Name;
	
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
			if (Key.find("Muon_")	     != std::string::npos) fxNuMuon = new Flux(Name);
			if (Key.find("Tau_")	     != std::string::npos) fxNuTau = new Flux(Name);

			//histogram details (number of bins, start of range, end of range)
			//
			//if (Key.find("BinNumber")    != std::string::npos) BinNumber = std::strtoul(Name.c_str(), NULL, 10);
			//if (Key.find("RangeStart")   != std::string::npos) RangeStart = std::strtod(Name.c_str(), NULL);
			//if (Key.find("RangeEnd")     != std::string::npos) RangeEnd = std::strtod(Name.c_str(), NULL);

			//phase space functions for 3body decays are loaded here
			//hopefully this is not needed
			/*
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
			*/

			//correction for charm to tau decay
			//
			if (Key.find("Modifier") != std::string::npos)
			{
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

		if (fxNuElectron)
			MakeElecComponent(fxHeavyElectron = new Flux(*fxNuElectron), N);

		if (fxNuMuon)
			MakeMuonComponent(fxHeavyMuon = new Flux(*fxNuMuon), N);

		if (fxNuTau)
			MakeTauComponent(fxHeavyTau = new Flux(*fxNuTau), N);

		return true;
	}
}

//Make electronic components, requires Neutrino (T if neutrino, F if antineutrino), the Flux object, the mass
void FluxDriver::MakeElecComponent(Flux *fxFlux, Neutrino *N)
{
	TH1D *hPoint, *hTotal;

	//pi+ -> e+ nu_e
	fxFlux->Scale(Flux::Pion, N->ProductionScale("PionE"));

	//K+ -> pi0 e+ nu_e	(5.07 %)	//I just assume that both decays are equally probable 
	//K+ -> e+ nu_e		(1.582e-3)	//for each different energy, so it is a linear combination
	double KaonFactor = 1.582e-3/(1.582e-3+5.07) * N->ProductionScale("KaonE") +
			    5.07/(1.582e-3+5.07) * N->ProductionScale("KaonCE");
	fxFlux->Scale(Flux::Kaon, KaonFactor);

	//K0 -> pi+ e+ nu_e
	fxFlux->Scale(Flux::Kaon0, N->ProductionScale("Kaon0E"));

	//mu+ -> nu_mu_bar e+ nu_e
	fxFlux->Scale(Flux::Muon, N->ProductionScale("MuonE"));

	//Ds+ -> e+ nu_e
	fxFlux->Scale(Flux::Charm, N->ProductionScale("CharmE"));

	fxFlux->Add();
}

//Make muonic components, requires Neutrino (T if neutrino, F if antineutrino), the Flux object, the mass
void FluxDriver::MakeMuonComponent(Flux *fxFlux, Neutrino* N)
{
	//pi+ -> mu+ nu_mu
	fxFlux->Scale(Flux::Pion, N->ProductionScale("PionM"));

	//K+ -> mu+ nu_mu	(63.56%)
	//K+ -> pi0 mu+ nu_mu	(3.53%)
	double KaonFactor = 63.56/(63.56+3.35) * N->ProductionScale("KaonM") +
			    3.35/(63.56+3.35) * N->ProductionScale("KaonCM");
	fxFlux->Scale(Flux::Kaon, KaonFactor);

	//K0 -> pi- mu+ nu_mu
	fxFlux->Scale(Flux::Kaon0, N->ProductionScale("Kaon0M"));

	//mu- -> nu_mu e- nu_e_bar
	fxFlux->Scale(Flux::Muon, N->ProductionScale("MuonM"));

	//Ds+ -> mu+ nu_mu
	fxFlux->Scale(Flux::Charm, N->ProductionScale("CharmM"));

	fxFlux->Add();
}

//Make tauonic components, requires Neutrino (T if neutrino, F if antineutrino), the Flux object, the mass
void FluxDriver::MakeTauComponent(Flux *fxFlux, Neutrino *N)
{
	//Ds+ -> tau+ nu_tau
	//"manual" modifier from empirical observation
	//the histogram is stretched and pulled to match the MC simulation spectrum
	double Modifier = 0;
	TH1D *hPoint;
	if ((hPoint = fxFlux->Get(Flux::Charm)) && 
	    vMdir.size() && 
	    (N->Mass() < Const::fMCharm - Const::fMTau) );
	{
		TH1D *hTemp = dynamic_cast<TH1D*> (hPoint->Clone());
		hPoint->Reset("ICES");

		double xdir, ydir;
		Modifier = Modify(xdir ,ydir, N->Mass());

		double Start, End;
		double EnStep = RangeBin(Start, End)/50;
		for (double Energy = Start; Energy < End; Energy += EnStep)
		{
			double Flux = hTemp->GetBinContent(hTemp->FindBin(Energy));
			hPoint->Fill(Energy*xdir, Flux/50.0);	//fix end point
		}

		delete hTemp;
	}
	fxFlux->Scale(Flux::Charm, N->ProductionScale("CharmT") * Modifier);

	//tau+ -> pi+ nu_tau
	fxFlux->Scale(Flux::Pion, N->ProductionScale("TauPion"));

	//tau+ -> pi+ pi0 nu_tau	//crossing simmetries
	fxFlux->Scale(Flux::PPion, N->Mass() > Const::fMTau - Const::fMPion - Const::fMPion0 ? 0.0 : 1.0);

	//tau+ -> nu_tau_bar e+ nu_e
	fxFlux->Scale(Flux::TauE, N->ProductionScale("TauEE"));

	//tau+ -> nu_tau_bar mu+ nu_mu
	fxFlux->Scale(Flux::TauM, N->ProductionScale("TauMT"));	//Three body

	fxFlux->Add();
}


//return the intensity of the flux at given energy
//

double FluxDriver::Intensity(Neutrino *N)	//Return flux intensity, given energy, simple linear interpolation
{
	double Energy = N->EnergyKin();

	double Intensity = 0;
	if (fxHeavyElectron)
		Intensity += N->Ue(2) * InterpolateIntensity(fxHeavyElectron->Get(Flux::Total), Energy);
	if (fxHeavyMuon)
		Intensity += N->Ue(2) * InterpolateIntensity(fxHeavyMuon->Get(Flux::Total), Energy);
	if (fxHeavyTau)
		Intensity += N->Ue(2) * InterpolateIntensity(fxHeavyTau->Get(Flux::Total), Energy);

	return Intensity;
}

double FluxDriver::InterpolateIntensity(TH1D* Hist, double Energy)
{
	int Bin = Hist->FindBin(Energy);
	double I1 = Hist->GetBinContent(Bin);
	double E1 = Hist->GetBinCenter(Bin);
	double I2, E2;

	double Ret = 0.0;
	if (Bin > 1 && Bin < BinNumber())
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
	if(fxHeavyElectron)
		fxHeavyElectron->Scale(Flux::Total, X);
	if(fxHeavyMuon)
		fxHeavyMuon->Scale(Flux::Total, X);
	if(fxHeavyTau)
		fxHeavyTau->Scale(Flux::Total, X);
}

bool FluxDriver::IsChanged(Neutrino *N)
{
	bool Ret = (fabs(N->Mass() - Mass_prev) > 1e-9) ||
	           (N->Helicity() != Helicity_prev) ||
		   (N->IsParticle() != Particle_prev);

	return Ret;
}

//Modificator for charm to tau flux
double FluxDriver::Modify(double &xdir, double &ydir, double M_Sterile)
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

	return xdir*ydir;
}

double FluxDriver::RangeBin()
{
	double Start, End;
	return RangeBin(Start, End);
}

double FluxDriver::RangeBin(double &Start, double &End)
{
	return Range(Start, End) / BinNumber();
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
	return End - Start;
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
