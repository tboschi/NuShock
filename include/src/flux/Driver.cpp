#include "Driver.h"

Driver::Driver(std::string fluxConfig, bool FHC) :
	fxNuElectron(0),
	fxNuMuon(0),
	fxNuTau(0),
	fxHeavyElectron(0),
	fxHeavyMuon(0),
	fxHeavyTau(0)
{
	CardDealer fc(fluxConfig);//, true);

	std::string name;
	if (FHC)
	{
		//std::cout << "D is FHC" << std::endl;
		if (fc.Get("Electron_", name))
		{
			fxNuElectron = new Flux(name);
			//std::cout << "E\tflux name " << name << std::endl;
		}
		if (fc.Get("Muon_", name))
		{
			//std::cout << "M\tflux name " << name << std::endl;
			fxNuMuon = new Flux(name);
		}
		if (fc.Get("Tau_", name))
		{
			//std::cout << "T\tflux name " << name << std::endl;
			fxNuTau = new Flux(name);
		}
	}
	else
	{
		//std::cout << "D is RHC" << std::endl;
		if (fc.Get("ElectronBar_", name))
			fxNuElectron = new Flux(name);
		if (fc.Get("MuonBar_", name))
			fxNuMuon = new Flux(name);
		if (fc.Get("TauBar_", name))
			fxNuTau = new Flux(name);
	}

	std::map<std::string, std::string> mMod;
	std::map<std::string, std::string>::iterator id;
	if (fc.Get("Mod_", mMod))
	{
		for (id = mMod.begin(); id != mMod.end(); ++id)
		{
			std::vector<double> *vM, *vA, *vB, *vP;
			if (id->first  == "Charm")
			{
				vM = &vM_Charm;
				vA = &vA_Charm;
				vB = &vB_Charm;
				vP = &vP_Charm;
			}
			else if (id->first == "TauE")
			{
				vM = &vM_TauE;
				vA = &vA_TauE;
				vB = &vB_TauE;
				vP = &vP_TauE;
			}
			else if (id->first == "TauM")
			{
				vM = &vM_TauM;
				vA = &vA_TauM;
				vB = &vB_TauM;
				vP = &vP_TauM;
			}
			else if (id->first == "PPion")
			{
				vM = &vM_PPion;
				vA = &vA_PPion;
				vB = &vB_PPion;
				vP = &vP_PPion;
			}
			else if (id->first == "Pion")
			{
				vM = &vM_Pion;
				vA = &vA_Pion;
				vB = &vB_Pion;
				vP = &vP_Pion;
			}

			std::ifstream inMod(id->second);
			std::string line;
			double mass, A, B, P;
			while (getline(inMod, line))
			{
				if (line[0] == '#')
					continue;

				std::stringstream ssl(line);
				ssl >> mass >> A >> B >> P;

				vM->push_back(mass);
				vA->push_back(A);
				vB->push_back(B);
				vP->push_back(P);
			}
			inMod.close();
		}
	}

	Mass_prev = -1.0;
}

//deconstructor
Driver::~Driver()
{
	delete fxNuElectron;
	delete fxNuMuon;
	delete fxNuTau;

	delete fxHeavyElectron;
	delete fxHeavyMuon;
	delete fxHeavyTau;
}

//clones the histograms in the kinematic files
void Driver::CloneCopy(TH1D*& T, TObject* X)
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
bool Driver::MakeFlux(Neutrino &N)
{
	if (!IsChanged(N))	//compute only if particle is changed
	{
		std::cerr << "WARNING: recopmuting flux with same neutrino. "
			  << "This should not be happening!" << std::endl;
		return false;
	}
	else
	{
		//std::cout << "Making flux in driver" << std::endl;
		delete fxHeavyElectron;
		delete fxHeavyMuon;
		delete fxHeavyTau;
		
		fxHeavyElectron = 0;
                fxHeavyMuon = 0;
                fxHeavyTau = 0;

		//std::cout << "mixings " << N.Ue() << ", " << N.Um() << ", " << N.Ut() << std::endl;
		if (fxNuElectron)// && N.Ue() > 0)
		{
			MakeElecComponent(fxHeavyElectron = new Flux(*fxNuElectron), N);
			//std::cout << "making E component " << fxHeavyElectron << std::endl;
		}

		if (fxNuMuon)// && N.Um() > 0)
		{
			MakeMuonComponent(fxHeavyMuon = new Flux(*fxNuMuon), N);
			//std::cout << "making M component " << fxHeavyMuon << std::endl;
		}

		if (fxNuTau)// && N.Ut() > 0)
		{
			MakeTauComponent(fxHeavyTau = new Flux(*fxNuTau), N);
			//std::cout << "making T component " << fxHeavyTau << std::endl;
		}

		return true;
	}
}

//Make electronic components, requires Neutrino (T if neutrino, F if antineutrino), the Flux object, the mass
void Driver::MakeElecComponent(Flux *fxFlux, Neutrino &N)
{
	//pi+ -> e+ nu_e
	fxFlux->Scale(Flux::Pion, N.ProductionScale("PionE"));

	//K+ -> pi0 e+ nu_e	(5.07 %)	//I just assume that both decays are equally probable 
	//K+ -> e+ nu_e		(1.582e-3)	//for each different energy, so it is a linear combination
	double KaonFactor = 1.582e-3/(1.582e-3+5.07) * N.ProductionScale("KaonE") +
			    5.07/(1.582e-3+5.07) * N.ProductionScale("KaonCE");
	fxFlux->Scale(Flux::Kaon, KaonFactor);

	//K0 -> pi+ e+ nu_e
	fxFlux->Scale(Flux::Kaon0, N.ProductionScale("Kaon0E"));

	//mu+ -> nu_mu_bar e+ nu_e
	fxFlux->Scale(Flux::Muon, N.ProductionScale("MuonE"));

	//Ds+ -> e+ nu_e
	fxFlux->Scale(Flux::Charm, N.ProductionScale("CharmE"));

	fxFlux->Add();
}

//Make muonic components, requires Neutrino (T if neutrino, F if antineutrino), the Flux object, the mass
void Driver::MakeMuonComponent(Flux *fxFlux, Neutrino &N)
{
	//pi+ -> mu+ nu_mu
	fxFlux->Scale(Flux::Pion, N.ProductionScale("PionM"));

	//K+ -> mu+ nu_mu	(63.56%)
	//K+ -> pi0 mu+ nu_mu	(3.53%)
	double KaonFactor = 63.56/(63.56+3.35) * N.ProductionScale("KaonM") +
			    3.35/(63.56+3.35) * N.ProductionScale("KaonCM");
	fxFlux->Scale(Flux::Kaon, KaonFactor);

	//K0 -> pi- mu+ nu_mu
	fxFlux->Scale(Flux::Kaon0, N.ProductionScale("Kaon0M"));

	//mu- -> nu_mu e- nu_e_bar
	fxFlux->Scale(Flux::Muon, N.ProductionScale("MuonM"));

	//Ds+ -> mu+ nu_mu
	fxFlux->Scale(Flux::Charm, N.ProductionScale("CharmM"));

	fxFlux->Add();
}

//Make tauonic components, requires Neutrino (T if neutrino, F if antineutrino), the Flux object, the mass
void Driver::MakeTauComponent(Flux *fxFlux, Neutrino &N)
{
	//Ds+ -> tau+ nu_tau
	//if (vM_Charm.size())
	if (false)
	{
		double Sx, Ex;
		double Mod = Modify(Flux::Charm, Sx, Ex, N.Mass());

		if (fxFlux->Stretch(Flux::Charm, Sx, Ex))
			fxFlux->Scale(Flux::Charm, Mod);
		else
			fxFlux->Scale(Flux::Charm, 0.0);
	}
	fxFlux->Scale(Flux::Charm, N.ProductionScale("CharmT"));

	//tau+ -> pi+ nu_tau
	//if (vM_Pion.size())
	if (false)
	{
		double Sx, Ex;
		double Mod = Modify(Flux::Pion, Sx, Ex, N.Mass());

		if (fxFlux->Stretch(Flux::Pion, Sx, Ex))
			fxFlux->Scale(Flux::Pion, Mod);
		else
			fxFlux->Scale(Flux::Pion, 0.0);
	}
	fxFlux->Scale(Flux::Pion, N.ProductionScale("TauPI"));

	//tau+ -> pi+ pi0 nu_tau	//crossing simmetries
	//if (vM_PPion.size())
	if (false)
	{
		double Sx, Ex;
		double Mod = Modify(Flux::PPion, Sx, Ex, N.Mass());

		if (fxFlux->Stretch(Flux::PPion, Sx, Ex))
			fxFlux->Scale(Flux::PPion, Mod);
		else
			fxFlux->Scale(Flux::PPion, 0.0);
	}
	fxFlux->Scale(Flux::PPion, N.ProductionScale("Tau2PI"));

	//tau+ -> nu_tau_bar e+ nu_e
	//if (vM_TauE.size())
	if (false)
	{
		double Sx, Ex;
		double Mod = Modify(Flux::TauE, Sx, Ex, N.Mass());

		if (fxFlux->Stretch(Flux::TauE, Sx, Ex))
			fxFlux->Scale(Flux::TauE, Mod);
		else
			fxFlux->Scale(Flux::TauE, 0.0);
	}
	fxFlux->Scale(Flux::TauE, N.ProductionScale("TauET"));

	//tau+ -> nu_tau_bar mu+ nu_mu
	//if (vM_TauM.size())
	if (false)
	{
		double Sx, Ex;
		double Mod = Modify(Flux::TauM, Sx, Ex, N.Mass());

		if (fxFlux->Stretch(Flux::TauM, Sx, Ex))
			fxFlux->Scale(Flux::TauM, Mod);
		else
			fxFlux->Scale(Flux::TauM, 0.0);
	}
	fxFlux->Scale(Flux::TauM, N.ProductionScale("TauMT"));	//Three body

	fxFlux->Add();
}


//return the intensity of the flux at given energy
//

double Driver::Intensity(Neutrino &N)	//Return flux intensity, given energy, simple linear interpolation
{
	double Energy = N.EnergyKin();

	//std::cout << fxHeavyElectron << ", " << fxHeavyMuon << ", " << fxHeavyTau << std::endl;
	double Intensity = 0;
	if (fxHeavyElectron)
		Intensity += N.Ue(2) * InterpolateIntensity(fxHeavyElectron->Get(Flux::Total), Energy);
	if (fxHeavyMuon)
		Intensity += N.Um(2) * InterpolateIntensity(fxHeavyMuon->Get(Flux::Total), Energy);
	if (fxHeavyTau)
		Intensity += N.Ut(2) * InterpolateIntensity(fxHeavyTau->Get(Flux::Total), Energy);

	return Intensity;
}

double Driver::InterpolateIntensity(TH1D* Hist, double Energy)
{
	if (!Hist)
		return 0.0;
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

void Driver::Scale(double X)
{
	if(fxHeavyElectron)
		fxHeavyElectron->Scale(Flux::Total, X);
	if(fxHeavyMuon)
		fxHeavyMuon->Scale(Flux::Total, X);
	if(fxHeavyTau)
		fxHeavyTau->Scale(Flux::Total, X);
}

bool Driver::IsChanged(Neutrino &N)
{
	bool Ret = (fabs(N.Mass() - Mass_prev) > 1e-9) ||
	           (N.Helicity() != Helicity_prev) ||
		   (N.IsParticle() != Particle_prev);

	if (Ret)
	{
		Mass_prev = N.Mass();
		Helicity_prev = N.Helicity();
		Particle_prev = N.IsParticle();
	}

	return Ret;
}

//Modificator for charm to tau flux
double Driver::Modify(Flux::Hist Name, double &A, double &B, double Mass)
{
	std::vector<double> *vM, *vA, *vB, *vP;
	switch (Name)
	{
		case Flux::Charm:
			vM = &vM_Charm;
			vA = &vA_Charm;
			vB = &vB_Charm;
			vP = &vP_Charm;
			break;
		case Flux::TauE:
			vM = &vM_TauE;
			vA = &vA_TauE;
			vB = &vB_TauE;
			vP = &vP_TauE;
			break;
		case Flux::TauM:
			vM = &vM_TauM;
			vA = &vA_TauM;
			vB = &vB_TauM;
			vP = &vP_TauM;
			break;
		case Flux::Pion:
			vM = &vM_Pion;
			vA = &vA_Pion;
			vB = &vB_Pion;
			vP = &vP_Pion;
			break;
		case Flux::PPion:
			vM = &vM_PPion;
			vA = &vA_PPion;
			vB = &vB_PPion;
			vP = &vP_PPion;
			break;
	}

	for (unsigned int i = 0; i < vM->size(); ++i)
	{
		if (vM->at(i) > Mass || fabs(Mass - vM->at(i)) < 1e-9)
		{
			A = vA->at(i);
			B = vB->at(i);
			return 1.0/vP->at(i);
		}
	}

	return 0.0;
}

double Driver::RangeWidth()
{
	double Start, End;
	return RangeWidth(Start, End);
}

double Driver::RangeWidth(double &Start, double &End)
{
	return Range(Start, End) / BinNumber();
}

double Driver::Range()
{
	double Start, End;
	return Range(Start, End);
}

double Driver::Range(double &Start, double &End)
{
	Start = RangeStart();
	End = RangeEnd();
	return End - Start;
}

double Driver::RangeStart()
{
	if (fxNuElectron)
		return fxNuElectron->RangeStart();
	else if (fxNuMuon)
		return fxNuMuon->RangeStart();
	else if (fxNuTau)
		return fxNuTau->RangeStart();
	else
		return -1.0;
}

double Driver::RangeEnd()
{
	if (fxNuElectron)
		return fxNuElectron->RangeEnd();
	else if (fxNuMuon)
		return fxNuMuon->RangeEnd();
	else if (fxNuTau)
		return fxNuTau->RangeEnd();
	else
		return -1.0;
}

int Driver::BinNumber()
{
	if (fxNuElectron)
		return fxNuElectron->BinNumber();
	else if (fxNuMuon)
		return fxNuMuon->BinNumber();
	else if (fxNuTau)
		return fxNuTau->BinNumber();
	else
		return -1.0;
}
