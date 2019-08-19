#include <iostream>
#include <fstream>
#include <vector>
#include <cstring>
#include <getopt.h>

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TRandom3.h"

#include "tools.h"
#include "physics.h"
#include "detector.h"
#include "flux.h"

void Usage(char* argv0);

int main(int argc, char** argv)
{
	const struct option longopts[] = 
	{
		{"detconfig", 	required_argument,	0, 'd'},
		{"fluxconfig", 	required_argument,	0, 'f'},
		{"channel", 	required_argument,	0, 'c'},
		{"mass", 	required_argument,	0, 'm'},
		{"number", 	required_argument,	0, 'n'},
		{"output", 	required_argument,	0, 'o'},
		{"left", 	no_argument,		0, 'L'},
		{"right", 	no_argument,		0, 'R'},
		{"dirac", 	no_argument,		0, 'D'},
		{"majorana", 	no_argument,		0, 'j'},
		{"help", 	no_argument,		0, 'h'},
		{0,	0, 	0,	0},
	};

	int index; 
	int iarg = 0;
	opterr = 1;
	
	std::string channel, outName, detConfig, fluxConfig;
	unsigned int nevent = 10000;
	double mass = 0.0;
	bool left = false, right = false;	//default unpolarised
	bool dirac = true;	//default unpolarised
	double ue = 0, um = 0, ut = 0, mix = 0;

	std::string Channel = "ALL", module;
	
	while((iarg = getopt_long(argc,argv, "d:l:f:m:n:c:o:LRDjE:M:T:x:h", longopts, &index)) != -1)	
	{
		switch(iarg)
		{
			case 'd':
				detConfig.assign(optarg);
				break;
			case 'l':
				module.assign(optarg);
				break;
			case 'f':
				fluxConfig.assign(optarg);
				break;
			case 'm':
				mass = strtod(optarg, NULL);
				break;
			case 'n':
				nevent = strtod(optarg, NULL);
				break;
			case 'c':
				channel.assign(optarg);
				break;
			case 'o':
				outName.assign(optarg);
				break;
			case 'L':
				left = true;
				right = false;
				break;
			case 'R':
				left = false;
				right = true;
				break;
			case 'D':
				dirac = true;
				break;
			case 'j':
				dirac = false;
				break;
			case 'E':
				ue = std::sqrt(std::strtod(optarg, NULL));
				break;
			case 'M':
				um = std::sqrt(std::strtod(optarg, NULL));
				break;
			case 'T':
				ut = std::sqrt(std::strtod(optarg, NULL));
				break;
			case 'X':
				mix = std::sqrt(std::strtod(optarg, NULL));
				break;
			case 'h':
				Usage(argv[0]);
				return 1;
			default:
				break;
		}
	}

	//std::ostream &Out = (OutFile.is_open()) ? OutFile : std::cout;
	
	Neutrino testNu(mass);
	double mass0 = testNu.DecayThreshold(channel);
	double massE = testNu.ProductionThreshold();
	if (mass0 > massE)
	{
		double tmp = massE;
		massE = mass0;
		mass0 = tmp;
	}

	if (mass < mass0 || mass > massE)
	{
		std::cerr << "Decay " << channel << " is not allowed for a mass of " << mass << " GeV" << std::endl;
		return 1;
	}

	Tracker *theBox = new Tracker(detConfig, module);
	Engine *theEngine = new Engine(fluxConfig);	//creating 1FHC and 1RHC fluxedrivers

	//creating neutrino objs and loading them in the Engine
	if (dirac)
	{
		Neutrino nu0_L(mass, Neutrino::Left  | Neutrino::Dirac );
		Neutrino nu0_R(mass, Neutrino::Right | Neutrino::Dirac );
		Neutrino nuB_L(mass, Neutrino::Left  | Neutrino::Dirac | Neutrino::Antiparticle);
		Neutrino nuB_R(mass, Neutrino::Right | Neutrino::Dirac | Neutrino::Antiparticle);

		nu0_L.SetDecayChannel(channel);
		nu0_R.SetDecayChannel(channel);
		nuB_L.SetDecayChannel(channel);
		nuB_R.SetDecayChannel(channel);

		//Binding will make a copy of neutrino objects
		theEngine->BindNeutrino("nu0_L", nu0_L, Engine::FHC);
		theEngine->BindNeutrino("nu0_R", nu0_R, Engine::FHC);
		theEngine->BindNeutrino("nuB_L", nuB_L, Engine::RHC);
		theEngine->BindNeutrino("nuB_R", nuB_R, Engine::RHC);

		outName.insert(outName.find(".root"), "_dirac_");
	}
	else //majorana
	{
		Neutrino nu_L(mass, Neutrino::Left  | Neutrino::Majorana );
		Neutrino nu_R(mass, Neutrino::Right | Neutrino::Majorana );

		nu_L.SetDecayChannel(channel);
		nu_R.SetDecayChannel(channel);

		theEngine->BindNeutrino("nu_L", nu_L, Engine::FHC);
		theEngine->BindNeutrino("nu_R", nu_R, Engine::RHC);

		outName.insert(outName.find(".root"), "_major");
	}

	std::stringstream ssm;
	ssm << std::setfill('0') << std::setw(4) << mass*1000;
	outName.insert(outName.find(".root"), ssm.str());

	theEngine->MakeFlux();		//create flux scaling by mass
	theEngine->ScaleToDetector(theBox);

	std::map<std::string, double> weights, ratios;	//obtain MC weights
	std::map<std::string, double>::iterator iw;
	double total, step = theEngine->RangeWidth();
	if (mix > 0)
		total = theEngine->MakeSampler(theBox, weights, mix, mix, mix);
	else
		total = theEngine->MakeSampler(theBox, weights, ue, um, ut);
	std::cout << "total number of events for " << outName << " is " << total << std::endl;

	for (iw = weights.begin(); iw != weights.end(); ++iw)
	{
		std::cout << "weights " << iw->first << " " << iw->second << std::endl;
		ratios[iw->first] = iw->second / total;
		iw->second /= nevent * step;		//MC weights are normalised
	}

	//open output file and load tree
	TFile *outFile = new TFile(outName.c_str(), "RECREATE");

	double True, W, R;
	bool H, P;
	int PdgA, PdgB;
	double E_A, P_A, T_A, TheA, PhiA, M_A, In_A, Out_A;
	double e_A, p_A, t_A, theA, phiA;
	double E_B, P_B, T_B, TheB, PhiB, M_B, In_B, Out_B;
	double e_B, p_B, t_B, theB, phiB;
	double E_0, P_0, T_0, The0, Phi0, M_0;
	double Angle;

	//CREATING TREE
	//
	//
	TTree *data = new TTree("Data", "Particle tracks");

	//true information
	data->Branch("True",   &True,   "fTrue/D");
	data->Branch("W", &W, "fW/D");
	data->Branch("R", &R, "fR/D");
	data->Branch("H", &H, "bH/O");
	data->Branch("P", &P, "bP/O");

	//particle A
	data->Branch("PdgA", &PdgA, "fPdgA/I");
	data->Branch("E_A", &E_A, "fEnergyA/D");
	data->Branch("P_A", &P_A, "fMomentA/D");
	data->Branch("T_A", &T_A, "fTransvA/D");
	data->Branch("TheA", &TheA, "fThetaA/D");
	data->Branch("PhiA", &PhiA, "fPhiA/D");
	//data->Branch("M_A", &M_A, "fMassA/D");
	data->Branch("In_A", &In_A, "fLengthA/D");
	data->Branch("Out_A", &Out_A, "fLengthoA/D");
	//data->Branch("e_A", &e_A, "fenergyA/D");
	//data->Branch("p_A", &p_A, "fmomentA/D");
	//data->Branch("t_A", &t_A, "ftransvA/D");
	//data->Branch("theA", &theA, "fthetaA/D");
	//data->Branch("phiA", &phiA, "fphiA/D");

	//Particle B
	data->Branch("PdgB", &PdgB, "fPdgB/I");		//true PDG code
	data->Branch("E_B", &E_B, "fEnergyB/D");
	data->Branch("P_B", &P_B, "fMomentB/D");
	data->Branch("T_B", &T_B, "fTransvB/D");
	data->Branch("TheB", &TheB, "fThetaB/D");
	data->Branch("PhiB", &PhiB, "fPhiB/D");
	//data->Branch("M_B", &M_B, "fMassB/D");
	data->Branch("In_B", &In_B, "fLengthB/D");
	data->Branch("Out_B", &Out_B, "fLengthoB/D");
	//data->Branch("e_B", &e_B, "fenergyB/D");
	//data->Branch("p_B", &p_B, "fmomentB/D");
	//data->Branch("t_B", &t_B, "ftransvB/D");
	//data->Branch("theB", &theB, "fthetaB/D");
	//data->Branch("phiB", &phiB, "fphiB/D");

	data->Branch("Angle", &Angle, "fAngle/D");

	//Particle 0 = A + B
	data->Branch("E_0", &E_0, "fEnergy0/D");
	data->Branch("P_0", &P_0, "fMoment0/D");
	data->Branch("T_0", &T_0, "fTransv0/D");
	data->Branch("The0", &The0, "fTh0ta0/D");
	data->Branch("Phi0", &Phi0, "fPhi0/D");
	data->Branch("M_0", &M_0, "fMass0/D");
	//
	//
	//END TREE

	TRandom3 *ran = new TRandom3(0);
	
	int ND = 0, ID;
	for (ID = 0; ID < nevent; ++ND)
	{
		//sampling 2/4 times at the same times
		//(over all the neutrinos loaded)
		std::map<std::string, double> energy, intensity;
		theEngine->SampleEnergy(energy, intensity);

		for (iw = energy.begin(); iw != energy.end(); ++iw)
		{
			if (iw->second < 0)
				continue;

			True = iw->second;	//true neutrino energy from sampling
			W = weights[iw->first];
			R = ratios[iw->first];

			//using uuid to determing helicity and particle/antiparticle nature
			if (iw->first.find("L") != std::string::npos)
				H = false;
			else if (iw->first.find("R") != std::string::npos)
				H = true;
			if (iw->first.find("0") != std::string::npos)
				P = true;
			else if (iw->first.find("B") != std::string::npos)
				P = false;

			Neutrino nu = theEngine->GetNeutrino(iw->first);

			nu.SetEnergy(True);
			theBox->Vertex(nu);
			theBox->Focus(nu);
			std::vector<Particle> product = nu.DecayPS(), particle;
			std::vector<Particle>::iterator ip;

			for (ip = product.begin(); ip != product.end(); ++ip)
			{
				int pdg = std::abs(ip->Pdg());
				ip->SetPosition(nu.X(), nu.Y(), nu.Z()); //set vertex

				//skip invisibles neutrinos
				if (pdg == 12 || pdg == 14 || pdg == 16)
					continue;
				else if (pdg == 111)		//pi0, must decay rn
				{
					Particle pA, pB;
					theBox->Pi0Decay(*ip, pA, pB);
					particle.push_back(pA);
					particle.push_back(pB);
				}
				else
					particle.push_back(*ip);
			}

			if (particle.size() == 2 &&
				theBox->Reconstruct(particle.at(0)) &&
				theBox->Reconstruct(particle.at(1)))
			{
				PdgA = particle.at(0).Pdg();
				E_A = particle.at(0).Energy();
				P_A = particle.at(0).Momentum();
				T_A = particle.at(0).Transverse();
				TheA = particle.at(0).Theta();
				PhiA = particle.at(0).Phi();
				//M_A = particle.at(0).Mass();
				In_A = particle.at(0).TrackIn();
				Out_A = particle.at(0).TrackOut();

				PdgB = particle.at(0).Pdg();
				E_B = particle.at(1).Energy();
				P_B = particle.at(1).Momentum();
				T_B = particle.at(1).Transverse();
				TheB = particle.at(1).Theta();
				PhiB = particle.at(1).Phi();
				//M_B = particle.at(1).Mass();
				In_B = particle.at(1).TrackIn();
				Out_B = particle.at(1).TrackOut();

				Angle = particle.at(0).Direction().Angle(particle.at(1).Direction());

				TLorentzVector Reco = particle.at(0).FourVector() + particle.at(1).FourVector();

				E_0 = Reco.E();
				P_0 = Reco.P();
				T_0 = Reco.Pt();
				The0 = Reco.Theta();
				Phi0 = Reco.Phi();
				M_0 = Reco.M();

				//TLorentzVector vA = particle.at(0).FourVector();
				//TLorentzVector vB = particle.at(1).FourVector();
				//TVector3 bst = - Reco.BoostVector();
				//vA.Boost(bst);
				//vB.Boost(bst);

				//e_A  = vA.E();
				//p_A  = vA.P();
				//t_A  = vA.Pt();
				//theA = vA.Theta();
				//phiA = vA.Phi();

				//e_B  = vB.E();
				//p_B  = vB.P();
				//t_B  = vB.Pt();
				//theB = vB.Theta();
				//phiB = vB.Phi();

				data->Fill();

				++ID;
			}

			if (ID % 10000 == 0)
			{
				outFile->cd();
				data->Write("", TObject::kWriteDelete);
			}
		}
	}

	int size = weights.size() / (dirac ? 1 : 2);
	std::cout << "Above detection thresholds there are " << (ID * 100.0 / size)/ND << "\% of simulated particles" << std::endl;

	outFile->cd();
	data->Write("", TObject::kWriteDelete);
	outFile->Close();

	return 0;
}

void Usage(char* argv0)
{
	std::cout << "Description" << std::endl;
	std::cout << "Usage : " << std::endl;
	std::cout << argv0 << " [OPTIONS]" << std::endl;
	std::cout <<"\n  -d,  --detconfig" << std::endl;
	std::cout << "\t\tDetector configuration file" << std::endl;
	std::cout <<"\n  -f,  --fluxconfig" << std::endl;
	std::cout << "\t\tFlux configuration file" << std::endl;
	std::cout <<"\n  -m,  --mass" << std::endl;
	std::cout << "\t\tNeutrino mass" << std::endl;
	std::cout <<"\n  -c,  --channel" << std::endl;
	std::cout << "\t\tDecay channel, defaul ALL" << std::endl;
	std::cout <<"\n  -n,  --number" << std::endl;
	std::cout << "\t\tNumber of events" << std::endl;
	std::cout <<"\n  -o,  --output" << std::endl;
	std::cout << "\t\tOutput root file" << std::endl;
	std::cout <<"\n  -L,  --left" << std::endl;
	std::cout << "\t\tNeutrino polarisation left, default unpolarised" << std::endl;
	std::cout <<"\n  -R,  --right" << std::endl;
	std::cout << "\t\tNeutrino polarisation right, default unpolarised" << std::endl;
	std::cout <<"\n  -h,  --help" << std::endl;
	std::cout << "\t\tPrint this message and exit" << std::endl;
}
