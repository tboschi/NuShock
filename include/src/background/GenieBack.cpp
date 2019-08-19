#include "GenieBack.h"

GenieBack::GenieBack(std::string backFile, std::string outFile, bool verb) :
	kVerbose(verb)
{
	InitInTree(backFile);
	InitOutTree(outFile);
	checkPt = 0;
}

GenieBack::~GenieBack()
{
	inBack->Close();
	data->Write("", TObject::kWriteDelete);
	outBack->Close();
}

int GenieBack::Charge(int pdg)
{
	if (chargeID)
		return pdg;
	else
		return std::abs(pdg);
}

void GenieBack::InitInTree(std::string backFile)
{
 	inBack = new TFile(backFile.c_str(), "READ");
	genie = new gst(static_cast<TTree*> (inBack->Get("gst")));
 	if(!genie)
		std::cerr << "No gst tree found in genie file" << std::endl;

	genie->fChain->SetBranchStatus("*",    0);
	genie->fChain->SetBranchStatus("cc",   1);	//Is CC?
	genie->fChain->SetBranchStatus("nc",   1);	//Is NC?
	genie->fChain->SetBranchStatus("qel",  1);	//Is quasi elastic?
	genie->fChain->SetBranchStatus("res",  1);	//Is resonance?
	genie->fChain->SetBranchStatus("dis",  1);	//Is DIS?
	genie->fChain->SetBranchStatus("coh",  1);	//Is coherent?
	genie->fChain->SetBranchStatus("neu",  1);	//neutrino PDG
	genie->fChain->SetBranchStatus("Ev",   1);	//neutrino E
	genie->fChain->SetBranchStatus("pxv",  1);	//neutrino px
	genie->fChain->SetBranchStatus("pyv",  1);	//neutrino py
	genie->fChain->SetBranchStatus("pzv",  1);	//neutrino pz
	genie->fChain->SetBranchStatus("El",   1);	//lepton E
	genie->fChain->SetBranchStatus("pxl",  1);	//lepton px
	genie->fChain->SetBranchStatus("pyl",  1);	//lepton py
	genie->fChain->SetBranchStatus("pzl",  1);	//lepton pz
	genie->fChain->SetBranchStatus("nf",   1);	//number of final states
	//vectors of size nf
	genie->fChain->SetBranchStatus("pdgf", 1);	//final pdg
	genie->fChain->SetBranchStatus("Ef",   1);	//final E
	genie->fChain->SetBranchStatus("pxf",  1);	//final px
	genie->fChain->SetBranchStatus("pyf",  1);	//final py
	genie->fChain->SetBranchStatus("pzf",  1);	//final pz
}

//this can be highly customised!
//
void GenieBack::InitOutTree(std::string outFile)
{
	outBack = new TFile(outFile.c_str(), "RECREATE");
	data = new TTree("Data", "All data from GENIE");	//Tree with candidate events

	//hist = new TH1D("Enu", "true energy neutrino", 100, 0, 20);

	/*
	data->Branch("ID", &ID, "iID/I");	//global event ID

	data->Branch("np", &np, "inp/I");	//number of particle

	data->Branch("E",   energy, "energy[np]/D");
	data->Branch("P",   moment, "moment[np]/D");
	data->Branch("T",   transv, "transv[np]/D");
	data->Branch("The", theta,  "theta[np]/D");
	data->Branch("Phi", phi,    "theta[np]/D");
	data->Branch("M",   mass,   "mass[np]/D");
	data->Branch("In",  lenOut, "lenOut[np]/D");
	data->Branch("Out", lenIn,  "lenIn[np]/D");

	data->Branch("nr", &nr, "inr/I");	//number of couples

	data->Branch("E0",   r_energy, "energy0[nr]/D");
	data->Branch("P0",   r_moment, "moment0[nr]/D");
	data->Branch("T0",   r_transv, "transv0[nr]/D");
	data->Branch("The0", r_theta,  "theta0[nr]/D");
	data->Branch("Phi0", r_phi,    "phi0[nr]/D");
	data->Branch("M0",   r_mass,   "mass0[nr]/D");
	data->Branch("Sep",  r_angle,  "angle0[nr]/D");	//separating angle
	*/

	//data->Branch("True", &True, "fTrue/D");
	data->Branch("PdgA", &PdgA, "fPdgA/I");		//true PDG code
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

	//and other variables
}

//Load tree
void GenieBack::LoadTree(const std::vector<Particle> &particle,
			 const std::vector<Particle> &original)
{
	/*
	np = particle.size();
	for (int i = 0; i < particle.size(); ++i)
	{
		energy[i] = particle.at(i).Energy();
		moment[i] = particle.at(i).Momentum();
		transv[i] = particle.at(i).Transverse();
		theta [i] = particle.at(i).Theta();
		phi   [i] = particle.at(i).Phi();
		mass  [i] = particle.at(i).Mass();
		lenIn [i] = particle.at(i).TrackIn();
		lenOut[i] = particle.at(i).TrackOut();
	}

	nr = particle.size();
	nr *= (nr - 1) / 2.0;
	for (int i = 0; i < particle.size(); ++i)
	{
		for (int j = i+1; j < particle.size(); ++j)
		{
			Particle reco(0, particle[i].Energy()+particle[j].Energy(),
					 particle[i].MomentumX()+particle[j].MomentumX(),
					 particle[i].MomentumY()+particle[j].MomentumY(),
					 particle[i].MomentumZ()+particle[j].MomentumZ(),
					 0, 0, 0);
			r_energy[i] = reco.Energy();
			r_moment[i] = reco.Momentum();
			r_transv[i] = reco.Transverse();
			r_theta [i] = reco.Theta();
			r_phi   [i] = reco.Phi();
			r_mass  [i] = reco.Mass();
			r_angle [i] = particle[i].Direction().Angle(particle[j].Direction());
		}
	}
	*/

	if (particle.size() == 2)
	{
		PdgA = original.at(0).Pdg();
		E_A = particle.at(0).Energy();
		P_A = particle.at(0).Momentum();
		T_A = particle.at(0).Transverse();
		TheA = particle.at(0).Theta();
		PhiA = particle.at(0).Phi();
		M_A = particle.at(0).Mass();
		In_A = particle.at(0).TrackIn();
		Out_A = particle.at(0).TrackOut();

		PdgB = original.at(1).Pdg();
		E_B = particle.at(1).Energy();
		P_B = particle.at(1).Momentum();
		T_B = particle.at(1).Transverse();
		TheB = particle.at(1).Theta();
		PhiB = particle.at(1).Phi();
		M_B = particle.at(1).Mass();
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
	}
}

TTree *GenieBack::FindBackground(Tracker *theTrack,
			const std::map<int, int> &process, int save)
{
	//useful if saving externally every other n = #save entries
	//the class wiremember the last event processed and will restart from there

	int batch;
	//if save negative, restart from 0, i.e. process/save everything
	if (save < 0)
	{
		checkPt = 0;
		batch = genie->fChain->GetEntries();
	}
	else
		batch = genie->fChain->GetEntries() / save;

	std::cout << "start at " << checkPt << std::endl;
	//for (ID = checkPt; ID < checkPt + batch && ID < genie->fChain->GetEntries(); ++ID)
	for (ID = 0; ID < batch; ++ID)
	//for (ID = 0; ID < genie->fChain->GetEntries(); ++ID)
	{
		bool hadronicActivity = false;
		genie->GetEntry(ID);	//get event from ID
		//hist->Fill(genie->Ev);

		//neutrino probe
		Particle nu(genie->neu, genie->pxv, genie->pyv, genie->pzv, genie->Ev);
		theTrack->Vertex(nu);
		theTrack->Focus(nu);

		std::vector<Particle> particle;

		//outgoing lepton from neutrino
		//if NC event, the pdg is the same of the probe
		//if CC eventm the pdg is the respective charged lepton
		int pdgl = genie->cc ? (genie->neu > 0 ? genie->neu - 1 : genie->neu + 1) : genie->neu;
		Particle p(pdgl, genie->pxl, genie->pyl, genie->pzl, genie->El, 
				 nu.X(), nu.Y(), nu.Z());
		if (theTrack->Reconstruct(p))
			particle.push_back(p);

		if (kVerbose)
		{
			std::cout << "Event " << ID << " : ";
			if (genie->cc)
				std::cout << "charged current,\t";
			else if (genie->nc)
				std::cout << "neutral current,\t";
			if (genie->qel)
				std::cout << "quasi elastic\n";
			else if (genie->res)
				std::cout << "resonance\n";
			else if (genie->coh)
				std::cout << "coherent\n";
			else if (genie->dis)
				std::cout << "deep inelastic\n";
			std::cout << "\tfinal state particles: " << particle.size() + genie->nf << std::endl;
			if (particle.size())
				std::cout << "\t" << pdgl << " (" << p.EnergyKin() << "), ";
		}

		for (int i = 0; i < genie->nf; ++i)
		{
			p = Particle(genie->pdgf[i],
				     genie->pxf[i], genie->pyf[i], genie->pzf[i],
				     genie->Ef[i]);
			if (kVerbose)
				std::cout << "\t" << p.Pdg() << " (" << p.EnergyKin() << "), ";

			/*
			if (std::abs(p.Pdg()) == 130 ||
					(std::abs(p.Pdg()) > 300 && std::abs(p.Pdg()) < 400))
				p.SetPdg(300);		//Kaons
			else if (std::abs(p.Pdg()) > 400 && std::abs(p.Pdg()) < 500)
				p.SetPdg(400);		//Charms
			else if (p.std::abs() != 0 && std::abs(p.Pdg()) > 1000)	
				p.SetPdg(1000);		//Hadrons
			*/

			if (std::abs(p.Pdg()) == 111)	//special treatments for pi0
			{				//almost 100% into 2photons
				Particle pA, pB;
				theTrack->Pi0Decay(p, pA, pB);

				if (theTrack->Reconstruct(pA))
					particle.push_back(pA);
				if (theTrack->Reconstruct(pB))
					particle.push_back(pB);
			}
			else if (std::abs(p.Pdg()) < 1e9)	//no nucleus
				//if true particle can be detected I will add it to the vector
				if (theTrack->Reconstruct(p))
					particle.push_back(p);
				//sould contains muons, electron, pions, protons, kaons and other strange and charmed kaons

			if (particle.size() && IsHadron(particle.back()))
			{
				hadronicActivity = true;
				break;
			}

		}
		std::cout << std::endl;

		if (!hadronicActivity && particle.size())
		{
			//make some simple misidentification of events
			std::vector<Particle> original = particle;
			MisIdentify(particle, theTrack);

			if (kVerbose)
			{
				std::cout << std::endl;
				std::cout << "\tparticles detected: " << particle.size() << std::endl;
				std::cout << "\toriginal: ";
				for (int i = 0; i < original.size(); ++i)
					std::cout << "\t" << original[i].Pdg() << ",";
				std::cout << std::endl;
				std::cout << "\tmis ID:   ";
				for (int i = 0; i < particle.size(); ++i)
					std::cout << "\t" << particle[i].Pdg() << ",";
				std::cout << std::endl;
			}

			//go through the particle and count if there is
			//the right number of particle we expect for process
			if (Identify(particle, process))
			{
				LoadTree(particle, original);
				if (kVerbose)
					std::cout << "TRIGGER!\t" << data->GetEntries() << std::endl;
			}
		}

		if (ID % batch == 0)
		{
			std::cout << "saving at " << ID << std::endl;
			data->Write("", TObject::kWriteDelete);
		}
	}

	data->Write("", TObject::kWriteDelete);
	//hist->Write();

	checkPt = ID;	//ID should be +1 
	return data;
}

/***** Identification of single particle and counting *****/
bool GenieBack::MisIdentify(std::vector<Particle> &particle, Tracker *theTrack)
{
	std::vector<Particle> electrons, photons;
	std::vector<Particle>::iterator ip, ie;

	for (ip = particle.begin(); ip != particle.end(); )
	{
		if (!theTrack->IsDetectable(*ip))
		{
			std::cout << "erasin"<< std::endl;
			ip = particle.erase(ip);
		}
		else
		{
			bool elec = true;
			switch (std::abs(ip->Pdg()))
			{
				case 13:
					break;
				case 211:
					if (!ip->IsShower() &&	//no shower & long track
					     ip->TrackIn() > theTrack->Get("LengthMIP"))
					{	//it is a muon
						if (kVerbose)
							std::cout << "Pion MisID!" << std::endl;
						ip->SetPdg(13);
						ip->SetMass(Const::MMuon);
						--ip;			//must recheck (for thr for instance)
					}
					break;
				case 11:
					elec = true;
					for (ie = electrons.begin(); ie != electrons.end(); )
					{
						double angle = ip->Direction().Angle(ie->Direction());
						if (angle < theTrack->Get("PairAngle") / Const::Deg)
						{		//two track are too close, it can be a pair production
							TVector3 mom = ip->Direction() + ie->Direction();
							double dist = sqrt(pow(ip->TrackIn(), 2) + pow(ie->TrackIn(), 2)
									+ 2*ip->TrackIn()*ie->TrackIn()*cos(angle));
							TLorentzVector gammaReco(mom, ip->Energy() + ie->Energy());

							*ip = Particle(22, gammaReco);
							ip->SetMass(0);
							ip->SetTrackIn(dist);

							--ip;			//must recheck (for thr for instance)

							if (kVerbose)
								std::cout << "Electron MisID!" << std::endl;
							ie = electrons.erase(ie);
							elec = false;
							break;
						}
						else ++ie;
					}
					if (elec)
						electrons.push_back(*ip);
					break;
				case 22:
					if(ip->TrackIn() < theTrack->Get("ConversionEM"))	//short displacement > 2cm
					{	//conversion before 3cm 
						if (kVerbose)
							std::cout << "Photon MisID!" << std::endl;
						ip->SetPdg(11);
						ip->SetMass(Const::MElectron);

						--ip;			//must recheck (for thr for instance)
					}
					break;
			}

			++ip;
		}
	}
}

//verify that particle contains the particle sought for, in mProc
bool GenieBack::Identify(const std::vector<Particle> &particle,
			 const std::map<int, int> &process)
{
	std::map<int, int> mpart;
	std::map<int, int>::const_iterator im, cm;
	std::vector<Particle>::const_iterator ip;
	for (ip = particle.begin(); ip != particle.end(); ++ip)
		mpart[std::abs(ip->Pdg())]++;

	if (process.size() != mpart.size())
		return false;

	for (im = process.begin(), cm = mpart.begin(); im != process.end(); ++im, ++cm)
		if (!(im->first == cm->first && im->second == cm->second))
			return false;
	return true;
}

bool GenieBack::IsHadron(const Particle &p)
{
	switch (std::abs(p.Pdg()))
	{
		case 11:
		case 13:
		case 15:
		case 22:
		case 211:
		case 111:
			return false;
		default:
			return true;
	}
}
