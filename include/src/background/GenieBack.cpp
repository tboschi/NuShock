#include "GenieBack.h"

GenieBack::GenieBack(std::string backFile)
{
	InitInTree(backFile);
	InitOutTree();
	checkPt = 0;
}

GenieBack::~GenieBack()
{
	inBack->Close();
}

void GenieBack::InitInTree(std::string backFile)
{
 	inBack = new TFile(backFile.c_str(), "OPEN");
	genie = new gst(static_cast<TTree*> (inBack->Get("gst")));
 	if(!genie)
		std::cerr << "No gst tree found in genie file" << std::endl;

	genie->fChain->SetBranchStatus("*",    0);
	genie->fChain->SetBranchStatus("neu",  1);	//neutrino PDG
	genie->fChain->SetBranchStatus("cc",   1);	//Is CC?
	genie->fChain->SetBranchStatus("nc",   1);	//Is NC?
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
void GenieBack::InitOutTree()
{
	data = new TTree("data", "All data from GENIE");	//Tree with candidate events

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

	//and other variables
}

//Load tree
void GenieBack::LoadTree(const std::vector<Particle> &particle)
{
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

	data->Fill();
}

TTree *GenieBack::FindBackground(Tracker *theTrack, std::map<int, int> &process, int save)
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

	for (ID = checkPt; ID < checkPt + batch && ID < genie->fChain->GetEntries(); ++ID)
	{
		genie->GetEntry(ID);	//get event from ID

		std::vector<Particle> particle;

		//outgoing lepton from neutrino
		//if NC event, the pdg is the same of the probe
		//if CC eventm the pdg is the respective charged lepton
		int pdgl = genie->cc ? (genie->neu > 0 ? genie->neu - 1 : genie->neu + 1) : genie->neu;
		Particle p(pdgl, genie->El, genie->pxl, genie->pyl, genie->pzl);

		for (int i = 0; i < genie->nf; ++i)
		{
			p = Particle(genie->pdgf[i], genie->Ef[i], genie->pxf[i], genie->pyf[i], genie->pzf[i]);

			if (abs(p.Pdg()) == 130 || (abs(p.Pdg()) > 300 && abs(p.Pdg()) < 400))
				p.SetPdg(300);		//Kaons
			else if (abs(p.Pdg()) > 400 && abs(p.Pdg()) < 500)
				p.SetPdg(400);		//Charms
			else if (p.Charge() != 0 && abs(p.Pdg()) > 1000)	
				p.SetPdg(1000);		//Hadrons

			if (abs(p.Pdg()) == 111)		//special treatments for pi0
			{					//almost 100% into 2photons
				Particle pA, pB;
				theTrack->Pi0Decay(p, pA, pB);

				if (theTrack->Reconstruct(pA))
					particle.push_back(pA);
				if (theTrack->Reconstruct(pB))
					particle.push_back(pB);
			}
			else if (abs(p.Pdg()) < 1e9)	//no nucleus
				//if true particle can be detected I will add it to the vector
				if (theTrack->Reconstruct(p))
					particle.push_back(p);
				//sould contains muons, electron, pions, protons, kaons and other strange and charmed kaons
		}

		//make some simple misidentification of events
		MisIdentify(particle, theTrack);

		//go through the particle and count if there is the number of particle we expect for process
		if (Identify(particle, process))
			LoadTree(particle);
	}
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
			ip = particle.erase(ip);
		else
		{
			bool elec = true;
			switch (abs(ip->Pdg()))
			{
				case 13:
					/*if (!MuonOrPion(*ip))
					{
						(*ip)->SetPdg(211);
						(*ip)->SetMass(M_Pion);
						--ip;			//must recheck (for thr for instance)
					}
					*/
					break;
				case 211:
					if (!ip->IsShower() &&						//no shower & long track
					     ip->TrackIn() > theTrack->Get("LengthMIP"))		//it is a muon
					{								
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
						{				//two track are too close, it can be a pair production
							TVector3 mom = ip->Direction() + ie->Direction();
							double dist = sqrt(pow(ip->TrackIn(), 2) + pow(ie->TrackIn(), 2)
									+ 2*ip->TrackIn()*ie->TrackIn()*cos(angle));
							TLorentzVector GammaReco(mom, ip->Energy() + ie->Energy());

							*ip = Particle(22, GammaReco);
							ip->SetMass(0);
							ip->SetTrackIn(dist);

							--ip;			//must recheck (for thr for instance)
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
					{		//conversion before 3cm 
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
bool GenieBack::Identify(const std::vector<Particle> &particle, const std::map<int, int> mProc)
{
	std::map<int, int> mCount;
	std::vector<Particle>::const_iterator ip;
	for (ip = particle.begin(); ip != particle.end(); ++ip)
	{
		mCount[ip->Pdg()]++;
	}

	return mCount.size() == mProc.size() && std::equal(mCount.begin(), mCount.end(), mProc.begin());
}
