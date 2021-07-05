#include "montecarlo/GENIEback.h"

namespace GENIEback {

	// more efficient to check multiple channels at once
	std::map<Decay::Channel, std::shared_ptr<hnl> > GenerateBackground(const Tracker &box,
						std::vector<Decay::Channel> channels,
						std::string file, double weight, //bool chargeID,
						bool kVerbose)
	{
		// output of type hnl for each channel studied
		std::map<Decay::Channel, std::shared_ptr<hnl> > hnls;
		for (auto & c : channels)
			hnls[std::move(c)] = std::shared_ptr<hnl>(new hnl);
		channels.clear();

		TFile *inb = new TFile(file.c_str());
		if (inb->IsZombie())
			throw std::invalid_argument("GENIEback: file \"" + file + "\" does not exist");
		if (!inb->Get("gst"))
			std::cerr << "GENIEback: no gst tree found in \"" << file << "\"\n";

		// input
		gst* genie = new gst(static_cast<TTree*> (inb->Get("gst")));

		weight /= genie->GetEntries();

		int wid = std::log10(genie->GetEntries()-1) + 1;
		std::string spacer(wid+9, ' ');

		for (long int id = 0; id < genie->GetEntries(); ++id)
		{
			genie->GetEntry(id);
			// repeat analysis for all modules, to populate both detectors
			for (const auto &mod : box.Modules()) {
	
				// module ratio for MC weighting
				double rati = box.Weight(mod) / box.Weight();

				//neutrino probe
				Particle nu(genie->neu, genie->Ev, genie->pxv, genie->pyv, genie->pzv);

				// define a rotation to align the final state particles
				Tracker::Event nu_evt = box.GenerateEvent(mod, std::move(nu));

				TRotation rz;
				rz.SetZAxis(nu_evt.first.Vect());

				// pdg of outgoing lepton is determined by interaction type
				// if CC then pdg is ±1 the neutrino one
				// if NC then outgoing particle is same as incoming neutrino
				int pdgl = genie->cc
					? (genie->neu > 0 ? genie->neu - 1 : genie->neu + 1)
					: genie->neu;

				if (kVerbose) {
					std::cout << "Event " << std::setw(wid) << id << " : "
						<< "probe (" << nu_evt.first.E() << ") in "
						<< mod << "\t";
					if (genie->cc)
						std::cout << "charged current,\t";
					else if (genie->nc)
						std::cout << "neutral current,\t";
					if (genie->qel)
						std::cout << "quasi elastic";
					else if (genie->mec)
						std::cout << "meson exchange";
					else if (genie->res)
						std::cout << "resonance";
					else if (genie->coh)
						std::cout << "coherent";
					else if (genie->dis)
						std::cout << "deep inelastic";
					std::cout << "\n" << spacer << "final state particles: "
						<< 1 + genie->nf << "\n";
					//std::cout << spacer << "\t" << pdgl;
					//for (int i = 0; i < genie->nf; ++i)
					//	std::cout << ", " << genie->pdgf[i];
					//std::cout << "\n";
				}

				std::vector<Tracker::Event> events;
				events.reserve(genie->nf+1);

				//outgoing lepton from neutrino
				//if NC event, the pdg is the same of the probe
				//if CC eventm the pdg is the respective charged lepton
				Particle lep(pdgl, genie->El, genie->pxl, genie->pyl, genie->pzl);
				lep.Transform(rz);		// align with neutrino probe
				Tracker::Event lep_evt(std::move(lep), nu_evt.second);

				if (box.Reconstruct(lep_evt))
					events.push_back(std::move(lep_evt));

				// other particles
				bool background = true;
				for (int i = 0; i < genie->nf; ++i) {
					Particle part(genie->pdgf[i], genie->Ef[i],
						genie->pxf[i], genie->pyf[i], genie->pzf[i]);
					// rotate along neutrino direction
					part.Transform(rz);
					Tracker::Event part_evt(std::move(part), nu_evt.second);
					//special treatments for pi0
					if (std::abs(part_evt.first.Pdg() == 111))
					{	//almost 100% into 2photons
						auto decay_res = box.Pi0Decay(std::move(part_evt));
						Tracker::Event gA_evt = std::move(decay_res[0]);
						Tracker::Event gB_evt = std::move(decay_res[1]);

						if (box.Reconstruct(gA_evt))
							events.push_back(std::move(gA_evt));
						if (box.Reconstruct(gB_evt))
							events.push_back(std::move(gB_evt));
					}
					// do not process nuclei
					//else if (std::abs(part_evt.first.Pdg()) < 1e9)
					else if (box.Reconstruct(part_evt)) {
						if (!Process::IsDetectable(part_evt.first)) {
							background = false;
							break;
						}
						events.push_back(std::move(part_evt));
					}

					// events contains muons, electron, pions, protons,
					// kaons and other strange and charmed kaons

					// check if last particle is a potential background event
					//if (events.size() && !Process::IsDetectable(events.back().first)) {
						//background = false;
						//break;
					//}
				}


				// no good for potential background candidate
				if (!background || events.size() < 2) {
					if (kVerbose)
						std::cout << "\n";
					continue;
				}

				if (kVerbose) {
					std::string sep = spacer + "recos -> ";
					for (const auto &t : events) {
						std::cout << sep << t.first.RealPdg();
						//<< t.first.EKin() << ")";
						sep = ", ";
					}
					std::cout << "\n";
				}

				//make some simple misidentification of events
				events = box.MisIdentify(std::move(events));

				if (kVerbose) {
					//std::cout << spacer << "particles detected: "
					//<< events.size() << "\n";
					std::string sep = spacer + "detec -> ";
					for (const auto &t : events) {
						std::cout << sep << t.first.Pdg();
						sep = ", ";
					}
					std::cout << "\n";
				}

				//go through the particles and count if there is
				//the right number of particle we expect for process
				for (auto & dc : hnls) {
					auto id = Process::Identify(dc.first, events);
					if (id == Process::Match::no_id)
						continue;

					dc.second->Fill(weight * genie->wght * rati,
							(id == Process::Match::charge_id),
							nu_evt, events[0], events[1]);

					if (kVerbose) {
						std::cout << "\nTRIGGER "
							<< Decay::toString(dc.first)
							<< " is " << id
							<< "\t" << dc.second->GetEntries() << "\n";
						std::chrono::milliseconds stop(500);
						std::this_thread::sleep_for(stop);
					}
				}
			}	// loop on detector modules end
		}	// loop on entires end

		inb->Close();
		return hnls;
	}
}
