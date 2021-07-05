#ifndef GENIEBACK_H
#define GENIEBACK_H

#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <chrono>
#include <thread>

#include "TTree.h"
#include "TFile.h"

#include "tools/RNG.h"

#include "physics/Const.h"
#include "physics/Decays.h"
#include "physics/Particle.h"
#include "detector/Track.h"
#include "detector/Tracker.h"

#include "montecarlo/gst.h"
#include "montecarlo/hnl.h"
#include "montecarlo/Process.h"

namespace GENIEback
{
	std::map<Decay::Channel, std::shared_ptr<hnl> >
		GenerateBackground(const Tracker &box, std::vector<Decay::Channel> chan,
				std::string file, double weight = 1., //bool chargeID = false,
				bool kVerbose = false);
	Process::Match Identify(const std::vector<Tracker::Event> &events, Decay::Channel chan);
			//bool chargeID = false);
};

#endif
