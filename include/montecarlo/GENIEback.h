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
#include "physics/Channels.h"
#include "physics/Particle.h"
#include "detector/Track.h"
#include "detector/Tracker.h"

#include "montecarlo/gst.h"
#include "montecarlo/hnl.h"

namespace GENIEback
{
	std::map<Channel::Name, std::shared_ptr<hnl> > GenerateBackground(const Tracker &box,
						std::vector<Channel::Name> chan, std::string file,
						double weight = 1., bool chargeID = false, bool kVerbose = false);
	bool Identify(const std::vector<Tracker::Event> &events, Channel::Name chan, bool chargeID = false);
	bool IsPotentialBackground(const Particle &p);
};

#endif
