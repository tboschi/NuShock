#ifndef PROCESS_H
#define PROCESS_H

#include "physics/Channels.h"
#include "tools/Const.h"

namespace Process {
	std::vector<double> Mass(Channel::Name chan);
	std::vector<int> Pdg(Channel::Name chan);
	std::vector<std::pair<double, int> > MassPdg(Channel::Name chan);
}

#endif
