/*
 * PhaseSpace generator for simulation of sterile neutrino decays (it is not needed for production)
 * It also returns correct 
 *
 * Author: Tommaso Boschi
 */

#ifndef DECAYRATES_H
#define DECAYRATES_H

#include <iostream>
#include <fstream>
#include <vector>
#include <cstring>

//ROOT include
#include "TLorentzVector.h"
#include "TGenPhaseSpace.h"

#include "Tools/Tools.h"
#include "Decay/ThreeBody.h"

#include "Physics/Neutrino.h"


class PhaseSpace
{
	public:
		PhaseSpace();

	private:
};

#endif
