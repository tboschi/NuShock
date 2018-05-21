/*
 * Neutrino class
 * Author: Tommaso Boschi
 */

#ifndef NEUTRINO_H
#define NEUTRINO_H

#include <iostream>
#include <fstream>
#include <vector>
#include <cstring>

//ROOT include
#include "TLorentzVector.h"
#include "TGenPhaseSpace.h"

#include "Tools/Tools.h"
#include "Decay/ThreeBody.h"

enum Neut
{
	Particle     = 0;
	Antiparticle = 1;
	Left         = 2;
	Right        = 4;
};

class Neutrino:
{
	public:
		enum Neut
		{
			Particle     = 0;
			Antiparticle = 1;
			Left         = 2;
			Right        = 4;
		};

		Neutrino(double Mass, Neut);
	private:
};

#endif
