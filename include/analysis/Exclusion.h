/*
 * Author: Tommaso Boschi
 */

#ifndef EXCLUSION_H
#define EXCLUSION_H

#include <iostream>
#include <vector>
#include <cmath>

#include "tools.h"
#include "flux.h"
#include "physics.h"
#include "detector.h"

class Exclusion
{
	public:
		Exclusion(Engine* TE, Detector *TB, Engine::Current type, 
				bool ue = false, bool um = false, bool ut = false, double thr = 2.44, double mod = 1);

		double Bisect(double S, double E, double &n);
		bool FindInterval(double S, double &M, double E);

		double Zero(double lu2);
		double NumberEvents(double lu2);
		void SetThreshold(double T);

	private:
		Engine *engine;
		Detector *box;
		Engine::Current horn;

		bool UeFlag, UmFlag, UtFlag;
		double threshold, modifier;
};

#endif
