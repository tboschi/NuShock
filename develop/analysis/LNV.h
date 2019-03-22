/*
 * Author: Tommaso Boschi
 */

#ifndef EXCLUSION_LNV_H
#define EXCLUSION_LNV_H

#include <iostream>
#include <vector>
#include <cmath>
#include <list>

#include "Tools.h"
#include "Flux.h"
#include "Physics.h"
#include "Detector.h"

class Exclusion_LNV
{
	public:
		Exclusion_LNV(Engine* TE, Engine::Current HornType,
			  Detector *TB,
			  std::vector<char> &vF, double Threshold);
		~Exclusion_LNV();
		double Bisect(double S, double E);
		bool FindInterval(double S, double &M, double E);
		double IsZero(std::list<double> &ll);
		void Split(std::list<double> &ll);
		double Function(double lu2);
		void SetMix(double lu2);
		void SetThr(double T);

	private:
		Engine *TheEngine;
		Engine::Current Horn;
		Detector *TheBox;
		std::vector<char> vFlag;
		double Thr;
};

#endif
