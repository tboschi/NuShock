/*
 * Author: Tommaso Boschi
 */

#ifndef EXCLUSION_H
#define EXCLUSION_H

#include <iostream>
#include <vector>
#include <cmath>
#include <list>

#include "Tools.h"
#include "Flux.h"
#include "Physics.h"
#include "Detector.h"

class Exclusion
{
	public:
		Exclusion(Engine* TE, Detector *TB, bool Efficiency, std::vector<char> &vF, double Threshold);
		~Exclusion();
		double Bisect(double S, double E);
		bool FindInterval(double S, double &M, double E);
		double IsZero(std::list<double> &ll);
		void Split(std::list<double> &ll);
		double Function(double lu2);
		void SetMix(double lu2);

	private:
		Engine *TheEngine;
		Detector *TheBox;
		bool Eff;
		std::vector<char> vFlag;
		double Thr;
};

#endif
