/*
 * Author: Tommaso Boschi
 */

#ifndef EXCLUSIONLNV_H
#define EXCLUSIONLNV_H

#include <iostream>
#include <vector>
#include <cmath>
#include <list>

#include "Tools.h"
#include "Flux.h"
#include "Physics.h"
#include "Detector.h"

class ExclusionLNV
{
	public:
		ExclusionLNV(Engine* TE, Detector *TB, std::vector<char> &vF);
		~ExclusionLNV();
		double Bisect(double S, double E, double &nFunction);
		bool FindInterval(double S, double &M, double E);
		double IsZero(std::list<double> &ll);
		void Split(std::list<double> &ll);
		double Function(double lu2);
		double ZeroCross(double lu2);
		void SetMix(double lu2);

		double fcLUT(double bb);

	private:
		Engine *TheEngine;
		Engine::Current Horn;
		Detector *TheBox;
		std::vector<char> vFlag;
		std::vector<double> vss;
		double ithr;
};

#endif
