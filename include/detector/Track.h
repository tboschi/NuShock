/*
 * Track class, collection of information of particle track
 * such as starting point, length, sagitta
 *
 * Author: Tommaso Boschi
 */

#ifndef TRACK_H
#define TRACK_H

#include <iostream>
#include <map>
#include <numeric>

#include "TVector3.h"

//#include "Tools.h"

class Track : public TVector3
{
	public:
		// pass starting point
		Track(double x = 0., double y = 0., double z = 0.);
		Track(const TVector3 &v);

		double Length() const;
		double Length(const std::string &mod) const;
		double EnergyDeposited() const;
		double EnergyDeposited(const std::string &mod) const;
		double Importance() const;
		double Importance(const std::string &mod) const;
		bool IsShower() const;

		
		void SetLength(const std::string &mod, double track = 0.);
		void SetShower(bool shower = true);
		void SetEnergyDeposited(const std::string &mod, double energy);

		void SetPosition(double x, double y, double z);
		void SetPosition(const TVector3 &v);

	protected:
		// first is length, second is energy deposited in GeV
		std::map<std::string, double> _tracks;
		std::map<std::string, double> _calors;
		bool kShower;

		friend std::ostream & operator<<(std::ostream &os, const Track &t);
};

#endif
